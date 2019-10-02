import os
from datetime import datetime
from pathlib import Path
from functools import partial

import tkinter as tk
from tkinter import filedialog

import numpy as np
import pandas as pd

import geopandas as gpd
import pyproj
from shapely.geometry import MultiPolygon
from shapely.ops import transform, unary_union
from cartopy.geodesic import Geodesic


def set_env_vars(env_name):
    user_dir = os.path.expanduser('~')
    conda_dir = Path(user_dir).joinpath('AppData', 'Local', 'Continuum', 'anaconda3')
    env_dir = conda_dir / 'envs' / env_name
    share_dir = env_dir / 'Library' / 'share'
    script_path = conda_dir / 'Scripts'
    gdal_data_path = share_dir / 'gdal'
    proj_lib_path = share_dir

    if script_path.name not in os.environ["PATH"]:
        os.environ["PATH"] += os.pathsep + str(script_path)
    os.environ["GDAL_DATA"] = str(gdal_data_path)
    os.environ["PROJ_LIB"] = str(proj_lib_path)


def print_region_header(region, sep):
    width = 54
    sep1_num = int(width / 2 - int(len(region) / 2))
    sep2_num = int(width - len(region) - sep1_num)
    sep1 = sep * sep1_num
    sep2 = sep * sep2_num
    print('\n{}{}{}'.format(sep1, region.upper(), sep2))    


def define_args():
    reference_gpd = Path(filedialog.askdirectory(title='Select Reference Shoreline Geodatabase'))
    reference = reference_gpd / 'Whole_US'
    cusp_gdb = Path(filedialog.askdirectory(title='Select CUSP Geodatabase'))
    cusp = cusp_gdb / 'shoreline'
    out_dir = Path(filedialog.askdirectory(title='Select Output Directory'))
    simp = 0.002  # degrees
    buff = 500  # meters
    return reference, cusp, out_dir, simp, buff


def print_splash():
    splash = r"""    
    ==================================================================
                       NOAA Remote Sensing Division's
      _____ _    _  _____ _____    _____  _    _ _      ______ _____  
     / ____| |  | |/ ____|  __ \  |  __ \| |  | | |    |  ____|  __ \ 
    | |    | |  | | (___ | |__) | | |__) | |  | | |    | |__  | |__) |
    | |    | |  | |\___ \|  ___/  |  _  /| |  | | |    |  __| |  _  / 
    | |____| |__| |____) | |      | | \ \| |__| | |____| |____| | \ \ 
     \_____|\____/|_____/|_|      |_|  \_\\____/|______|______|_|  \_\
    ==================================================================
    """

    print(splash)


def lat_lon_buffer(line, radius):
    _lon, _lat = line.centroid.coords[0]
    aeqd = pyproj.Proj(proj='aeqd', ellps='WGS84', datum='WGS84', lat_0=_lat, lon_0=_lon)
    projected_line = transform(partial(pyproj.transform, wgs84_globe, aeqd), line)
    transformed_buffer = transform(partial(pyproj.transform, aeqd, wgs84_globe), 
                                   projected_line.buffer(radius))

    try:
        if not transformed_buffer.is_valid:
            _transformed_buffer = transformed_buffer.buffer(0.0001)
            assert _transformed_buffer.geom_type == 'Polygon'
            assert _transformed_buffer.is_valid
            transformed_buffer = _transformed_buffer
        return transformed_buffer
    except Exception as e:
        print(r"ERROR: POLYGON GEOMETRY ERROR (polygon will be excluded from analysis):")
        print(transformed_buffer)


def buffer_lat_lon_multiline(multiline, radius):
    print(' (buffering line segments...)')
    buffers = [lat_lon_buffer(line, radius) for line in multiline]
    print(' (creating union of buffered line segments...)')
    return unary_union([b for b in buffers if b])


def distance(pt1, pt2):  # from https://pelson.github.io/2018/coast-path/
    result = np.array(geod.inverse(np.asanyarray(pt1), np.asanyarray(pt2)))
    return result[:, 0]


def linestring_distance(geom):  # from https://pelson.github.io/2018/coast-path/
    if hasattr(geom, 'geoms'):
        return sum(linestring_distance(subgeom) for subgeom in geom.geoms)
    else:
        points = np.array(geom.coords)
        return distance(points[:-1, :2], points[1:, :2]).sum()


def format_output_file(output_file):
    temp_file = output_file.replace('.txt', '_.txt')
    with open(temp_file, 'w') as new_f:
        with open(output_file, 'r') as old_f:
            for line in old_f:
                new_line = line.rstrip() + ';' + os.linesep
                print(new_line)
                new_f.write(new_line)
    os.remove(output_file)
    os.rename(temp_file, output_file)


if __name__ == '__main__':
    tic = datetime.now()
    print_splash()

    root = tk.Tk()
    root.withdraw()

    geod = Geodesic()

    reference, cusp, out_dir, simp, buff = define_args()
    env_name = 'cusp'
    set_env_vars(env_name)

    wgs84_epsg = {'init': 'epsg:4326'}
    wgs84_globe = pyproj.Proj(proj='latlong', ellps='WGS84')
    
    print('reading {}...'.format(reference))
    ref_gdb = str(reference.parent)
    ref_gdf = gpd.read_file(ref_gdb, layer=reference.name).to_crs(wgs84_epsg)

    print('reading {}...'.format(cusp))
    cusp_gdb = str(cusp.parent)
    cusp_layer = cusp.name
    cusp_gdf = gpd.read_file(cusp_gdb, layer=cusp_layer, crs=wgs84_epsg)

    results = []
    rounding = {'StateLeng': 0, 'ClippedLeng': 0, 'Percentage': 2}
    types = {'StateLeng': 'int64', 'ClippedLeng': 'int64'}

    for region in cusp_gdf.NOAA_Regio.unique():
        region_id = ''.join([c.capitalize() for c in region.split(' ')])
        cusp_gpkg = out_dir / 'CUSP_{}.gpkg'.format(region_id)
        cusp_shp = out_dir / 'CUSP_{}.shp'.format(region_id)
        print_region_header(region, '-')
        print('extracting CUSP data...')
        cusp_region = cusp_gdf[cusp_gdf.NOAA_Regio == region]

        print('extracting reference shoreline...')
        ref_region = ref_gdf[ref_gdf.NOAA_REGIO == region].copy()
        print(' (measuring geodesic lengths...)')
        geodesic_lengths = [linestring_distance(g) for g in ref_region['geometry']]
        ref_region['km_total'] = geodesic_lengths

        print('simplifying CUSP data...')
        cusp_region_simp = cusp_region.copy()
        cusp_region_simp['geometry'] = cusp_region.geometry.simplify(tolerance=simp, preserve_topology=True)
        print(' (saving to geopackage...)')
        layer = 'CUSP_{}_simplified'.format(region_id)
        cusp_region_simp['geometry'].to_file(cusp_gpkg, layer=layer, driver='GPKG')

        print('buffering simplified CUSP data)...')
        cusp_region_simp_buff = buffer_lat_lon_multiline(cusp_region_simp.geometry, buff)
        cusp_buffer_gdf = gpd.GeoDataFrame(geometry=[cusp_region_simp_buff]).explode()
        print(' (saving to geopackage...)')
        layer = 'CUSP_{}_buffered'.format(region_id)
        cusp_buffer_gdf.to_file(cusp_gpkg, layer=layer, driver='GPKG')

        print('clipping reference shoreline with simplified CUSP buffers...')
        ref_region_clipped = ref_region.copy()
        ref_region_clipped.geometry = ref_region.intersection(cusp_region_simp_buff)
        print(' (measuring geodesic lengths...)')
        geodesic_lengths = [linestring_distance(g) for g in ref_region_clipped['geometry']]
        ref_region_clipped['km_mapped'] = geodesic_lengths

        print('summing regional stats...')
        ref_region_lengths = ref_region.groupby('State')['km_total'].sum()
        ref_region_clipped_lengths = ref_region_clipped.groupby('State')['km_mapped'].sum()

        df = pd.DataFrame({'StateLeng': ref_region_lengths / 1000,
                           'ClippedLeng': ref_region_clipped_lengths / 1000,
                           'Percentage': ref_region_clipped_lengths / ref_region_lengths,
                           'noaaRegion': [region] * ref_region_lengths.shape[0]})

        df.index.names = ['stateName']
        col_order = ['noaaRegion', 'ClippedLeng', 'StateLeng', 'Percentage']
        print(df[col_order].round(rounding).astype(types))
        results.append(df[col_order])

    results_df = pd.concat(results).round(rounding).astype(types)
    print_region_header('ALL PROCESSED REGIONS', '=')
    print(results_df)

    output_file = '{}\CUSP_Progress.txt'.format(out_dir)
    results_df.to_csv(output_file, sep=',')
    format_output_file(output_file)
    print('TOTAL TIME: {}'.format(datetime.now() - tic))