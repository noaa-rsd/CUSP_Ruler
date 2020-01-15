import os
import json
from datetime import datetime
from pathlib import Path
from functools import partial
import numpy as np
import pandas as pd
import geopandas as gpd
import pyproj
import arcpy
from shapely.ops import transform, unary_union
from cartopy.geodesic import Geodesic


def set_env_vars(env_name):
    user_dir = os.path.expanduser('~')
    conda_dir = Path(user_dir).joinpath('AppData', 'Local', 
                                        'Continuum', 'anaconda3')
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
    arcpy.AddMessage('\n{}{}{}'.format(sep1, region.upper(), sep2))    


def define_args(config, arc_params):
    ref_shoreline = Path(str(arc_params[0].value))
    cusp_dir = Path(str(arc_params[1].value))
    out_dir = Path(str(arc_params[2].value))

    #reference_gdb_default = config['reference_gdb']    
    #cups_gdb_default = config['cusp_gdb']
    #out_dir_default = config['output_dir']
    
    simp = 0.002  # degrees
    buff = 500  # meters

    return ref_shoreline, cusp_dir, out_dir, simp, buff


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

    arcpy.AddMessage(splash)


def lat_lon_buffer(line, radius):
    _lon, _lat = line.centroid.coords[0]
    aeqd = pyproj.Proj(proj='aeqd', ellps='WGS84', datum='WGS84', 
                       lat_0=_lat, lon_0=_lon)

    wgs84 = pyproj.Proj(proj='latlong', ellps='WGS84')
    projected_line = transform(partial(pyproj.transform, wgs84, aeqd), line)
    transformed_buffer = transform(partial(pyproj.transform, aeqd, wgs84), 
                                   projected_line.buffer(radius))

    try:
        if not transformed_buffer.is_valid:
            _transformed_buffer = transformed_buffer.buffer(0.0001)
            assert _transformed_buffer.geom_type == 'Polygon'
            assert _transformed_buffer.is_valid
            transformed_buffer = _transformed_buffer
        return transformed_buffer
    except Exception as e:
        arcpy.AddMessage(r'GEOMETRY ERROR: polygon will be excluded from analysis')
        arcpy.AddMessage(transformed_buffer)


def buffer_lat_lon_multiline(multiline, radius):
    buffers = [lat_lon_buffer(line, radius) for line in multiline]
    return unary_union([b for b in buffers if b])


def distance(pt1, pt2):  # from https://pelson.github.io/2018/coast-path/
    result = np.array(geod.inverse(np.asanyarray(pt1), np.asanyarray(pt2)))
    return result[:, 0]


def calc_line_length(geom):  # from https://pelson.github.io/2018/coast-path/
    if hasattr(geom, 'geoms'):
        return sum(calc_line_length(subgeom) for subgeom in geom.geoms)
    else:
        points = np.array(geom.coords)
        return distance(points[:-1, :2], points[1:, :2]).sum()


def format_output_file(output_file):
    temp_file = output_file.replace('.txt', '_.txt')
    with open(temp_file, 'w') as new_f:
        with open(output_file, 'r') as old_f:
            for line in old_f:
                new_line = line.rstrip() + ';' + os.linesep
                arcpy.AddMessage(new_line)
                new_f.write(new_line)
    os.remove(output_file)
    os.rename(temp_file, output_file)


def load_default_paths(default_paths_path):
    with open(default_paths_path) as f:
        paths = json.load(f)
    return paths


def save_config(reference, cusp, out_dir, default_paths_path):
    default_paths = {'reference_gdb': str(list(reference.parents)[1]), 
                     'cusp_gdb': str(list(cusp.parents)[1]),
                     'output_dir': str(list(out_dir.parents)[0])}

    with open(default_paths_path, 'w') as f:
        json.dump(default_paths, f)


if __name__ == '__main__':
    tic = datetime.now()
    print_splash()

    cusp_ruler_dir = Path(os.path.dirname(os.path.realpath(__file__)))
    default_paths_path = cusp_ruler_dir / 'default_paths.txt'

    geod = Geodesic()

    default_paths = load_default_paths(default_paths_path)
    arc_params = arcpy.GetParameterInfo()
    reference, cusp_dir, out_dir, simp, buff = define_args(default_paths, arc_params)
    #save_config(reference, cusp, out_dir, default_paths_path)  TODO: account for multiple FYs
    env_name = 'cusp'
    set_env_vars(env_name)

    wgs84_epsg = {'init': 'epsg:4326'}
    
    arcpy.AddMessage('reading {}...'.format(reference))
    ref_gdb = str(reference.parent)
    ref_gdf = gpd.read_file(ref_gdb, layer=reference.name).to_crs(wgs84_epsg)

    cusps = {
        'FY11': '20110908contemporary_shoreline.gdb',
        'FY12': '20120814contemporary_shoreline.gdb',
        'FY13': '20130829contemporary_shoreline_1.gdb',
        'FY14': '20140930contemporary_shoreline.gdb',
        'FY15': '20150928bcontemporary_shoreline.gdb',
        'FY16': '20160927contemporary_shoreline.gdb',
        'FY17': '20170913_contemporary_shoreline.gdb',
        'FY18': '20180918_contemporary_shoreline.gdb',
        'FY19': '20190925contemporary_shoreline.gdb',
        }

    results = []
    rounding = {'StateLeng': 0, 'ClippedLeng': 0, 'Percentage': 2}
    types = {'StateLeng': 'int64', 'ClippedLeng': 'int64'}

    for fy, gdb in cusps.items():
        arcpy.AddMessage('\n' + ' ***** '.join([fy] * 5))

        layer = 'shoreline'
        cusp_gdb = cusp_dir / gdb
        arcpy.AddMessage('reading {}...'.format(cusp_gdb))
        cusp_gdf = gpd.read_file(str(cusp_gdb), layer=layer, crs=wgs84_epsg)

        if 'NOAA_Regio' in cusp_gdf.columns:
            gc_indices = cusp_gdf.SOURCE_ID.str[0:2] == 'GC'
            cusp_gdf = cusp_gdf[gc_indices]

            for region in ['Alaska']:
                region_id = ''.join([c.capitalize() for c in region.split(' ')])
                cusp_gpkg = out_dir / 'CUSP_{}.gpkg'.format(region_id)
                cusp_shp = out_dir / 'CUSP_{}.shp'.format(region_id)
                print_region_header(region, '-')
                arcpy.AddMessage('extracting CUSP data...')
                cusp = cusp_gdf[cusp_gdf.NOAA_Regio == region]

                arcpy.AddMessage('\nextracting reference shoreline...')
                ref = ref_gdf[ref_gdf.NOAA_REGIO == region].copy()
                geodesic_lengths = [calc_line_length(g) for g in ref['geometry']]
                ref['km_total'] = geodesic_lengths

                arcpy.AddMessage('\nsimplifying CUSP data...')
                cusp['geometry'] = cusp.geometry.simplify(tolerance=simp,
                                                          preserve_topology=False)
                layer = '{}_CUSP_{}_simplified'.format(fy, region_id)
                cusp['geometry'].to_file(cusp_gpkg, layer=layer, driver='GPKG')

                arcpy.AddMessage('\nbuffering simplified CUSP data...')
                cusp_simp_buff = buffer_lat_lon_multiline(cusp.geometry, buff)
                cusp_buffer_gdf = gpd.GeoDataFrame(geometry=[cusp_simp_buff]).explode()
                layer = '{}_CUSP_{}_buffered'.format(fy, region_id)
                cusp_buffer_gdf.to_file(cusp_gpkg, layer=layer, driver='GPKG')

                arcpy.AddMessage('\nclipping reference shoreline with CUSP buffers...')
                ref_clipped = ref.copy()
                ref_clipped.geometry = ref.intersection(cusp_simp_buff)
                geodesic_lengths = [calc_line_length(g) for g in ref_clipped['geometry']]
                ref_clipped['km_mapped'] = geodesic_lengths

                arcpy.AddMessage('\nsumming regional stats...')
                ref_lengths = ref.groupby('State')['km_total'].sum()
                ref_clipped_lengths = ref_clipped.groupby('State')['km_mapped'].sum()

                df = pd.DataFrame({
                    'StateLeng': ref_lengths / 1000,
                    'ClippedLeng': ref_clipped_lengths / 1000,
                    'Percentage': ref_clipped_lengths / ref_lengths,
                    'noaaRegion': [region] * ref_lengths.shape[0],
                    'FY': [fy] * ref_lengths.shape[0]
                    })

                df.index.names = ['stateName']
                col_order = ['FY', 'noaaRegion', 'ClippedLeng', 'StateLeng', 'Percentage']
                arcpy.AddMessage(df[col_order].round(rounding).astype(types))
                results.append(df[col_order])

        else:
            arcpy.AddMessage('WARNING: no NOAA_Regio attribute')

    results_df = pd.concat(results).round(rounding).astype(types)
    print_region_header('ALL PROCESSED REGIONS', '=')
    arcpy.AddMessage(results_df)

    output_file = '{}\CUSP_Progress.txt'.format(out_dir)
    results_df.to_csv(output_file, sep=',')
    format_output_file(output_file)
    arcpy.AddMessage('TOTAL TIME: {}'.format(datetime.now() - tic))
