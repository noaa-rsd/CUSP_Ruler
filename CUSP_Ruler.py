import os
from datetime import datetime
from pathlib import Path
from functools import partial
import numpy as np
import pandas as pd
import geopandas as gpd
import pyproj
from shapely.geometry import MultiPolygon
from shapely.ops import transform, unary_union


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
    reference = Path(r'Z:\CUSP_progress\70K_shoreline_corrected\70K_Shoreline_4CUSP.gdb\Whole_US')
    cusp = Path(r'Z:\CUSP_progress\20190717_contemporary_shoreline.gdb\shoreline')
    out_dir = Path(r'Z:\CUSP_progress\Results')
    simp = 0.001  # degrees
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


def haversine(segments):
    segments = np.radians(segments)
    lon1, lat1 = segments[:, 0], segments[:, 1]
    lon2, lat2 = segments[:, 2], segments[:, 3]
    dlon = lon2 - lon1 
    dlat = lat2 - lat1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    c = 2 * np.arcsin(np.sqrt(a))
    r = 6371000  # meters (assumed Earth radius)
    return c * r


def calc_multiline_length(multiline):
    multiline_length = 0
    for line in multiline:
        a = line.coords[:]
        line_segments = np.hstack((a[0:-1], a[1:]))
        line_length = np.sum(haversine(line_segments))
        multiline_length += line_length
    return multiline_length


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


if __name__ == '__main__':
    tic = datetime.now()
    print_splash()

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
    rounding = {'km_total': 3, 'km_mapped': 3, 'pct_mapped': 3}
    types = {'km_total': 'int64', 'km_mapped': 'int64'}

    for region in cusp_gdf.NOAA_Regio.unique():
        region_id = ''.join([c.capitalize() for c in region.split(' ')])
        cusp_gpkg = out_dir / 'CUSP_{}.gpkg'.format(region_id)
        cusp_shp = out_dir / 'CUSP_{}.shp'.format(region_id)
        print_region_header(region, '-')
        print('extracting CUSP data...')
        cusp_region = cusp_gdf[cusp_gdf.NOAA_Regio == region]

        print('extracting reference shoreline...')
        ref_region = ref_gdf[ref_gdf.NOAA_REGIO == region].copy()
        print(' (measuring great-circle lengths...)')
        great_circle_lengths = [calc_multiline_length(g) for g in ref_region['geometry']]
        ref_region['km_total'] = great_circle_lengths

        print('simplifying CUSP data...')
        cusp_region_simp = cusp_region.copy()
        cusp_region_simp['geometry'] = cusp_region.geometry.simplify(tolerance=simp, 
                                                                     preserve_topology=True)
        print(' (saving to geopackage...)')
        layer = 'CUSP_{}_simplified'.format(region_id)
        cusp_region_simp['geometry'].to_file(cusp_gpkg, layer=layer, driver='GPKG')

        print('buffering simplified CUSP data...')
        cusp_region_simp_buff = buffer_lat_lon_multiline(cusp_region_simp.geometry, buff)
        cusp_buffer_gdf = gpd.GeoDataFrame(geometry=[cusp_region_simp_buff]).explode()
        print(' (saving to geopackage...)')
        layer = 'CUSP_{}_buffered'.format(region_id)
        cusp_buffer_gdf.to_file(cusp_gpkg, layer=layer, driver='GPKG')

        print('clipping reference shoreline with buffered simplified CUSP data...')
        ref_region_clipped = ref_region.copy()
        ref_region_clipped.geometry = ref_region.intersection(cusp_region_simp_buff)
        print(' (measuring great-circle lengths...)')
        great_circle_lengths = [calc_multiline_length(g) for g in ref_region_clipped['geometry']]
        ref_region_clipped['km_mapped'] = great_circle_lengths

        print('summing regional stats...')
        ref_region_lengths = ref_region.groupby('State')['km_total'].sum()
        ref_region_clipped_lengths = ref_region_clipped.groupby('State')['km_mapped'].sum()

        df = pd.DataFrame({'km_total': ref_region_lengths / 1000,
                           'km_mapped': ref_region_clipped_lengths / 1000,
                           'pct_mapped': ref_region_clipped_lengths / ref_region_lengths,
                           'region': [region] * ref_region_lengths.shape[0]})

        print(df.round(rounding).astype(types))
        results.append(df)

    results_df = pd.concat(results).round(rounding).astype(types)
    print_region_header('ALL PROCESSED REGIONS', '=')
    print(results_df)
    results_df.to_csv('{}\CUSP_Progress.txt'.format(out_dir), sep='\t')
    print('TOTAL TIME: {}'.format(datetime.now() - tic))
