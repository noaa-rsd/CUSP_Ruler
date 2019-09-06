import os
import argparse
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime


def set_env_vars():
    user_dir = os.path.expanduser('~')

    script_path = Path(user_dir).joinpath('AppData', 'Local', 'Continuum', 
                                          'anaconda3', 'Scripts')

    gdal_data = Path(user_dir).joinpath('AppData', 'Local', 'Continuum', 
                                        'anaconda3', 'envs', 'cusp', 
                                        'Library', 'share', 'gdal')

    proj_lib = Path(user_dir).joinpath('AppData', 'Local', 'Continuum', 
                                       'anaconda3', 'envs', 'cusp', 
                                       'Library', 'share')

    if script_path.name not in os.environ["PATH"]:
        os.environ["PATH"] += os.pathsep + str(script_path)

    os.environ["GDAL_DATA"] = str(gdal_data)
    os.environ["PROJ_LIB"] = str(proj_lib)


def print_region_header(region, sep):
    width = 54
    sep1_num = int(width / 2 - int(len(region) / 2))
    sep2_num = int(width - len(region) - sep1_num)
    sep1 = sep * sep1_num
    sep2 = sep * sep2_num
    print('\n{}{}{}'.format(sep1, region.upper(), sep2))    


def get_args():
    reference = input('\nEnter 70k shoreline feature class path:\n>'.upper())
    cusp = input('\nEnter CUSP feature class path:\n>'.upper())
    out_dir = input('\nEnter output directory:\n>'.upper())
    simp = input('\nEnter CUSP simplification tolerance (meters)\n>'.upper())
    buff = input('\nEnter CUSP buffer value (meters)\n>'.upper())
    
    return Path(reference.strip()), Path(cusp.strip()), Path(out_dir.strip()), int(simp), int(buff)


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


if __name__ == '__main__':

    tic = datetime.now()
    print_splash()

    reference, cusp, out_dir, simp, buff = get_args()
    set_env_vars()

    web_mercator_epsg = {'init': 'epsg:3857'}
    wgs84_epsg = {'init': 'epsg:4326'}
    
    print()
    print('reading {}...'.format(reference))
    ref_gdb = str(reference.parent)
    ref_layer = reference.name
    ref_gdf = gpd.read_file(ref_gdb, layer=ref_layer, crs=wgs84_epsg)

    print('reading {}...'.format(cusp))
    cusp_gdb = str(cusp.parent)
    cusp_layer = cusp.name
    cusp_gdf = gpd.read_file(cusp_gdb, layer=cusp_layer, crs=wgs84_epsg)

    results = []
    rounding = {'km_total': 0, 'km_mapped': 0, 'pct_mapped': 2}
    types = {'km_total': 'int64', 'km_mapped': 'int64'}

    for region in cusp_gdf.NOAA_Regio.unique():
        region_id = ''.join([c.capitalize() for c in region.split(' ')])
        cusp_gpkg = out_dir / 'CUSP_{}.gpkg'.format(region_id)
        print_region_header(region, '-')
        print('extracting CUSP data...')
        cusp_region = cusp_gdf[cusp_gdf.NOAA_Regio == region].to_crs(web_mercator_epsg)

        print('extracting 70k shoreline...')
        ref_region = ref_gdf[ref_gdf.NOAA_REGIO == region].to_crs(web_mercator_epsg)
        ref_region['km_total'] = ref_region.geometry.length

        print('simplifying CUSP data...', end='')
        cusp_region_simp = cusp_region.copy()
        cusp_region_simp['geometry'] = cusp_region.geometry.simplify(tolerance=simp, preserve_topology=True)
        print('saving to geopackage...')
        cusp_simp_shp = out_dir / 'CUSP.gpkg'.format(region_id)
        layer = 'CUSP_{}_simplified'.format(region_id)
        cusp_region_simp['geometry'].to_file(cusp_gpkg, layer=layer, driver='GPKG')

        print('buffering simplified CUSP data...', end='')
        cusp_region_simp_buff = cusp_region_simp.buffer(buff).unary_union
        cusp_buffer_plot = gpd.GeoDataFrame(geometry=[cusp_region_simp_buff])
        print('saving to geopackage...')
        layer = 'CUSP_{}_buffered'.format(region_id)
        cusp_buffer_plot.to_file(cusp_gpkg, layer=layer, driver='GPKG')

        print('clipping 70k shoreline with bufferd simplified CUSP data...')
        ref_region_clipped = ref_region.copy()
        ref_region_clipped.geometry = ref_region.intersection(cusp_region_simp_buff)
        ref_region_clipped['km_mapped'] = ref_region_clipped.geometry.length
      
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
    print()
    print('TOTAL TIME: {}'.format(datetime.now() - tic))
