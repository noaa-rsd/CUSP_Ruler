import os
import argparse
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime
#import arcpy


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


def plot_region():
    ax = ref_gdf.plot(color='red')
    cusp_region.plot(ax=ax, color='blue')
    cusp_region_simp.plot(ax=ax, color='cyan')
    cusp_buffer_plot.plot(ax=ax, color='cyan', alpha=0.5)
    ref_region_clipped.plot(ax=ax, color='gray')
    plt.show()


def export_map_results():
    layer = 'CUSP_simp_buff'
    gdb_path = 'Z:\CUSP_progress\CUSP_Progress.gdb'
    print('saving {}...'.format(gdb_path))
    cusp_buffer_plot.to_file(gdb_path, layer=layer, driver='FileGDB')

    
def print_region_header(region, sep):
    width = 54
    sep1_num = int(width / 2 - int(len(region) / 2))
    sep2_num = int(width - len(region) - sep1_num)
    sep1 = sep * sep1_num
    sep2 = sep * sep2_num
    print('\n{}{}{}'.format(sep1, region.upper(), sep2))    


#def get_tool_parameters():
#    reference = Path(arcpy.GetParameterAsText(0))
#    cusp = Path(arcpy.GetParameterAsText(1))
#    out_dir = Path(arcpy.GetParameterAsText(2))
#    simp = int(arcpy.GetParameterAsText(3))
#    buff = int(arcpy.GetParameterAsText(4))
#    is_arc = True
#    return reference, cusp, out_dir, simp, buff, is_arc


#def define_tool_parameters():
#    reference = Path(r'Z:\CUSP_progress\70K_shoreline_corrected\70K_Shoreline_4CUSP.gdb\Whole_US')
#    cusp = Path(r'Z:\CUSP_progress\20190717_contemporary_shoreline.gdb\shoreline')
#    out_dir = Path(r'Z:\CUSP_progress\Results')
#    simp = 500
#    buff = 500
#    is_arc = False
#    return reference, cusp, out_dir, simp, buff, is_arc


#def print_msg(msg):
#    if not is_arc:
#        print(msg)
#    else:
#        arcpy.AddMessage(msg)


def get_args_argsparse():
    parser = argparse.ArgumentParser(description='Description of your program')
    
    parser.add_argument('-r','--reference', help='70k feature class path', required=True)
    parser.add_argument('-c','--cusp', help='CUSP feature class path', required=True)
    parser.add_argument('-o','--outdir', help='output directory', required=True)
    parser.add_argument('-t','--tolerance', help='simplify tolerance (m)', required=False, default=500)
    parser.add_argument('-b','--buffer', help='buffer (m)', required=False, default=500)

    args = vars(parser.parse_args())

    ref = Path(args['reference'])
    cusp = Path(args['cusp'])
    out = Path(args['outdir'])
    simp = int(args['tolerance'])
    buff = int(args['buffer'])

    return ref, cusp, out, simp, buff


def get_args():
    reference = input('Enter 70k shoreline feature class path (in quotes if path has spaces):  ')
    cusp = input('Enter CUSP feature class path (in quotes if path has spaces):  ')
    out_dir = input('Enter output directory (in quotes if path has spaces):  ')
    simp = input('Enter CUSP simplification tolerane (m):  ')
    buff = input('Enter CUSP buffer value (m):  ')
    
    return Path(reference), Path(cusp), Path(out_dir), int(simp), int(buff)


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

    #reference, cusp, out_dir, simp, buff, is_arc = get_tool_parameters()
    #reference, cusp, out_dir, simp, buff, is_arc = define_tool_parameters()

    web_mercator_epsg = {'init': 'epsg:3857'}
    wgs84_epsg = {'init': 'epsg:4326'}
    
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
        print_region_header(region, '-')
        print('extracting CUSP data...')
        cusp_region = cusp_gdf[cusp_gdf.NOAA_Regio == region].to_crs(web_mercator_epsg)

        print('extracting 70k shoreline...')
        ref_region = ref_gdf[ref_gdf.NOAA_REGIO == region].to_crs(web_mercator_epsg)
        ref_region['km_total'] = ref_region.geometry.length

        print('simplifying CUSP data...')
        cusp_region_simp = cusp_region.copy()
        cusp_region_simp['geometry'] = cusp_region.geometry.simplify(tolerance=simp, preserve_topology=True)
       
        print('buffering simplified CUSP data...')
        cusp_region_simp_buff = cusp_region_simp.buffer(buff).unary_union
        cusp_buffer_plot = gpd.GeoDataFrame(geometry=[cusp_region_simp_buff])

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
