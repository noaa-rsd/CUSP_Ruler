import os
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import datetime


def set_env_vars():
    user_dir = os.path.expanduser('~')

    script_path = Path(user_dir).joinpath('AppData', 'Local', 'Continuum', 
                                          'anaconda3', 'Scripts')

    gdal_data = Path(user_dir).joinpath('AppData', 'Local', 'Continuum', 
                                        'anaconda3', 'envs', 'qchecker', 
                                        'Library', 'share', 'gdal')

    proj_lib = Path(user_dir).joinpath('AppData', 'Local', 'Continuum', 
                                       'anaconda3', 'envs', 'qchecker', 
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
    width = 50
    sep1_num = int(width / 2 - int(len(region) / 2))
    sep2_num = int(width - len(region) - sep1_num)
    sep1 = sep * sep1_num
    sep2 = sep * sep2_num
    print('\n{}{}{}'.format(sep1, region.upper(), sep2))    


def get_tool_parameters():
    reference = arcpy.GetParameterAsText(0)
    cusp = arcpy.GetParameterAsText(1)
    out_dir = arcpy.GetParameterAsText(2)
    simp = int(arcpy.GetParameterAsText(3))
    buff = int(arcpy.GetParameterAsText(4))
    return reference, cusp, out_dir, simp, buff


if __name__ == '__main__':

    tic = datetime.datetime.now()

    set_env_vars()

    reference, cusp, out_dir, simp, buff = get_tool_parameters()

    web_mercator_epsg = {'init': 'epsg:3857'}
    wgs84_epsg = {'init': 'epsg:4326'}
    
    print('reading {}...'.format(reference))
    layer = reference.split(os.sep)[-1]
    ref_gdf = gpd.read_file(reference, layer=layer, crs=wgs84_epsg).to_crs(web_mercator_epsg)

    print('reading {}...'.format(cups))
    layer = cusp.split(os.sep)[-1]
    cusp_gdf = gpd.read_file(cusp, layer=layer, crs=wgs84_epsg)

    results = []

    for region in cusp_gdf.NOAA_Regio.unique():

        print_region_header(region, '-')
        print('extracting CUSP data...')
        cusp_region = cusp_gdf[cusp_gdf.NOAA_Regio == region].to_crs(web_mercator_epsg)

        print('extracting 70k data...')
        ref_region = ref_gdf[ref_gdf.NOAA_REGIO == region].to_crs(web_mercator_epsg)
        ref_region['length_km_total'] = ref_region.geometry.length

        print('simplifying CUSP data...')
        cusp_region_simp = cusp_region.copy()
        cusp_region_simp['geometry'] = cusp_region.geometry.simplify(tolerance=tolerance, preserve_topology=True)
        
        print('buffering simplified CUSP data...')
        cusp_region_simp_buff = cusp_region_simp.buffer(buffer).unary_union
        cusp_buffer_plot = gpd.GeoDataFrame(geometry=[cusp_region_simp_buff])

        print('clipping ref with bufferd simplified CUSP data...')
        ref_region_clipped = ref_region.copy()
        ref_region_clipped.geometry = ref_region.intersection(cusp_region_simp_buff)
        ref_region_clipped['length_km_mapped'] = ref_region_clipped.geometry.length
        
        print('summing regional stats...')
        ref_region_lengths = ref_region.groupby('State')['length_km_total'].sum()
        ref_region_clipped_lengths = ref_region_clipped.groupby('State')['length_km_mapped'].sum()

        df = pd.DataFrame({
            'km_total': ref_region_lengths / 1000,
            'km_mapped': ref_region_clipped_lengths / 1000,
            'pct_mapped': ref_region_clipped_lengths / ref_region_lengths,
            'region': [region] * ref_region_lengths.shape[0],
            })
        print(df)
        results.append(df)

    rounding = {
        'km_total': 0,
        'km_mapped': 0,
        'pct_mapped': 2
        }

    types = {
        'km_total': 'int64',
        'km_mapped': 'int64'
        }

    results_df = pd.concat(results).round(rounding).astype(types)
    print()
    print_region_header('ALL PROCESSED REGIONS', '=')
    print(results_df)
    results_df.to_csv('{}\CUSP_Progress.txt'.format(out_dir), sep='\t')
    print()
    print('TOTAL TIME: {}'.format(datetime.datetime.now() - tic))