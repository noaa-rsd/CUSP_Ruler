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
    arcpy.AddMessage('\n{}{}{}'.format(sep1, region, sep2))    


def get_args(arc_params):
    ref_shoreline = str(arc_params[0].value)
    cusps = str(arc_params[1].value)
    data_sources = str(arc_params[2].value)
    out_dir = str(arc_params[3].value)
    
    simp = 0.002  # degrees
    buff_radius = 500  # meters

    return ref_shoreline, cusps, data_sources, out_dir, simp, buff_radius


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


def load_default_paths(defaults_path):
    with open(defaults_path) as f:
        paths = json.load(f)
    return paths


def save_config(reference, cusp_paths, data_sources, out_dir, simp, buff_radius):
    default_paths = {'reference_fc': reference, 
                     'cusp_fcs': cusp_paths.split(';'),
                     'data_sources': data_sources,
                     'output_dir': out_dir,
                     'simplification': simp,
                     'buff_radius': buff_radius}

    with open(defaults_path, 'w') as f:
        json.dump(default_paths, f, indent=1)


if __name__ == '__main__':
    tic = datetime.now()
    print_splash()

    cusp_ruler_dir = Path(os.path.dirname(os.path.realpath(__file__)))
    os.chdir(cusp_ruler_dir)
    defaults_path = cusp_ruler_dir / 'default_paths.txt'

    geod = Geodesic()

    default_paths = load_default_paths(defaults_path)
    arc_params = arcpy.GetParameterInfo()
    reference, cusp_paths, data_sources, out_dir, simp, buff_radius = get_args(arc_params)
    save_config(reference, cusp_paths, data_sources, out_dir, simp, buff_radius)

    support_dir = Path(r'../')
    support_gpkg = support_dir / 'CUSP_Ruler.gpkg'

    env_name = 'shorex2'
    set_env_vars(env_name)

    wgs84_epsg = {'init': 'epsg:4326'}
    
    ref_path = Path(reference)
    cusp_paths = [Path(c) for c in cusp_paths.split(';')]
    out_dir = Path(out_dir)

    arcpy.AddMessage(f'reading {str(ref_path)}...')
    ref_gdb = str(ref_path.parent)
    layer = ref_path.name.replace('main.', '')  # geopackage puts 'main.' in layer name
    ref_gdf = gpd.read_file(ref_gdb, layer=layer).to_crs(wgs84_epsg)

    results = []
    rounding = {'StateLeng': 0, 'ClippedLeng': 0, 'Percentage': 2}
    types = {'StateLeng': 'int64', 'ClippedLeng': 'int64'}
    
    cols_to_drop = ['level_0', 'level_1']  # artifacts from explode()

    for cusp_path in cusp_paths:
        cusp_gdb = str(cusp_path.parent)
        cusp_dataset = cusp_gdb.split('\\')[-1]
        cusp_layer = cusp_path.name

        print_region_header(cusp_dataset, '*')

        arcpy.AddMessage(f'reading {str(cusp_path)}...')
        cusp_gdf = gpd.read_file(cusp_gdb, layer=cusp_layer, crs=wgs84_epsg)

        if data_sources == 'GC & DM':
            gc_idx = cusp_gdf.SOURCE_ID.str[0:2] == 'GC'
            dm_idx = cusp_gdf.SOURCE_ID.str[0:2] == 'DM'
            cusp_gdf = cusp_gdf[gc_idx | dm_idx]

        for region in cusp_gdf.NOAA_Regio.unique()[0:1]:
            region_id = ''.join([c.capitalize() for c in region.split(' ')])
            cusp_gpkg = out_dir / 'CUSP_{}.gpkg'.format(region_id)
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
            layer = '{}_CUSP_{}_simplified'.format(cusp_layer, region_id)
            cusp.to_file(cusp_gpkg, layer=layer, driver='GPKG')

            arcpy.AddMessage('\nbuffering CUSP data...')
            cusp_simp_buff = buffer_lat_lon_multiline(cusp.geometry, buff_radius)
            gdf = gpd.GeoDataFrame(geometry=[cusp_simp_buff], crs=wgs84_epsg)
            cusp_buffer_gdf = gdf.explode().reset_index().drop(cols_to_drop, axis=1)

            layer = 'BufferBlocks'
            bb_gdf = gpd.read_file(str(support_gpkg), layer=layer, crs=wgs84_epsg)
            cusp_buffer_gdf = gpd.overlay(cusp_buffer_gdf, bb_gdf, how='intersection')

            layer = '{}_CUSP_{}_buffered'.format(cusp_layer, region_id)
            cusp_buffer_gdf.to_file(cusp_gpkg, layer=layer, driver='GPKG')

            arcpy.AddMessage('\nclipping reference shoreline with CUSP buffers...')
            ref_to_clip = ref.copy()
            ref_sidx = ref_to_clip.sindex

            refs_clipped = []
            for buff in cusp_buffer_gdf.geometry:
                possible_ref_idx = list(ref_sidx.intersection(buff.bounds))
                possible_ref = ref_to_clip.iloc[possible_ref_idx]

                for j, row in possible_ref.iterrows():
                    row.geometry = row.geometry.intersection(buff)
                    refs_clipped.append(row)

            refs_clipped = gpd.GeoDataFrame(refs_clipped, crs=wgs84_epsg)

            geodesic_lengths = [calc_line_length(g) for g in refs_clipped['geometry']]
            refs_clipped['km_mapped'] = geodesic_lengths

            arcpy.AddMessage('\nsumming regional stats...')
            ref_lengths = ref.groupby('State')['km_total'].sum()
            ref_clipped_lengths = refs_clipped.groupby('State')['km_mapped'].sum()

            df = pd.DataFrame({
                'StateLeng': ref_lengths / 1000,
                'ClippedLeng': ref_clipped_lengths / 1000,
                'Percentage': ref_clipped_lengths / ref_lengths,
                'noaaRegion': [region] * ref_lengths.shape[0],
                #'CUSP': [cusp_dataset] * ref_lengths.shape[0]
                })

            df.index.names = ['stateName']
            col_order = ['noaaRegion', 'ClippedLeng', 'StateLeng', 'Percentage']
            arcpy.AddMessage(df[col_order].round(rounding).astype(types))
            results.append(df[col_order])

        results_df = pd.concat(results).round(rounding).astype(types)
        print_region_header('ALL PROCESSED REGIONS', '=')
        arcpy.AddMessage(results_df)

        run_datetime = datetime.now()
        output_path = f'{out_dir}\CUSP_Progress_{run_datetime:%m%d%Y}.txt'
        results_df.to_csv(output_path, sep=',')
        format_output_file(output_path)

        arcpy.AddMessage('creating CUSP Ruler metadata file...')
        meta_path = str(output_path).replace('.txt', '.meta')
        with open(meta_path, 'w') as f:
            json.dump(f, default_paths, indent=1)

    arcpy.AddMessage('TOTAL TIME: {}'.format(datetime.now() - tic))

