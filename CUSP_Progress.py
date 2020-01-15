from pathlib import Path
from shutil import copytree
import geopandas as gpd


cusp_dir = Path(r'W:\Versions\Archive')
out_dir = Path (r'Z:\CUSP_AK')

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

for fy, gdb in cusps.items():
    cusp_gdb = cusp_dir / gdb
    out_gdb = out_dir / gdb

    print('-' * 50)
    
    #print('copying {} to {}...'.format(cusp_gdb, out_gdb))
    #copytree(str(cusp_path), str(out_path))

    layer = 'shoreline'
    print('opening {}...'.format(out_gdb / layer))
    wgs84_epsg = {'init': 'epsg:4326'}
    cusp_gdf = gpd.read_file(out_gdb, layer=layer, crs=wgs84_epsg)
    print(cusp_gdf.columns)
    