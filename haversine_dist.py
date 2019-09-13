import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
from pathlib import Path
from datetime import datetime


def haversine(segments):
    segments = np.radians(segments)
    lon1, lat1 = segments[:, 0], segments[:, 1]
    lon2, lat2 = segments[:, 2], segments[:, 3]
    dlon = lon2 - lon1 
    dlat = lat2 - lat1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    c = 2 * np.arcsin(np.sqrt(a))
    r = 6371
    return c * r


def calc_multiline_length(multiline):
    multiline_length = 0
    for line in multiline:
        a = line.coords[:]
        line_segments = np.hstack((a[0:-1], a[1:]))
        line_length = np.sum(haversine(line_segments))
        multiline_length += line_length
    print('multiline_length = {}'.format(multiline_length))
    return multiline_length
    

tic = datetime.now()

wgs84_epsg = {'init': 'epsg:4326'}

reference = Path(r'Z:\CUSP_progress\70K_shoreline_corrected\70K_Shoreline_4CUSP.gdb\Whole_US')
print('reading {}...'.format(reference))
ref_gdb = str(reference.parent)
ref_layer = reference.name
ref_gdf = gpd.read_file(ref_gdb, layer=ref_layer, crs=wgs84_epsg)


ref_gdf['geometry'].apply(lambda g: calc_multiline_length(g))

print('TOTAL TIME: {}'.format(datetime.now() - tic))
print()

#ref_gdf['geometry'].apply(lambda x: ref_gdf['geometry'].apply(lambda y: haversine(x, y)))
#ref_gdf[['geometry']].apply(lambda x: ref_gdf.loc[:x.name, 'geometry'].apply(lambda y: haversine(x['geometry'], y)), axis=1)


