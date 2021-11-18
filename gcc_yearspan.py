""" Get the dates spanned by the evergreen GCC sites """
import os
import pandas as pd
import numpy as np
from glob import glob

path_root = os.path.join(os.environ['PROJDIR'], 'DATA', 'Vegetation', 'PhenoCam_V2_1674', 
                         'data', 'data_record_4')

sitelist = glob(os.path.join(path_root, '*_EN_*_1day.csv'))
sitenames = [s.split('/')[-1].split('_')[0] for s in sitelist]

startlist = []
endlist   = []
latlist   = []
lonlist   = []
for file in sitelist:
    data  = pd.read_csv(file, skiprows = 24, index_col = 0, parse_dates=True)['smooth_gcc_90']
    data  = data.dropna().sort_index()

    f     = open(file)
    found = 0
    while (found < 2):
        aline = f.readline().lower()
        if   'lat' in aline:
            lat = aline.replace(' ','').replace(':','').replace('#','').replace('\n','').replace('lat', '')
            found += 1
        elif 'lon' in aline:
            lon = aline.replace(' ','').replace(':','').replace('#','').replace('\n','').replace('lon', '')
            found += 1
    f.close()
    latlist.append(lat)
    lonlist.append(lon)

    start = data.index[ 0].strftime('%Y-%m-%d')
    end   = data.index[-1].strftime('%Y-%m-%d')
    startlist.append(start)
    endlist  .append(end)

result = pd.DataFrame({'site': sitenames, 'start': startlist, 'end': endlist, 
                       'lat': latlist, 'lon': lonlist})
result.to_csv('gcc_yearspan.csv', index = False)