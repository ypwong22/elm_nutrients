""" Extract the time series of multiple runs and variables. """
import itertools as it
import xarray as xr
import pandas as pd
import numpy as np
import os
from glob import glob
from optparse import OptionParser
from utils.constants import *


"""
parser = OptionParser()
parser.add_option("--sim_id", dest="sim_id", default = '', \
                  help = 'Index to the simulation to extract.')
parser.add_option('--plot_id', dest = 'plot_id', default = '', \
                  help = 'Index to the chamber to extract.')
(options,args) = parser.parse_args()

prefix  =  sims_prefix[int( options.sim_id)]
plot    = chamber_list[int(options.plot_id)]
"""


def read_vars(prefix, plot):
    path_data = os.path.join(os.environ['PROJDIR'], 'E3SM', 'output', f'{prefix}_plot{plot:02d}_US-SPR_ICB20TRCNPRDCTCBC', 'run')

    data_list = {}
    for i, varname in enumerate(var_list):
        print(varname, unit_list[i])

        varname1 = varname.split('_')[0]
        if len(varname.split('_')) > 1:
            varname2 = varname.split('_')[1]
        else:
            varname2 = ''

        if varname2 in ['lala', 'pima', 'shrub', 'moss'] or varname1 == 'BGNPP':
            flist = sorted(glob(os.path.join(path_data, '*.h2.*.nc')))[:-1]
            hr = xr.open_mfdataset(flist, decode_times = False)
        else:
            flist = sorted(glob(os.path.join(path_data, '*.h1.*.nc')))[:-1]
            hr = xr.open_mfdataset(flist, decode_times = False)

        for kk, loc in zip([0,1], ['hummock', 'hollow']):
            if varname2 == 'lala':
                var = hr[varname1][:,  3 + 17*kk]
            elif varname2 == 'pima': 
                var = hr[varname1][:,  2 + 17*kk]
            elif varname2 == 'shrub':
                var = hr[varname1][:, 11 + 17*kk]
            elif varname2 == 'moss':
                var = hr[varname1][:, 12 + 17*kk]
            elif varname1 == 'BGNPP':
                var = hr[varname1][:, 2 + 17*kk] * 0.36 + hr[varname1][:, 3 + 17*kk] * 0.14 + hr[varname1][:, 11 + 17*kk]*0.25
            elif varname1 == 'H2OSOI':
                subset = np.where(hr['levgrnd'].values <= 0.4)[0][-1]
                thickgrnd = soil_interfaces[:(subset+1)].copy()
                thickgrnd[-1] = 0.4
                thickgrnd = np.diff(thickgrnd)
                var = (hr[varname1][:, :subset, kk] * \
                       thickgrnd.reshape(1,-1)).sum(dim = 'levgrnd') / thickgrnd.sum()
            elif varname1 == 'ZWT':
                if loc == 'hummock':
                    var = - hr[varname1][:, kk] + 0.15 # add the hummock height of 15cm
                else:
                    # add surface water
                    var = - hr[varname1][:, kk] + hr['H2OSFC'][:, kk].mean() / 1000
            else:
                var = hr[varname1][:, kk]
            if unit_list[i] in ['gC m-2 day-1', 'mm day-1']:
                var = var * 86400 # s-1 to day-1
            data_list[(loc, varname)] = pd.Series(var, index = sim_tvec)

            #print(data_list[(loc, varname)])

    hr.close()
    return data_list


prefix = '20230212' # '20221212' # '20221231'
for plot in [4, 10, 11, 16, 19, 6, 8, 13, 17, 20]:
    path_out = os.path.join(os.environ['PROJDIR'], 'Phenology_ELM', 'output_elm')
    data_list = read_vars(prefix, plot)
    pd.DataFrame(data_list).to_csv(os.path.join(path_out, prefix, f'plot{plot}_ts_extract.csv'))
