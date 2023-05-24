""" Extract the time series of multiple runs and variables. """
import itertools as it
import xarray as xr
import pandas as pd
import numpy as np
import os
from glob import glob
from utils.constants import *
from utils.analysis import *
from utils.paths import *

"""
prefix = '20221212'
var_list = {}
var_list['pft'] = ['TSOI_AVG', 'SWC_AVG', 'ONSET_FLAG', 'FROOTC', 'FROOTC_STORAGE', 'DOWNREG']
var_list['col'] = ['TSOI_4', 'SWC_4', 'SMP_MAX', 'H2OSFC']
var_list['const'] = ['BSW', 'WATSAT', 'SUCSAT']

collection_ts, collection_const = extract_sim_one(prefix, var_list)

path_out_2 = os.path.join(path_out, 'extract', prefix)
if not os.path.exists(path_out_2):
    os.mkdir(path_out_2)
collection_ts.to_csv(os.path.join(path_out_2, 'for_root_ts.csv'))
collection_const.to_csv(os.path.join(path_out_2, 'for_root_const.csv'))
"""

for prefix in ['20230518']: # ['20221212', '20230120', '20230121', '20230122']:
    var_list = {}
    var_list['pft'] = ['TLAI', 'GPP', 'AGNPP', 'BGNPP', 'NPP', 'QVEGE', 'QVEGT', 'LEAFC', 'FROOTC', 'FROOTC_STORAGE', 'FROOTC_STORAGE_TO_XFER', 'FROOTC_XFER', 
                       'AR', 'FROOT_MR', 'CPOOL_FROOT_GR',  'CPOOL_FROOT_STORAGE_GR', 'TOTVEGC', 'FROOTC_ALLOC', 'LEAF_MR', 'CPOOL_LEAF_GR', 
                       'LEAFC', 'LEAFC_XFER', 'LEAFC_STORAGE', 'DOWNREG', 'ONSET_FLAG', 'ONSET_FLAG_ROOT', 'DORMANT_FLAG_ROOT', 
                       'ONSET_RATE_FROOT', 'COMPS_RATE_FROOT', 'BGLFR_FROOT', 'LEAFC_TO_LITTER', 'FROOTC_TO_LITTER']
    # 1 = near surface, 5 = 20 cm, 4 = 10 cm, 3 = 6 cm
    var_list['col'] = ['TBOT', 'TSOI_3', 'TSOI_4', 'H2OSOI_1', 'H2OSOI_5', 'SWC_1', 'SWC_5', 'HR', 'NEE']
    var_list['const'] = []

    collection_ts, collection_const = extract_sims(prefix, var_list)

    path_out_2 = os.path.join(path_out, 'extract', prefix)
    if not os.path.exists(path_out_2):
        os.mkdir(path_out_2)
    collection_ts.to_csv(os.path.join(path_out_2, 'analysis_ts.csv'))
    collection_const.to_csv(os.path.join(path_out_2, 'analysis_const.csv'))
