""" HR parameters that were found to work well in Xiaoying's simulations, plus a couple others """
import xarray as xr
import numpy as np
import os
from utils.paths import *
from utils.analysis import *

parm_best = os.path.join(os.environ['HOME'], 'models', 'OLMT', 'UQ_output', 'UQ_default_US-SPR_ICB20TRCNPRDCTCBC', 'MCMC_output', 'parms_best.txt')

prefix = os.path.join(os.environ['PROJDIR'], 'E3SM', 'inputdata', 'atm', 'datm7', 'CLM1PT_data', 'SPRUCE_data')

file = 'clm_params.nc_yang_dmr_02242021'
suffix_new = 'nc_yang_dmr_20230509'

# Optimized parameters https://bg.copernicus.org/articles/18/467/2021/#&gid=1&pid=1

hr = xr.open_dataset(os.path.join(prefix, file), decode_times = False)
hr['q10_hr'] = 3

# // moss doesn't seem sensitive to those, still
hr['slatop'][12] = 0.0286 * 2 # 20230434
hr['grperc'][12] = 0.01 # 20230433
hr['br_mr'][12] = 2.5838138262092e-06 # 20230434
hr['q10_mr'][12] = 2.5 # 20230433

# based on FRED data
hr['frootcn'][2] = 48
hr['frootcn'][3] = 40
hr['frootcn'][11] = 60

encoding = {}
for data_var in hr.data_vars:
    if '_FillValue' in hr[data_var].encoding.keys():
        continue
    else:
        encoding[data_var] = {'_FillValue': None}

suffix = file.split('.')[-1]

hr.to_netcdf(os.path.join(prefix, file.replace(suffix, suffix_new)), encoding = encoding, format = 'NETCDF3_CLASSIC')
hr.close()