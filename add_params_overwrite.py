# 20230120 - reduce the photosynthesis parameters
import xarray as xr
import numpy as np
import os
from utils.paths import *
from utils.analysis import *


parm_best = os.path.join(os.environ['HOME'], 'models', 'OLMT', 'UQ_output', 'UQ_default_US-SPR_ICB20TRCNPRDCTCBC', 'MCMC_output', 'parms_best.txt')

prefix = os.path.join(os.environ['PROJDIR'], 'E3SM', 'inputdata', 'atm', 'datm7', 'CLM1PT_data', 'SPRUCE_data')

suffix_new = '20230510_UQ_default_XGBClassifier_MLPRegressor'

for file in ['clm_params.nc_yang_dmr_02242021']: # ['clm_params.nc_yang_dmr_yw_20230212', 'clm_params.nc_yang_dmr_02242021', 'clm_params.nc_yang_dmr_yw_20221231']
    hr = xr.open_dataset(os.path.join(prefix, file), decode_times = False)

    f = open(parm_best, 'r')
    for s in f:
        name, pft, val = s.split()
        pft = int(pft)
        val = float(val)
        if pft == 0:
            if name == 'q10_mr':
                # Catch a pitfall: q10_mr is modified globally but is a PFT level parameter
                hr[name][:] = val
            else:
                hr[name] = val
        else:
            hr[name][pft] = val
    f.close()

    encoding = {}
    for data_var in hr.data_vars:
        if '_FillValue' in hr[data_var].encoding.keys():
            continue
        else:
            encoding[data_var] = {'_FillValue': None}

    suffix = file.split('_')[-1]

    """
    if suffix == '20230212':
        suffix_new = '20230120_root'
    elif suffix == '02242021':
        suffix_new = '20230120'
    elif suffix == '20221231':
        suffix_new = '20230120_leaf'
    """

    hr.to_netcdf(os.path.join(prefix, file.replace(suffix, suffix_new)), encoding = encoding, format = 'NETCDF3_CLASSIC')
    hr.close()