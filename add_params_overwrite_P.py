"""
Parameters for the updated nutrients competition routine

Note: it seems parameter difference at the level of 1e-27 can still cause results to diverge. For exact replication, please copy the parameter file directly from UQ folder. 
"""
import xarray as xr
import numpy as np
import os
from utils.paths import *

path_parameter = os.path.join(os.environ["E3SM_ROOT"], "inputdata", "atm", "datm7",
                              "CLM1PT_data", "SPRUCE_data")

orgfile = 'clm_params_SPRUCE_UQ_20240112_g01944.nc_npcompet_cost0'
newfile = 'clm_params_SPRUCE_UQ_20240112_g01944.nc_npcompet_cost0_P'

hr = xr.open_dataset(os.path.join(path_parameter, orgfile), decode_times=False)

hr['np_s1_new'] = 40
hr['np_s2_new'] = 40
hr['np_s3_new'] = 40
hr['np_s4_new'] = 40

encoding = {}
for data_var in hr.data_vars:
    if "_FillValue" in hr[data_var].encoding.keys():
        continue # somehow I cannot drop this line or E3SM throws error
    elif np.any(np.isnan(hr[data_var].values)):
        encoding[data_var] = {"_FillValue": -1e20}
    else:
        encoding[data_var] = {"_FillValue": None}
hr.to_netcdf(
    os.path.join(path_parameter, newfile), encoding=encoding, format="NETCDF3_CLASSIC"
)