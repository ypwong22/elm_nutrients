"""
cp /gpfs/wolf2/cades/cli185/scratch/ywo/E3SM/inputdata/atm/datm7/CLM1PT_data/SPRUCE_data/clm_params_SPRUCE_20231120_spruceroot.nc_root_process_phased_refact
or_alloc /gpfs/wolf2/cades/cli185/scratch/ywo/E3SM/output/20240303_US-SPR_ICB1850CNRDCTCBC_ad_spinup/run/clm_params.nc
"""

import xarray as xr
import numpy as np
import os
from utils.paths import *

path_parameter = os.path.join(os.environ["PROJDIR"], "E3SM", "inputdata", "atm", "datm7",
                              "CLM1PT_data", "SPRUCE_data")

orgfile = 'clm_params_SPRUCE_20231120_spruceroot.nc_root_process_phased_refactor'
#orgfile = 'clm_params_SPRUCE_20231120_spruceroot_eca.nc_root_process_phased_refactor'
newfile = 'clm_params_SPRUCE_20231120_spruceroot.nc_root_process_phased_refactor_alloc'

hr = xr.open_dataset(os.path.join(path_parameter, orgfile), decode_times=False)

# edit leaf-root ratio to achieve better photosynthesis
hr['froot_leaf'][3] = 0.01
hr['froot_leaf'][11] = 0.05
hr['froot_long'][3] = 3
hr['froot_long'][11] = 2

#
hr['compet_pft_sminn'] = xr.DataArray(
    [
        np.nan, np.nan, 1, 2, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 4, 1, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "gN gC-1 s-1", "long_name": "N uptake per unit fine root biomass"},
)

hr['compet_pft_sminp'] = xr.DataArray(
    [
        np.nan, np.nan, 1, 2, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 4, 1, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "gP gC-1 s-1", "long_name": "P uptake per unit fine root biomass"},
)

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