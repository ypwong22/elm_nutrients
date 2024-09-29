""" Update some parameters based on optimization against Paul's data """
import xarray as xr
import numpy as np
import os
from utils.paths import *

path_parameter = os.path.join(os.environ["E3SM_ROOT"], "inputdata", "atm", "datm7",
                              "CLM1PT_data", "SPRUCE_data")

orgfile = 'clm_params_SPRUCE_20231120_spruceroot.nc_CNP'
newfile = 'clm_params_SPRUCE_20231120_spruceroot.nc_CNP_optim'

hr = xr.open_dataset(os.path.join(path_parameter, orgfile), decode_times=False)

hr['mbbopt'][2] = 0.3
hr['mbbopt'][3] = 0.38
hr['mbbopt'][11] = 0.3

hr['vcmaxha'][2] = 0.3
hr['vcmaxha'][3] = 0.38
hr['vcmaxha'][11] = 0.3

hr['vcmaxhd'][2] = 0.3
hr['vcmaxhd'][3] = 0.38
hr['vcmaxhd'][11] = 0.3

hr['flnr'][2] = 0.3
hr['flnr'][3] = 0.38
hr['flnr'][11] = 0.3

hr['slatop'][2] = 0.3
hr['slatop'][3] = 0.38
hr['slatop'][11] = 0.3

hr['br_mr_pft'][2] = 0.3
hr['br_mr_pft'][3] = 0.38
hr['br_mr_pft'][11] = 0.3

hr['q10_mr_pft'][2] = 0.3
hr['q10_mr_pft'][3] = 0.38
hr['q10_mr_pft'][11] = 0.3

hr['froot_leaf'][2] = 0.3
hr['froot_leaf'][3] = 0.38
hr['froot_leaf'][11] = 0.3
