""" Update some CNP ratios based on observational data"""
import xarray as xr
import numpy as np
import os
from utils.paths import *

path_parameter = os.path.join(os.environ["PROJDIR"], "E3SM", "inputdata", "atm", "datm7",
                              "CLM1PT_data", "SPRUCE_data")

orgfile = 'clm_params_SPRUCE_20231120_spruceroot.nc'
newfile = 'clm_params_SPRUCE_20231120_spruceroot.nc_CNP'

hr = xr.open_dataset(os.path.join(path_parameter, orgfile), decode_times=False)

# (1) based on the observations, the ratios mustn't exceed 1 and is on avg. ~0.5
# 2004 Dissertation says Tamarak root:leaf ratio > Spruce root:leaf ratio
#     values taken from Fig. 6.4 for the second  growing season
#     also see abstract, tamarak allocates more to the root
# Islam, MD. A. (2004). Ecophysiological adaptations of black spruce (Picea mariana) and tamarack (Larix laricina) seedlings to flooding and nutrition stress. University of Alberta, Edmonton, Alberta, Canada.
hr['froot_leaf'][2] = 0.3 # 0.28 - 0.5
hr['froot_leaf'][3] = 0.38 # 0.38 - 0.5
hr['froot_leaf'][11] = 0.3 # try to start low to match observation
# "Observations" at SPRUCE assume Picea mariana needles live 5 years long
# Salmon, V. G., Brice, D. J., Bridgham, S., Childs, J., Graham, J., Griffiths, N. A., et al. (2021). Nitrogen and phosphorus cycling in an ombrotrophic peatland: a benchmark for assessing change. Plant and Soil, 466(1–2), 649–674. https://doi.org/10.1007/s11104-021-05065-x
hr['leaf_long'][2] = 5

# attempted parameter fixes
hr['br_mr_pft'][3] = 8e-06 # increase the MR of tamarack to bring down biomass
hr['flnr'][2] = 0.284240288420581 # increase the photosynthesis of shrub to increase AGNPP
hr['br_mr_pft'][3] = 8e-06 # increase the MR of shrub to bring down biomass

## reduce SPRUCE's base maintenance respiration but increase Q10
hr['leaf_long'][2] = 5 # according to Paul's assumption
hr['br_mr_pft'][2] = 2e-06
hr['q10_mr_pft'][2] = 2.5
## reduce tamarack as an association
## increase shrub's base maintenance respiration
hr['br_mr_pft'][11] = 8e-06

# (3) based on fine root CN ratio observations
# @ Iversen2013RootIngrowthCoresS1Bog_RootTraits.csv (.xlsx)
# https://doi.org/10.25581/spruce.091/1782483
hr['frootcn'][2] = 35
hr['frootcn'][3] = 40
hr['frootcn'][11] = 55

# (4) based on twig CN ratio observations (Phillips et al.)
# @ SPRUCE_Plant_Tissue_Analyses_2009-2013_20170330.csv
# http://dx.doi.org/10.3334/CDIAC/spruce.038
# (leaf CN already reflected in the default parameter file)
# (leaf CP isn't complete, but look pretty uniform across species; use Pima's for lala)
hr['livewdcn'][2] = 90
hr['livewdcn'][3] = 60
hr['livewdcn'][11] = 85
hr['leafcp'][2] = 655.78 # But the default values look adjusted
hr['leafcp'][3] = 655.78
hr['leafcp'][11] = 594.72 # use Chaemodaphe because it's the largest fraction

# (5) Match the 3rd & 4th pool to the deep peat CN ratio. The
# 1st and 2nd pool should have higher CN ratios. Set to 25.
hr['cn_s1'] = 25 # 21-25
hr['cn_s2'] = 25 # 21-25
hr['cn_s3'] = 20
hr['cn_s4'] = 20

hr['rf_s3s4'] = 0.83 # increase the fraction lost to CO2 to better match HR

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