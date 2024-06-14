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
# 1st and 2nd pool should have higher CN ratios. Set to 22.
hr['cn_s1'] = 22 # 21-25
hr['cn_s2'] = 22 # 21-25
hr['cn_s3'] = 20
hr['cn_s4'] = 20

# (6) Edits based on parameter optimization outcomes 20231113_2a and 20231113_2b
hr['flnr'][3] = 1.541807362373865953e-01
hr['br_mr_pft'][3] = 8.275866311697777393e-06
hr['froot_leaf'][3] = 4.891957834556074358e-01
hr['stem_leaf'][3] = -3.076457295690046001e-01
hr['croot_stem'][3] = 7.389691614435485167e-01

hr['flnr'][11] = 2.083813227331323681e-01
hr['br_mr_pft'][11] = 4.995466990066244512e-06
hr['froot_leaf'][11] = 1.228784812082149847e-01
hr['stem_leaf'][11] = 1.509871452260563018e-01
hr['croot_stem'][11] = 1.678106396067260975e-01


# (7) make mortality a PFT-specific parameter
hr['r_mort'] = xr.DataArray(
    [
        0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.12, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "yr-1", "long_name": "Whole-plant mortality"},
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
