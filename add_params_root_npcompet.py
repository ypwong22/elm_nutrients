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

orgfile = 'clm_params_SPRUCE_20231120_spruceroot.nc_CNP'
newfile = 'clm_params_SPRUCE_20231120_spruceroot.nc_npcompet'
#orgfile = 'clm_params_SPRUCE_UQ_20231113_g00682.nc'
#newfile = 'clm_params_SPRUCE_UQ_20231113_g00682.nc_npcompet'

hr = xr.open_dataset(os.path.join(path_parameter, orgfile), decode_times=False)

# Wait on this for parameter optimization
"""# increase the productivity of spruce and let shrub be more than larch (of course that should be?)
# https://www.nature.com/articles/s41467-021-25163-9/figures/1
# in the boreal zone, flnr = 10-20%, but evergreen needleleaf forest sometimes reach 30%
# leaf mass per area is strongly negatively correlated with flnr
# so let's do a gradiet
hr['flnr'][2] = 0.10
hr['flnr'][3] = 0.28
hr['flnr'][11] = 0.28
"""

# sensitivity of fine root to leaf ratio to water table
# shrub = 1.4
hr['zwt_froot_a'] = xr.DataArray(
    [
        0, np.nan, 0, 0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 1.4, 0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "", "long_name": "sensitivity of froot_leaf to water table depth"},
)

hr['froot_radius'] = xr.DataArray(
    [
        0, np.nan, 0.02, 0.02, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 0.01, 0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "cm", "long_name": "fine root radius"},
)

hr['froot_density'] = xr.DataArray(
    [
        0, np.nan, 0.03, 0.03, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 0.03, 0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "gC cm-3", "long_name": "fine root density"},
)

# calibrated to replace the plant's total N content of
# (100% * leaf + 100% * froot + r_mort*rest) in 1 year for the deciduous species
# (20% * leaf + 20% * froot + r_mort*rest) in 1 year for Picea mariana, since the root & leaf longevity are 5 years
# FROOTC:LEAFC = 22% for Pima, 49% for Lala, 12% for shfrub
# FROOTC:(TOTVEGC-LEAFC-FROOTC) = 77.4% for Pima, 9% for Lala, 4.3% for shrub
# r_mort = 0.02 for Pima & Lala, 0.12 for shrub
# vmax * frootc * 0.01 / froot_radius**2 / froot_density * (86400*365) = totvegc / CN
# hence, 
# vmax = np.array([0.2 / 0.22 + 0.2 + 0.02 / 0.774, 1/0.49 + 1 + 0.02 / 0.09, 1/0.12 + 1 + 0.12 / 0.043]) / 0.01 * np.array([0.02,0.02,0.01])**2 * np.array([0.03,0.03,0.03]) / 86400 / 365 / np.array([67.52, 73.80, 53.32])
#      = np.array([6.39e-13, 1.68e-12, 2.16e-12])
hr['vmax_froot_n'] = xr.DataArray(
    [
        0, np.nan, 2.2155167694702e-12, 2.39107323423353e-11, np.nan, np.nan, np.nan, np.nan, 
        np.nan, np.nan, np.nan, 1.58785467901505e-11, 0, np.nan, np.nan, np.nan, np.nan, 
        np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "gN cm-2 s-1", "long_name": "maximum N uptake rate per unit area of fine root"},
)

# use N:P ratio of each plant
# The ratios printed from code are as follows
#           1           3 C:N   67.515413567007371      C:P   946.81194843536082
#           1           4 C:N   73.799250851819622      C:P   1004.8973966937480
#           1          12 C:N   51.325518820304183      C:P   832.73157270429067
# vmax = np.array([0.2 / 0.22 + 0.2 + 0.02 / 0.774, 1/0.49 + 1 + 0.02 / 0.09, 1/0.12 + 1 + 0.12 / 0.043]) / 0.01 * np.array([0.02,0.02,0.01])**2 * np.array([0.03,0.03,0.03]) / 86400 / 365 / np.array([946.81, 1004.90, 832.73])
#      = array([4.56122134e-14, 1.23558890e-13, 1.38502456e-13])
# Yet, reality still needs adjustment. 
hr['vmax_froot_p'] = xr.DataArray(
    [
        0, np.nan, 2.31299374647962e-13, 9.00270690350306e-13, np.nan, np.nan, np.nan, np.nan, 
        np.nan, np.nan, np.nan, 9.07929749348862e-13, 0, np.nan, np.nan, np.nan, np.nan, np.nan,
        np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "gP cm-2 s-1", "long_name": "maximum P uptake rate per unit area of fine root"},
)

# It's actually meaningless to have both - NH4 is the vast majority of N, and the 
# model subroutines do not distinguish between NH4 and NO3 in its N demand
# literature value - 0.14, but is too high except for layer 3
hr['km_froot_n'] = xr.DataArray(
    [7], dims=["allpft"],
    attrs={"units": "gN m-3", "long_name": "half saturation point for NH4 uptake rate by fine root"},
)

# literature value - 0.7, also way to high
# Set to extremely low values to create stability
hr['km_froot_p'] = xr.DataArray(
    [0.00495538937141517], dims=["allpft"],
    attrs={"units": "gP m-3", "long_name": "half saturation point for PO4 uptake rate by fine root"},
)

hr['q10_upt'] = xr.DataArray(
    [2], dims=["allpft"],
    attrs={"units": "", "long_name": "temperature sensitivity of nutrient uptake by fine root"},
)

hr['swc_opt'] = xr.DataArray(
    [0.6], dims=["allpft"],
    attrs={"units": "m3 m-3", "long_name": "optimal volumetric soil moisture content for nutrient uptake by mycorrhizal fungi"},
)

hr['alpha_fpg'] = xr.DataArray(
    [1.5], dims=["allpft"],
    attrs={"units": "", "long_name": "sensitivity of nutrient uptake to nutrient limitation"},
)

hr['zwt_fungi_a'] = xr.DataArray(
    [
        0, np.nan, 0.7, 0.5, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 1, 0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "", "long_name": "sensitivity of fungi inhibition to nutrients"},
)

hr['zwt_fungi_b'] = xr.DataArray(
    [
        0, np.nan, 0, 0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, -0.04, 0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "", "long_name": "sensitivity of fungi inhibition to nutrients"},
)

hr['zwt_fungi_c'] = xr.DataArray(
    [
        0, np.nan, 0, 0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 0, 0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "", "long_name": "sensitivity of fungi inhibition to nutrients"},
)

# calculate from the baseline of vmax_froot_n
# r = np.array([0.02, 0.02, 0.01])
# rho = 0.03
# vmax_froot_n * 0.01 / r**2 / rho
# = array([5.51293064e-09, 2.82589436e-08, 1.46261823e-07])
hr['vmax_fungi_din'] = xr.DataArray(
    [
        0, np.nan, 4.5977396756146e-09, 3.48329403706601e-08, np.nan, np.nan, np.nan, np.nan, 
        np.nan, np.nan, np.nan, 3.32894097101296e-09, 0, np.nan, np.nan, np.nan, np.nan, np.nan, 
        np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "gN g-1 s-1", "long_name": "maximum inorganic N uptake rate per gram mycorrhizal root"},
)
# calculate from the baseline of vmax_froot_p
# vmax_froot_n * 0.01 / r**2 / rho
# = array([5.51293064e-09, 2.82589436e-08, 1.46261823e-07])
hr['vmax_fungi_dip'] = xr.DataArray(
    [
        0, np.nan, 1.222934959755e-10, 2.75666458920476e-10, np.nan, np.nan, np.nan, np.nan, 
        np.nan, np.nan, np.nan, 8.0368607446859e-10, 0, np.nan, np.nan, np.nan, np.nan, np.nan, 
        np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "gP g-1 s-1", "long_name": "maximum inorganic P uptake rate per gram mycorrhizal root"},
)

#hr['km_fungi_din'] = xr.DataArray(
#    [1e-2], dims=["allpft"],
#    attrs={"units": "gN m-3", "long_name": "half saturation point for inorganic N uptake by mycorrhizal root"},
#)

#hr['km_fungi_dip'] = xr.DataArray(
#    [1e-3], dims=["allpft"],
# pl   attrs={"units": "gP m-3", "long_name": "half saturation point for inorganic P uptake by mycorrhizal root"},
#)

hr['km_nsc'] = xr.DataArray(
    [2], dims=["allpft"],
    attrs={"units": "", "long_name": "parameter controlling the impact of nonstructural carbohydrate saturation on nutrient uptake"},
)

# baseline value should be identical to organic uptake
hr['vmax_fungi_son'] = xr.DataArray(
    [
        0, np.nan, 2.55004240634888e-09, 1.02089436e-08, np.nan, np.nan, np.nan, np.nan,
        np.nan, np.nan, np.nan, 3.57476004390863e-08, 0, np.nan, np.nan, np.nan, np.nan, np.nan,
        np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "gN gC-1 s-1", "long_name": "rate of mining organic N by mycorrhizal fungi"},
)

hr['vmax_fungi_sop'] = xr.DataArray(
    [
        0, np.nan, 1.52007484344376e-11, 1.07918248129927e-10, np.nan, np.nan, np.nan, np.nan,
        np.nan, np.nan, np.nan, 1.12650425480598e-10, 0, np.nan, np.nan, np.nan, np.nan, np.nan,
        np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "gP gC-1 s-1", "long_name": "rate of mining organic P by mycorrhizal fungi"},
)

# The effect of the cost variables seem to be mainly on reducing the hollow biomass (TLAI),
# limited effect on hummock
hr['fungi_cost_n'] = xr.DataArray(
    [20], dims=["allpft"],
    attrs={"units": "gC gN-1", "long_name": "carbon cost of fungal N uptake"},
)

hr['fungi_cost_p'] = xr.DataArray(
    [200], dims=["allpft"],
    attrs={"units": "gP m-3", "long_name": "carbon cost of fungal P uptake"},
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