"""
Parameters for the uptaed nutrients competition routine
"""
import xarray as xr
import numpy as np
import os
from utils.paths import *

path_parameter = os.path.join(os.environ["PROJDIR"], "E3SM", "inputdata", "atm", "datm7",
                              "CLM1PT_data", "SPRUCE_data")

orgfile = 'clm_params_SPRUCE_20231120_spruceroot.nc_rootpheno'
newfile = 'clm_params_SPRUCE_20231120_spruceroot.nc_rootpheno_npcompet'
#orgfile = 'clm_params_SPRUCE_20231120_spruceroot.nc'
#newfile = 'clm_params_SPRUCE_20231120_spruceroot.nc_npcompet'

hr = xr.open_dataset(os.path.join(path_parameter, orgfile), decode_times=False)

# (1) based on the observations, the ratios mustn't exceed 1 and is on avg. ~0.5
# (2) 2004 Dissertation says Tamarak root:leaf ratio > Spruce root:leaf ratio
#     values taken from Fig. 6.4 for the second  growing season
#     also see abstract, tamarak allocates more to the root
hr['froot_leaf'][2] = 0.28
hr['froot_leaf'][3] = 0.38
hr['froot_leaf'][11] = 0.5

# sensitivity of fine root to leaf ratio to nutrient limitation
# For ombrotrophic shrub, increasing mineral nutrient content means increasing
# ease of extracting them using roots intead of mycorrhizae. Hence they reduce the 
# belowground allocation. 
# 2004 Dissertation p25: flooding reduces froot_leaf because roots cannot grow under hypoxia
# But larch is slightly better at handling hypoxia so say larch slope = 1
hr['froot_leaf_slope'] = xr.DataArray(
    [
        np.nan, np.nan, 2, 1, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, -2, 0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "", "long_name": "Sensitivity of froot_leaf to Michaelis-Menten limitation"},
)

# multiplication scalars on the nutrients uptake
hr['compet_pft_sminn'] = xr.DataArray(
    [
        np.nan, np.nan, 15, 30, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 45, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "gN gC-1 s-1", "long_name": "N uptake per unit fine root biomass"},
)

hr['compet_pft_sminp'] = xr.DataArray(
    [
        np.nan, np.nan, 15, 30, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 45, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "gP gC-1 s-1", "long_name": "P uptake per unit fine root biomass"},
)

hr['cpool_pft_sminn'] = xr.DataArray(
    [
        np.nan, np.nan, 0.5, 1, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 2, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "", "long_name": "Multiplicative factor on the CPOOL-drive mycorrhizae N uptake"},
)

hr['cpool_pft_sminp'] = xr.DataArray(
    [
        np.nan, np.nan, 0.5, 1, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 2, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "", "long_name": "Multiplicative factor on the CPOOL-drive mycorrhizae P uptake"},
)

hr['alpha_fpg'] = xr.DataArray(
    [3], dims = ["allpft"], attrs = {"units": "", "long_name": "adjust the rate of decreasing dependence on mycorrhizae-driven uptake as soil N content increase"})

hr['alpha_fpg_p'] = xr.DataArray(
    [3], dims = ["allpft"], attrs = {"units": "", "long_name": "adjust the rate of decreasing dependence on mycorrhizae-driven uptake as soil P content increase"})

# Q10 equation from 
# Ghimire et al. 2016 Representing leaf and root physiological traits in CLM 
#                     improves global carbon and nitrogen cycling predictions
# q10_uptake & scale_uptake from the same paper
# tbase_uptake from Shao et al. 2023 Ericoid mycorrhizal fungi mediate the 
#                                    response of ombrotrophic peatlands to 
#                                    fertilization: a modeling study
hr['q10_uptake'] = xr.DataArray(
    [
        np.nan, np.nan, 3., 1, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 4, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "", "long_name": "Q10 in root nutrients uptake Q10 function"},
)
hr['tbase_uptake'] = xr.DataArray(
    [293.15], dims = ["allpft"], attrs = {"units": "K", "long_name": "Base temperature in root nutrients uptake Q10 function"})
hr['scale_uptake'] = xr.DataArray(
    [10], dims=["allpft"],
    attrs={"units": "", "long_name": "Scale factor in root nutrients uptake Q10 function"},
)

# Michaelis-Menten parameters
# N-value taken from Ghimire et al. 2016 Representing leaf and root physiological traits in CLM
#                    improves global carbon and nitrogen cycling predictions
# P-value taken from Wang et al. (2020). Coupling of Phosphorus Processes With Carbon and Nitrogen 
#                    Cycles in the Dynamic Land Ecosystem Model: Model Structure, Parameterization,
#                    and Evaluation in Tropical Forests. Journal of Advances in Modeling Earth 
#                    Systems, 12(10), e2020MS002123. https://doi.org/10.1029/2020MS002123
hr['kmin_nuptake'] = xr.DataArray(
    [
        np.nan, np.nan, 10, 5, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 10, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "gN m-2", "long_name": "Half saturation point of nitrogen uptake"},
)
hr['kmin_puptake'] = xr.DataArray(
    [
        np.nan, np.nan, 10, 5, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 10, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "gP m-2", "long_name": "Half saturation point of phosphorus uptake"},
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