"""
Parameters for the uptaed nutrients competition routine
"""
import xarray as xr
import numpy as np
import os
from utils.paths import *

path_parameter = os.path.join(os.environ["PROJDIR"], "E3SM", "inputdata", "atm", "datm7",
                              "CLM1PT_data", "SPRUCE_data")

orgfile = 'clm_params_SPRUCE_20231120_spruceroot.nc_CNP'
newfile = 'clm_params_SPRUCE_20231120_spruceroot.nc_npcompet'

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

# sensitivity of fine root to leaf ratio to nutrient limitation
# For ombrotrophic shrub, increasing mineral nutrient content means increasing ease 
# of extracting them using roots intead of mycorrhizae. Hence they increase the allocation
# to roots. 
hr['froot_leaf_slope'] = xr.DataArray(
    [
        0, np.nan, 1.4, 1, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 1.4, 0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "", "long_name": "Sensitivity of froot_leaf to Michaelis-Menten limitation"},
)

# Multiplication scalars on the nutrients uptake
# Those should be around 1, but we know larch is slightly more efficient at self-uptake than spruce.
# Shrub is probably bad at self-uptake but ericoid mycorrhizae fungi is much better than
# the ecto-fungi at uptake.
hr['compet_pft_sminn'] = xr.DataArray(
    [
        0, np.nan, 0.5, 0.5, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 1.5, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "gN gC-1 s-1", "long_name": "N uptake per unit fine root biomass"},
)

hr['compet_pft_sminp'] = xr.DataArray(
    [
        0, np.nan, 0.5, 0.5, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 1.5, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "gP gC-1 s-1", "long_name": "P uptake per unit fine root biomass"},
)

# Level of competitiveness that fungi gains for the trees in their uptake
# > 1 to be more competitive than the tree roots.
hr['cpool_pft_sminn'] = xr.DataArray(
    [
        0, np.nan, 1., 1.5, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 1.5, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "", "long_name": "Multiplicative factor on the CPOOL-drive mycorrhizae N uptake"},
)

hr['cpool_pft_sminp'] = xr.DataArray(
    [
        0, np.nan, 1., 1.5, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 1.5, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "", "long_name": "Multiplicative factor on the CPOOL-drive mycorrhizae P uptake"},
)

# Q10 equation from 
# Ghimire et al. 2016 Representing leaf and root physiological traits in CLM 
#                     improves global carbon and nitrogen cycling predictions
# q10_uptake & scale_uptake from the same pape○♥r
# tbase_uptake from Shao et al. 2023 Ericoid mycorrhizal fungi mediate the 
#                                    response of ombrotrophic peatlands to 
#                                    fertilization: a modeling study
hr['q10_uptake'] = xr.DataArray(
    [
        1, np.nan, 2, 2, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 2, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "", "long_name": "Q10 in root nutrients uptake Q10 function"},
)
hr['tbase_uptake'] = xr.DataArray(
    [283.15], dims = ["allpft"], attrs = {"units": "K", "long_name": "Base temperature in root nutrients uptake Q10 function"})
hr['scale_uptake'] = xr.DataArray(
    [10], dims=["allpft"],
    attrs={"units": "", "long_name": "Scale factor in root nutrients uptake Q10 function"},
)

# Michaelis-Menten parameters: set to about the average N&P concentration in my simulations
hr['kmin_nuptake'] = xr.DataArray(
    [
        1e20, np.nan, 1e-4, 1e-4, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 1, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "gN m-3", "long_name": "Half saturation point of nitrogen uptake"},
)
hr['kmin_puptake'] = xr.DataArray(
    [
        1e20, np.nan, 1e-9, 2e-6, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 1e-5, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "gP m-3", "long_name": "Half saturation point of phosphorus uptake"},
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