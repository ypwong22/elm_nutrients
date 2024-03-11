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

## edit leaf-root ratio to achieve better photosynthesis
hr['froot_leaf'][3] = 0.15
hr['froot_leaf'][11] = 0.15

# reduce the root N to boreal PFT levels in Ben's paper (2019) and should reduce the MR
hr['frootcn'][2] = 45
hr['frootcn'][3] = 45
hr['frootcn'][11] = 45

# The longevity of black spruce and tamarack both appear rather low
# Ruess, R.W., Hendrick, R.L., Burton, A.J., Pregitzer, K.S., Sveinbjornssön, B., Allen, M.F. and Maurer, G.E. (2003), COUPLING FINE ROOT DYNAMICS WITH ECOSYSTEM CARBON CYCLING IN BLACK SPRUCE FORESTS OF INTERIOR ALASKA. Ecological Monographs, 73: 643-662. https://doi.org/10.1890/02-4032
hr['froot_long'][2] = 1.5
hr['froot_long'][3] = 1.5
hr['froot_long'][11] = 1.5 # from MWM model, Siya Shao 2022 New Phytologist

# multiplication scalars on the nutrients uptake
hr['compet_pft_sminn'] = xr.DataArray(
    [
        np.nan, np.nan, 2.5, 1.5, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 6, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "gN gC-1 s-1", "long_name": "N uptake per unit fine root biomass"},
)

hr['compet_pft_sminp'] = xr.DataArray(
    [
        np.nan, np.nan, 2.5, 1.5, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 6, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "gP gC-1 s-1", "long_name": "P uptake per unit fine root biomass"},
)

hr['cpool_pft_sminn'] = xr.DataArray(
    [
        np.nan, np.nan, 1, 1, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 1, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "", "long_name": "Multiplicative factor on the CPOOL-drive mycorrhizae N uptake"},
)

hr['cpool_pft_sminp'] = xr.DataArray(
    [
        np.nan, np.nan, 1, 1, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 1, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "", "long_name": "Multiplicative factor on the CPOOL-drive mycorrhizae P uptake"},
)

hr['alpha_fpg'] = xr.DataArray(
    [1], dims = ["allpft"], attrs = {"units": "", "long_name": "adjust the rate of decreasing dependence on mycorrhizae-driven uptake as soil N content increase"})

hr['alpha_fpg_p'] = xr.DataArray(
    [1], dims = ["allpft"], attrs = {"units": "", "long_name": "adjust the rate of decreasing dependence on mycorrhizae-driven uptake as soil P content increase"})

# Q10 equation from 
# Ghimire et al. 2016 Representing leaf and root physiological traits in CLM 
#                     improves global carbon and nitrogen cycling predictions
# q10_uptake & scale_uptake from the same paper
# tbase_uptake from Shao et al. 2023 Ericoid mycorrhizal fungi mediate the 
#                                    response of ombrotrophic peatlands to 
#                                    fertilization: a modeling study
hr['q10_uptake'] = 1.5
hr['tbase_uptake'] = xr.DataArray(
    [293.15], dims = ["allpft"], attrs = {"units": "K", "long_name": "Base temperature in root nutrients uptake Q10 function"})
hr['scale_uptake'] = xr.DataArray(
    [
        np.nan, np.nan, 10, 10, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 10, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "", "long_name": "Scale factors in root nutrients uptake Q10 function"},
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
        np.nan, np.nan, 0.5, 0.5, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 0.5, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "gN m-2", "long_name": "Half saturation point of nitrogen uptake"},
)
hr['kmin_puptake'] = xr.DataArray(
    [
        np.nan, np.nan, 0.5, 0.5, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 0.5, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
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
