import xarray as xr
import numpy as np
import os
from utils.paths import *

path_parameter = os.path.join(
    os.environ["SCRATCHDIR"],
    "E3SM",
    "inputdata",
    "atm",
    "datm7",
    "CLM1PT_data",
    "SPRUCE_data",
)

orgfile = (
    "clm_params.nc_yang_dmr_20230518_UQ_default2_optimized_XGBClassifier_MLPRegressor"
)
# newfile = "clm_params.nc_yang_dmr_20230518_UQ_default2_optimized_XGBClassifier_MLPRegressor_root_phased_refactored"
newfile = "clm_params.nc_yang_dmr_20230518_UQ_default2_optimized_XGBClassifier_MLPRegressor_root_phased_refactored_NP"

hr = xr.open_dataset(os.path.join(path_parameter, orgfile), decode_times=False)


hr["rf_scale"] = -0.0441
hr["rf_tbase"] = xr.DataArray(
    [6.848 + 273.15],
    dims=["allpft"],
    attrs={
        "units": "K",
        "long_name": "Base temperature in uniforc model",
    },
)
hr["crit_onset_rf"] = 66.943
hr["ndays_on"] = xr.DataArray(
    [
        np.nan,
        np.nan,
        30,
        30,
        np.nan,
        np.nan,
        30,
        30,
        30,
        np.nan,
        30,
        30,
        np.nan,
        30,
        30,
        30,
        30,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "days", "long_name": "Number of days to complete leaf onset"},
)

hr["gdd_tbase"] = xr.DataArray(
    [
        np.nan,
        np.nan,
        np.nan,
        279.5,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        279.05,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "K", "long_name": "Base temperature for GDD accumulation"},
)

hr["crit_chil1"] = xr.DataArray(
    [
        np.nan,
        np.nan,
        np.nan,
        9.0,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        33,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "", "long_name": "Parameter in alternating model"},
)

hr["crit_chil2"] = xr.DataArray(
    [
        np.nan,
        np.nan,
        np.nan,
        2112,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        1388,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "", "long_name": "Parameter in alternating model"},
)

hr["crit_chil3"] = xr.DataArray(
    [
        np.nan,
        np.nan,
        np.nan,
        -0.04,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        -0.02,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "", "long_name": "Parameter in alternating model"},
)


hr["maxday_off"] = xr.DataArray(
    [286],
    coords={"allpfts": hr["allpfts"]},
    dims=["allpfts"],
    attrs={"units": "days", "long_name": "Day of maximum litterfall"},
)

hr["ndays_off"] = xr.DataArray(
    [
        np.nan,
        np.nan,
        48,
        15,
        np.nan,
        np.nan,
        15,
        15,
        15,
        np.nan,
        15,
        15,
        np.nan,
        15,
        15,
        15,
        15,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={
        "units": "days",
        "long_name": "Number of days to complete leaf litterfall (for evergreen this is only the relevant parameter not exact)",
    },
)

hr["crit_dayl"] = xr.DataArray(
    [
        np.nan,
        np.nan,
        np.nan,
        39300,
        np.nan,
        np.nan,
        36000,
        39300,
        39300,
        np.nan,
        36000,
        39300,
        np.nan,
        36000,
        36000,
        36000,
        36000,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "seconds", "long_name": "Critical day length for senescence"},
)

hr["crit_onset_rf_root"] = 42

hr["nmin_scale"] = xr.DataArray(
    [1],
    dims=["allpft"],
    attrs={
        "units": "",
        "long_name": "Parameter in soil moisture control on root growth and mortality",
    },
)

hr["wt_scale"] = xr.DataArray(
    [100],
    dims=["allpft"],
    attrs={
        "units": "mm",
        "long_name": "Parameter in soil moisture control on root growth and mortality",
    },
)

hr["mort_tsoi"] = xr.DataArray(
    [10 + 273.15],
    dims=["allpft"],
    attrs={
        "units": "K",
        "long_name": "Parameter in soil moisture control on root growth and mortality",
    },
)

hr["mort_a"] = 0.00175
hr["mort_b"] = 0.75

hr["mort_psi"] = xr.DataArray(
    [35],
    dims=["allpft"],
    attrs={
        "units": "",
        "long_name": "Parameter in soil moisture control on root growth and mortality",
    },
)

hr["mort_h2o"] = xr.DataArray(
    [0.6],
    dims=["allpft"],
    attrs={
        "units": "",
        "long_name": "Parameter in soil moisture control on root growth and mortality",
    },
)

hr["mort_d"] = 0.05

hr["ndays_on_root"] = xr.DataArray(
    [30],
    dims=["allpft"],
    attrs={
        "units": "days",
        "long_name": "Number of days to complete root onset",
    },
)

hr["ndays_off_fcur"] = xr.DataArray(
    [138],
    dims=["allpft"],
    attrs={
        "units": "days",
        "long_name": "Number of days to decrease fcur_dyn from 1 to 0",
    },
)

hr["hardiness_root"] = xr.DataArray(
    [
        np.nan,
        np.nan,
        273.15 - 7.5,
        273.15 - 7.5,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        273.15 - 7.5,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={
        "units": "K",
        "long_name": "Base temperature for GDD accumulation for root",
    },
)

hr["crit_chil2_root"] = xr.DataArray(
    [
        np.nan,
        np.nan,
        np.nan,
        1200.0,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        400.0,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={
        "units": "K",
        "long_name": "Scale for the critical onset gdd of the alternating model",
    },
)


# additional parameters

hr["froot_long"][2] = 0.6
hr["froot_long"][3] = 0.6
hr["froot_long"][11] = 0.6

hr["fcur"][2] = 0.2
hr["fcur"][3] = 0.0
hr["fcur"][11] = 0.0

hr["off_pstart"] = xr.DataArray(
    [
        np.nan,
        np.nan,
        np.nan,
        46800.0,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        54600.0,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "", "long_name": "Paramter for deciduous offset"},
)

hr["off_pend"] = xr.DataArray(
    [
        np.nan,
        np.nan,
        np.nan,
        1750.0,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        1600.0,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "", "long_name": "Paramter for deciduous offset"},
)

hr["off_tbase"] = xr.DataArray(
    [
        np.nan,
        np.nan,
        np.nan,
        294.5,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        290.15,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
    ],
    coords={"pft": hr["pft"]},
    dims=["pft"],
    attrs={"units": "", "long_name": "Paramter for deciduous offset"},
)

if newfile.endswith("_NP"):
    """Relative demand allocation by fine root biomass"""
    hr["VMAX_PLANT_N_RD"] = xr.DataArray(
        [
            np.nan,
            np.nan,
            1e-7,
            1e-7,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            1e-7,
            1e-7,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ],
        coords={"pft": hr["pft"]},
        dims=["pft"],
        attrs={
            "units": "gN gC-1 s-1",
            "long_name": "Paramter for plant uptake of nitrogen proportional to fine root biomass",
        },
    )
    hr["VMAX_PLANT_P_RD"] = xr.DataArray(
        [
            np.nan,
            np.nan,
            1e-7,
            1e-7,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            1e-7,
            1e-7,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ],
        coords={"pft": hr["pft"]},
        dims=["pft"],
        attrs={
            "units": "gP gC-1 s-1",
            "long_name": "Paramter for plant uptake of phosphorus proportional to fine root biomass",
        },
    )

encoding = {}
for data_var in hr.data_vars:
    if "_FillValue" in hr[data_var].encoding.keys():
        continue
    elif np.any(np.isnan(hr[data_var].values)):
        encoding[data_var] = {"_FillValue": -1e20}
    else:
        encoding[data_var] = {"_FillValue": None}
hr.to_netcdf(
    os.path.join(path_parameter, newfile), encoding=encoding, format="NETCDF3_CLASSIC"
)
