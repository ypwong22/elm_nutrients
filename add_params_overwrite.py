# 20230120 - reduce the photosynthesis parameters
import xarray as xr
import numpy as np
import os
from utils.paths import *
from utils.analysis import *


prefix = os.path.join(os.environ["PROJDIR"],"E3SM","inputdata","atm",
                      "datm7","CLM1PT_data","SPRUCE_data")

#experiment = "UQ_default2_1_optimized"
experiment = "UQ_default3"

parm_best = os.path.join(os.environ["HOME"],"models","OLMT",
    "UQ_output",
    f"{experiment}_US-SPR_ICB20TRCNPRDCTCBC",
    "MCMC_output",
    "parms_best.txt",
)

# suffix_new = '20230518_UQ_default2_optimized_XGBClassifier_MLPRegressor'
# suffix_new = '20230524_UQ_default2_MLPRegressor'
suffix_new = f"20230720_{experiment}_XGBClassifier_MLPRegressor_NP"

# for file in ['clm_params.nc_yang_dmr_02242021']: # ['clm_params.nc_yang_dmr_yw_20230212', 'clm_params.nc_yang_dmr_02242021', 'clm_params.nc_yang_dmr_yw_20221231']
file = "clm_params.nc_yang_dmr_20230509"

hr = xr.open_dataset(os.path.join(prefix, file), decode_times=False)

f = open(parm_best, "r")
for s in f:
    name, pft, val = s.split()
    pft = int(pft)
    val = float(val)
    if pft == 0:
        if name == "q10_mr":
            # Catch a pitfall: q10_mr is modified globally but is a PFT level parameter
            hr[name][:] = val
        else:
            hr[name] = val
    else:
        if not name in hr.variables:
            hr[name] = xr.DataArray(np.full(len(hr['pft']), np.nan), dims = ['pft'], coords = {'pft': hr['pft']})
        hr[name][pft] = val
f.close()


if suffix_new.endswith("_NP") and not ("VMAX_PLANT_N_RD" in hr.variables):
    """Relative demand allocation by fine root biomas
    Only the relative magnitudes across PFTs matter.
    Therefore one should "anchor" against one of the PFTs (always = 1)
    """
    hr["VMAX_PLANT_N_RD"] = xr.DataArray(
        [
            np.nan,
            np.nan,
            1,
            1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            2,
            1,
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
            1,
            1,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            2,
            1,
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

suffix = file.split("_")[-1]

"""
if suffix == '20230212':
    suffix_new = '20230120_root'
elif suffix == '02242021':
    suffix_new = '20230120'
elif suffix == '20221231':
    suffix_new = '20230120_leaf'
"""

hr.to_netcdf(
    os.path.join(prefix, file.replace(suffix, suffix_new)),
    encoding=encoding,
    format="NETCDF3_CLASSIC",
)
hr.close()
