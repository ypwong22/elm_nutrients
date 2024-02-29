""" Extract the time series of multiple runs and variables. """
import itertools as it
import xarray as xr
import pandas as pd
import numpy as np
import os
from glob import glob
from utils.constants import *
from utils.analysis import *
from utils.paths import *

"""
prefix = '20221212'
var_list = {}
var_list['pft'] = ['TSOI_AVG', 'SWC_AVG', 'ONSET_FLAG', 'FROOTC', 'FROOTC_STORAGE', 'DOWNREG']
var_list['col'] = ['TSOI_4', 'SWC_4', 'SMP_MAX', 'H2OSFC']
var_list['const'] = ['BSW', 'WATSAT', 'SUCSAT']

collection_ts, collection_const = extract_sim_one(prefix, var_list)

path_out_2 = os.path.join(path_out, 'extract', prefix)
if not os.path.exists(path_out_2):
    os.mkdir(path_out_2)
collection_ts.to_csv(os.path.join(path_out_2, 'for_root_ts.csv'))
collection_const.to_csv(os.path.join(path_out_2, 'for_root_const.csv'))
"""

for prefix in [
    "20231112"
]:  # ['20221212', '20230120', '20230121', '20230122', '20230526', '20230623', "20230720"]:
    var_list = {}
    var_list["pft"] = [
        "TLAI",
        "GPP",
        "AGNPP",
        "BGNPP",
        "NPP",
        "AR",
        "TOTVEGC",
        "DOWNREG",
    ]  # "ONSET_FLAG", "OFFSET_FLAG", "DORMANT_FLAG", "QVEGE", "QVEGT"
    # "PSNSUN", "PSNSHA",

    for pool in ["LEAF", "FROOT", "LIVESTEM", "DEADSTEM", "LIVECROOT", "DEADCROOT"]:
        var_list["pft"] = var_list["pft"] + [
            f"{pool}C",
            f"{pool}C_STORAGE",
            f"{pool}C_XFER",
        ]

    for pool in ["LEAF", "FROOT"]:
        var_list["pft"].extend(
            [
                f"{pool}C_ALLOC",
                f"{pool}_MR",
                f"CPOOL_TO_{pool}C",
                f"CPOOL_TO_{pool}C_STORAGE",
                f"{pool}C_XFER_TO_{pool}C",
            ]
        )

    var_list["pft"].extend(["GRESP_XFER", "GRESP_STORAGE"])

    if prefix in ["20230623"]:
        var_list["pft"] = var_list["pft"] + [
            "BGLFR_LEAF",
            "BGLFR_FROOT",
            "LEAFC_TO_LITTER",
            "FROOTC_TO_LITTER",
            "ONSET_FLAG_ROOT",
            "OFFSET_FLAG_ROOT",
            "DORMANT_FLAG_ROOT",
            "FCUR_DYN",
            "ONSET_FROOT_FNMIN",
            "ONSET_FROOT_FW",
            "LFR_FROOT_TD",
            "LFR_FROOT_WD",
        ]
    # 1 = near surface, 5 = 20 cm, 4 = 10 cm, 3 = 6 cm
    var_list["col"] = [
        "TBOT",
        "TSOI_3",
        "TSOI_4",
        "H2OSOI_1",
        "H2OSOI_5",
        "SWC_1",
        "SWC_5",
        "HR",
        "NEE",
        "NET_NMIN",
        "BTRAN",
        "ZWT",
    ]
    var_list["const"] = []

    collection_ts, collection_const = extract_sims(prefix, var_list)

    path_out_2 = os.path.join(path_out, "extract", prefix)
    if not os.path.exists(path_out_2):
        os.mkdir(path_out_2)
    collection_ts.to_csv(os.path.join(path_out_2, "analysis_ts.csv"))
    collection_const.to_csv(os.path.join(path_out_2, "analysis_const.csv"))
