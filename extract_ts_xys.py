import pandas as pd
import os
import numpy as np
from utils.constants import *
from utils.paths import *
from utils.plotting import *
from utils.analysis import *
from tqdm import tqdm

var_list = {}
var_list["pft"] = [
    "TLAI",
    "GPP",
    "AGNPP",
    "BGNPP",
    "NPP",
    "AR",
    "TOTVEGC",
]

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
    "ZWT",
]
var_list["const"] = []

collection_ts, collection_const = extract_xys(var_list)
collection_ts.to_csv(
    os.path.join(path_out, "extract", "20230414_spruceroot", "analysis_ts.csv")
)
collection_const.to_csv(
    os.path.join(path_out, "extract", "20230414_spruceroot", "analysis_const.csv")
)
