""" Extract diagnostics variables from the ELM simulations """
import itertools as it
import xarray as xr
import pandas as pd
import numpy as np
import os
from glob import glob
from utils.constants import *
from utils.analysis import *
from utils.paths import *

prefix = '20240316_2'
rootpheno = True
npcompet = True

#for prefix,rootpheno,npcompet in zip(['20231113', '20240311', '20240316'],
#                                     [False, False, True],  [False, True, True]):

var_list = {}
var_list['pft'] = ['TLAI','GPP','NPP','AGNPP','BGNPP',
                   'TOTVEGC','TOTVEGC_ABG','MR','GR','XR','LITFALL',
                    'SMINN_TO_NPOOL', 'SMINP_TO_PPOOL', 'CPOOL', 'NPOOL', 'PPOOL', 'BTRAN']
for pool in ['LEAF', 'FROOT', 'LIVESTEM', 'DEADSTEM', 'LIVECROOT', 'DEADCROOT']:
    var_list['pft'] = var_list['pft'] + [f'{pool}C',f'{pool}C_STORAGE',f'{pool}C_XFER',
                                         f'CPOOL_TO_{pool}C', f'{pool}C_XFER_TO_{pool}C']
    if not 'DEAD' in pool:
        var_list['pft'] = var_list['pft'] + [f'{pool}_MR']
for pool in ['LEAF', 'FROOT']:
    var_list['pft'].extend([f'{pool}C_TO_LITTER'])
if npcompet:
    var_list['pft'].extend(['FPG_PATCH', 'FPG_P_PATCH','PLANT_NDEMAND','PLANT_PDEMAND',
                            'PLANT_NABSORB','PLANT_PABSORB','PLANT_NFUNGI_PATCH',
                            'PLANT_PFUNGI_PATCH','ZWT_ROOT_PATCH'])
if rootpheno:
    var_list['pft'] = var_list['pft'] + [
        'ONSET_FLAG', 'OFFSET_FLAG', 'DORMANT_FLAG',
        'BGLFR_LEAF', 'BGLFR_FROOT', 'ONSET_FLAG_ROOT',
        'OFFSET_FLAG_ROOT', 'DORMANT_FLAG_ROOT',
        'FCUR_DYN', 'ONSET_FROOT_FNMIN', 'ONSET_FROOT_FW',
        'LFR_FROOT_TD', 'LFR_FROOT_WD']

var_list['col'] = [
    'TBOT', 'TSOI_30', 'H2OSOI_30', 'HR', 'NEE', 'FPG', 'FPG_P', 
    'ZWT', 'SMINN_30', 'SOLUTIONP_30', 'ACTUAL_IMMOB', 'ACTUAL_IMMOB_P',
    'FPI', 'FPI_P', 'RH2M']
if not npcompet:
    var_list['col'].extend(['FPG', 'FPG_P'])

var_list['const'] = []

collection_ts, collection_const = extract_sims(prefix, var_list)

path_out_2 = os.path.join(path_out, 'extract', prefix)
if not os.path.exists(path_out_2):
    os.mkdir(path_out_2)
print(path_out_2)
collection_ts.to_csv(os.path.join(path_out_2, 'analysis_ts.csv'))
collection_const.to_csv(os.path.join(path_out_2, 'analysis_const.csv'))
