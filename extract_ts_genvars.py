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

prefix = '20240101'
ensemble_id = None
#prefix  = "UQ_20240327"
#ensemble_id = 2527
rootpheno = False
npcompet = True


#for prefix,rootpheno,npcompet in zip(['20231113', '20240311', '20240316'],
#                                     [False, False, True],  [False, True, True]):

var_list = {}
var_list['pft'] = ['TLAI','GPP','NPP','AGNPP','BGNPP','FROOTC_ALLOC','LEAFC_ALLOC','TOTVEGC',
                   'TOTVEGC_ABG','MR','GR','XR','LITFALL','FROOTN','FROOTP','LEAFN','LEAFP',
                   'SMINN_TO_NPOOL', 'SMINP_TO_PPOOL', 'CPOOL', 'NPOOL', 'PPOOL', 'BTRAN',
                   'PLANT_NDEMAND','PLANT_PDEMAND', 'RETRANSN_TO_NPOOL','RETRANSP_TO_PPOOL',
                   'XSMRPOOL','AVAILC']
for pool in ['LEAF', 'FROOT', 'LIVESTEM', 'DEADSTEM', 'LIVECROOT', 'DEADCROOT']:
    var_list['pft'] = var_list['pft'] + [f'{pool}C',f'{pool}C_STORAGE',f'{pool}C_XFER',
                                         f'CPOOL_TO_{pool}C', f'CPOOL_TO_{pool}C_STORAGE',
                                         f'{pool}C_XFER_TO_{pool}C']
    if not 'DEAD' in pool:
        var_list['pft'] = var_list['pft'] + [f'{pool}_MR']
for pool in ['LEAF', 'FROOT']:
    var_list['pft'].extend([f'{pool}C_TO_LITTER'])

# root fraction weighted environmental variables - must be PFT level
var_list['pft'].extend([
    'TSOI_ROOTFR', 'H2OSOI_ROOTFR','SMINN_vr_ROOTFR','SOLUTIONP_vr_ROOTFR',
    'SMIN_NH4_vr_ROOTFR', 'SMIN_NO3_vr_ROOTFR'
])
# root fraction weighted organic nutrient availability - must be PFT level
var_list['pft'].extend([
    'LITR1C_vr_ROOTFR', 'LITR2C_vr_ROOTFR', 'LITR3C_vr_ROOTFR',
    'LITR1N_vr_ROOTFR', 'LITR2N_vr_ROOTFR', 'LITR3N_vr_ROOTFR',
    'LITR1P_vr_ROOTFR', 'LITR2P_vr_ROOTFR', 'LITR3P_vr_ROOTFR'
])
if npcompet:
    # nutrient limitation quantities
    var_list['pft'].extend(['FPG_PATCH', 'FPG_P_PATCH',
                            'PLANT_NDEMAND_POT','PLANT_PDEMAND_POT',
                            'FROOT_NDEMAND_POT','FROOT_PDEMAND_POT',
                            'FUNGI_NDEMAND_POT','FUNGI_PDEMAND_POT',
                            'FUNGI_INHIB_PATCH','ZWT_FROOT_PATCH',
                            'FUNGI_SOM_TO_NPOOL','FUNGI_SOM_TO_PPOOL',
                            'CPOOL_TO_FUNGI'])
    # nutrient uptake multipliers
    var_list['pft'].extend([
            'FFR_SRA_PATCH_ROOTFR', 'FFR_N_PATCH_ROOTFR', 'FFR_P_PATCH_ROOTFR', 
            'FFR_TSOI_PATCH_ROOTFR', 'FFR_SWC_PATCH_ROOTFR', 'FFN_NSC_PATCH', 
            'FFN_N_PATCH_ROOTFR', 'FFN_P_PATCH_ROOTFR', 'FFR_FPG_PATCH',
            'FFR_FPG_P_PATCH'])

if rootpheno:
    var_list['pft'] = var_list['pft'] + [
        'ONSET_FLAG', 'OFFSET_FLAG', 'DORMANT_FLAG',
        'BGLFR_LEAF', 'BGLFR_FROOT', 'ONSET_FLAG_ROOT',
        'OFFSET_FLAG_ROOT', 'DORMANT_FLAG_ROOT',
        'FCUR_DYN', 'ONSET_FROOT_FNMIN', 'ONSET_FROOT_FW',
        'LFR_FROOT_TD', 'LFR_FROOT_WD']

var_list['col'] = [
    'TBOT', 'TSOI_30', 'H2OSOI_30', 'HR', 'NEE', 'FPG', 'FPG_P', 'ZWT', 
    'ACTUAL_IMMOB', 'ACTUAL_IMMOB_P', 'SMINN_30', 'SOLUTIONP_30', 
    'FPI', 'FPI_P', 'RH2M', 'H2OSFC']
if not npcompet:
    var_list['col'].extend(['FPG', 'FPG_P'])
var_list['const'] = []

collection_ts, collection_const = extract_sims(prefix, var_list, ensemble_id)

path_out_2 = os.path.join(path_out, 'extract', prefix)
if not os.path.exists(path_out_2):
    os.mkdir(path_out_2)
print(path_out_2)
collection_ts.to_csv(os.path.join(path_out_2, 'analysis_ts.csv'))
collection_const.to_csv(os.path.join(path_out_2, 'analysis_const.csv'))
