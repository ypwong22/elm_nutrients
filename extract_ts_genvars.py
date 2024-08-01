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

prefix = '20240311_3_1'
rootpheno = False
npcompet = True

#for prefix,rootpheno,npcompet in zip(['20231113', '20240311', '20240316'],
#                                     [False, False, True],  [False, True, True]):

var_list = {}
var_list['pft'] = ['TLAI','GPP','NPP','AGNPP','BGNPP','FROOTC_ALLOC','LEAFC_ALLOC','TOTVEGC',
                   'TOTVEGC_ABG','MR','GR','XR','LITFALL','FROOTN','FROOTP','LEAFN','LEAFP',
                   'SMINN_TO_NPOOL', 'SMINP_TO_PPOOL', 'CPOOL', 'NPOOL', 'PPOOL', 'BTRAN',
                   'PLANT_NDEMAND','PLANT_PDEMAND', 'RETRANSN_TO_NPOOL','RETRANSP_TO_PPOOL',
                   'XSMRPOOL']
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
    'TSOI_ROOTFR', 'H2OSOI_ROOTFR','SMINN_vr_ROOTFR','SOLUTIONP_vr_ROOTFR'
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
    # root fraction weighted organic nutrient availability - must be PFT level
    var_list['pft'].extend([
        'LITR1N_vr_ROOTFR', 'LITR2N_vr_ROOTFR', 'LITR3N_vr_ROOTFR',
        'LITR1P_vr_ROOTFR', 'LITR2P_vr_ROOTFR', 'LITR3P_vr_ROOTFR'
    ])

if rootpheno:
    var_list['pft'] = var_list['pft'] + [
        'ONSET_FLAG', 'OFFSET_FLAG', 'DORMANT_FLAG',
        'BGLFR_LEAF', 'BGLFR_FROOT', 'ONSET_FLAG_ROOT',
        'OFFSET_FLAG_ROOT', 'DORMANT_FLAG_ROOT',
        'FCUR_DYN', 'ONSET_FROOT_FNMIN', 'ONSET_FROOT_FW',
        'LFR_FROOT_TD', 'LFR_FROOT_WD']

var_list['col'] = [
    'TBOT', 'TSOI_30', 'H2OSOI_30', 'HR', 'NEE', 'FPG', 'FPG_P', 'ZWT', 
    'ACTUAL_IMMOB', 'ACTUAL_IMMOB_P', # 'SMINN_30', 'SOLUTIONP_30', 
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



########################################################################################
# root fraction weighted layerwise variables
########################################################################################

def extract_sims_rootfr_weighted(prefix, var_list={"pft": [], "col": []}):
    tvec = pd.date_range("2015-01-01", "2021-12-31", freq="1D")
    tvec = tvec[(tvec.month != 2) | (tvec.day != 29)]

    pft_list = [2, 3, 11, 12]
    hol_add = 17

    collection_ts = {}
    collection_const = pd.DataFrame(
        np.nan, index=["hummock", "hollow"], columns=var_list["const"]
    )

    for plot in chamber_list_complete:
        path_data = os.path.join(os.environ["PROJDIR"], "E3SM",
            "output", f"{prefix}_US-SPR_ICB20TRCNPRDCTCBC",
            "spruce_treatments", f'plot{plot:02d}_US-SPR_ICB20TRCNPRDCTCBC', 'run')

        print(f'Plot {plot}', path_data)

        flist_pft = sorted(glob(os.path.join(path_data, "*.h2.*.nc")))[:-1]
        flist_col = sorted(glob(os.path.join(path_data, "*.h1.*.nc")))[:-1]

        var_list_pft = var_list["pft"]
        hr = xr.open_mfdataset(flist_pft, decode_times=False)
        for var in var_list_pft:
            for pft in pft_list:
                if var == 'TOTVEGC':
                    # need to subtract out CPOOL to be comparable to obs
                    collection_ts[(plot, var, pft, "hummock")] = hr[var][:, pft].values - \
                        hr['CPOOL'][:, pft].values
                    collection_ts[(plot, var, pft, "hollow")] = \
                        hr[var][:, pft + hol_add].values - \
                        hr['CPOOL'][:, pft + hol_add].values
                elif var == "TSOI_AVG":
                    rootfr = hr["ROOTFR"][:, :, pft].values
                    hr2 = xr.open_mfdataset(flist_col)
                    collection_ts[(plot, var, pft, "hummock")] = np.sum(
                        (hr2["TSOI"][:, :, 0].values - 273.15) * rootfr, axis=1
                    )
                    collection_ts[(plot, var, pft, "hollow")] = np.sum(
                        (hr2["TSOI"][:, :, 1].values - 273.15) * rootfr, axis=1
                    )
                    hr2.close()

                elif var == "SWC_AVG":
                    rootfr = hr["ROOTFR"][:, :, pft].values
                    hr2 = xr.open_mfdataset(flist_col)
                    collection_ts[(plot, var, pft, "hummock")] = np.sum(
                        hr2["H2OSOI"][:, :, 0].values * rootfr, axis=1)
                    collection_ts[(plot, var, pft, "hollow")] = np.sum(
                        hr2["H2OSOI"][:, :, 1].values * rootfr, axis=1)
                    hr2.close()

                else:
                    collection_ts[(plot, var, pft, "hummock")] = hr[var][:, pft].values
                    collection_ts[(plot, var, pft, "hollow")] = hr[var][:, pft + hol_add].values
        hr.close()

        var_list_col = var_list["col"]
        hr = xr.open_mfdataset(flist_col, decode_times=False)
        for var in var_list_col:
            if 'TSOI_' in var or 'H2OSOI_' in var or 'SMINN_' in var or 'SOLUTIONP_' in var:

                LEVGRND = np.array([0.007100635, 0.027925, 0.06225858, 0.1188651, 0.2121934,
                                    0.3660658, 0.6197585, 1.038027, 1.727635, 2.864607, 4.739157,
                                    7.829766, 12.92532, 21.32647, 35.17762])
                LEVGRND_I = np.append(np.insert(
                    (LEVGRND[1:] + LEVGRND[:-1])*0.5, 0, 0
                ), LEVGRND[-1] + 0.5 * (LEVGRND[-1] - LEVGRND[-2]))
                THICKNESS = np.diff(LEVGRND_I)

                depth = float(var.split('_')[1]) / 100.
                maxlayer = np.where(LEVGRND_I < depth)[0][-1]

                if 'TSOI_' in var:
                    thisvar = 'TSOI'
                elif 'H2OSOI_' in var:
                    thisvar = 'H2OSOI'
                elif 'SMINN_' in var:
                    thisvar = 'SMINN_vr'
                elif 'SOLUTIONP_' in var:
                    thisvar = 'SOLUTIONP_vr'

                data = np.zeros([hr[thisvar].shape[0], 2])
                for i in range(maxlayer):
                    data = data + hr[thisvar][:, i, :] * THICKNESS[i]
                last_depth = min(THICKNESS[i], depth - LEVGRND_I[maxlayer])
                data = data + hr[thisvar][:, maxlayer, :] * last_depth
                data = data / depth

                collection_ts[(plot, var, 0, "hummock")] = data.values[:, 0]
                collection_ts[(plot, var, 0, "hollow")] = data.values[:, 1]

                if 'TSOI' in var:
                    collection_ts[(plot, var, 0, "hummock")] -= 273.15
                    collection_ts[(plot, var, 0, "hollow")] -= 273.15

            elif var == "SMP_MAX":
                collection_ts[(plot, var, 0, "hummock")] = np.nanmax(
                    hr["SMP"][:, :, 0].values, axis=1
                )  # mm
                collection_ts[(plot, var, 0, "hollow")] = np.nanmax(
                    hr["SMP"][:, :, 1].values, axis=1
                )

            elif var == "ZWT":
                hr2 = xr.open_mfdataset(flist_col)
                collection_ts[(plot, var, 0, "hummock")] = (
                    0.15 - hr2["ZWT"][:, 0].values
                )
                collection_ts[(plot, var, 0, "hollow")] = (
                    hr2["H2OSFC"][:, 1] / 1000.0 - hr2["ZWT"][:, 1]
                )
                hr2.close()

            else:
                collection_ts[(plot, var, 0, "hummock")] = hr[var][:, 0].values
                collection_ts[(plot, var, 0, "hollow")] = hr[var][:, 1].values
                if var == 'TBOT':
                    collection_ts[(plot, var, 0, "hummock")] -= 273.15
                    collection_ts[(plot, var, 0, "hollow")] -= 273.15
        hr.close()

        # some time-constant variables are in the h0 files
        var_list_const = var_list["const"]
        hr = xr.open_mfdataset(flist_const, decode_times=False)
        for var in var_list_const:
            if var == "SUCSAT":
                collection_const.loc[:, "SUCSAT"] = pd.Series(
                    hr["SUCSAT"][2, :].values, index=["hummock", "hollow"]
                )
            elif var == "WATSAT":
                collection_const.loc[:, "WATSAT"] = pd.Series(
                    hr["WATSAT"][2, :].values, index=["hummock", "hollow"]
                )
            elif var == "BSW":
                collection_const.loc[:, "BSW"] = pd.Series(
                    hr["BSW"][2, :].values, index=["hummock", "hollow"]
                )
            else:
                raise Exception("Not implemented")
        hr.close()

    collection_ts = pd.DataFrame(collection_ts)
    collection_ts.index = tvec
    collection_ts.columns.names = ["plot", "variable", "pft", "topo"]

    return collection_ts, collection_const
