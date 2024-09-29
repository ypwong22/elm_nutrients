# Use netCDF4 because xarray is too slow
from netCDF4 import Dataset
from mpi4py import MPI
import itertools as it
import xarray as xr
import pandas as pd
import numpy as np
import os
from scipy.stats import t,linregress
from glob import glob
from utils.constants import *
from utils.paths import *
import sys

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

workdir = os.getcwd()

PREFIX = 'UQ_20240311_1'
if PREFIX == 'UQ_20240311_1':
    # Down-select the runs because it is impossible to finish processing
    # all 4000 runs in reasonable amount of time
    subset_ensemble = \
        [2472, 1476, 2954, 1082, 2537, 2879, 1741, 3207, 2788,  848,  598, 2055, 3426, 2552, 
         3346, 3144, 3550, 2339,  760, 2337,  638, 1198,  883, 2047, 2728, 3211, 3563, 2722, 
         2721, 2558,  743,  318, 2364, 3412, 1790,  769, 2717, 1568, 1409, 2714, 3000, 3512, 
         3333,   44, 3576, 1171,  323,  974, 2900, 1304, 1599,   52, 1397,  547, 3824, 2309, 
           57, 3503, 2708,  542, 3947, 3944, 1811, 1564, 3275, 1165, 3825, 2396, 2706, 2568, 
          826,  534, 2401, 1390,   75,   76, 1679, 1223, 1381,  526,  523, 1154,  663,   84, 
         3923, 2947, 3268, 3920, 3265,  331,   91, 2498, 2685, 1852, 3158, 1095, 2955,  667, 
         1563, 2267, 1375, 3760,  501, 3907,  498, 1535,  813, 2820, 3639, 1368, 1367, 2575, 
         3085, 2868, 2914, 2249,  228, 3896, 3327,  488, 2864, 2969,  905,  674,  225, 2590, 
         3650, 2594,  129,  808, 2663, 2221, 1994, 3063, 1278, 1059, 2596,  368, 1888, 1009, 
         1110,  370,  143, 1979, 1017,  688,  804,  447, 2651, 3841, 1916, 1637,  207, 2643, 
          435, 1057, 2205, 1336,  205, 1628, 1274, 3321, 2851, 3325, 1113, 3844, 1951,  168, 
         1548, 3392, 1953, 1955, 3456,  174,  175, 1321, 2194, 3714, 2187,  398, 1051,  394, 
         2905, 3223, 1464, 2174,  193, 2172, 1115, 2532,  632,  858, 3559, 1056, 1144,  932, 
          217, 2711, 1812, 2907]

elif PREFIX == 'UQ_20240311_2':
    subset_ensemble = [1312,981,3517,3473,1682,360,3452,432,461,3403,11,3393,343,3375,
                       1649,1645,17,3553,502,3332,541,3276,3275,316,582,603,3218,1605,
                       1774,3607,3958,640,3157,3103,35,1578,708,3054,3636,749,3650,1554,
                       2964,1857,2944,253,3656,801,1870,2924,820,1518,2901,845,2888,1912,
                       855,1495,2863,2835,2793,2786,1966,1465,1460,211,1451,931,935,938,
                       940,2725,2011,1425,1421,989,2688,1013,2606,3788,1088,2525,2521,2462,
                       2417,136,2355,133,89,2348,2269,1246,3844,1258,1293,2118,2175,2203,
                       3509,1228,2788,2776,1998,48,3733,3032,1289,108,709,990,4,3772,
                       3194,2548,2053,3306,3523,2482,2418,2386,3416,1694,221,3132,23,1651,
                       3718,2820,492,501,1639,3818,2074,2072,3811,874,1917,847,1302,2768,
                       2391,2889,841,3663,3729,2106,3454,959,1731,150,754,413,3314,3472,
                       2496,2507,1839,1627,1622,2713,604,2534,720,3752,1075,1762,1561,1679,
                       3381,984,2033,3494,1018,291,1833,3427,1975,86,2396,3348,411,2851,
                       525,329,1262,1853,2679,1979,1255,2363,3364,3318,3496,491,292,1208,
                       2750,1003,1166,8]

N = len(subset_ensemble)

# Assume folder already exists
#time.sleep(0.02*rank) # ensure the mkdir doesn't conflict with each other
#if not os.path.exists(os.path.join(path_out, 'extract', PREFIX)):
#    os.mkdir(os.path.join(path_out, 'extract', PREFIX))

# number of ensembles to save in each bin file
# this avoids having difficulty in dumping file
BLOCK = 4
if np.mod(N, BLOCK) != 0:
    raise Exception("N must be a multiply of BLOCK")
niter = int(N / BLOCK)

var_list_pft = ['FPG_PATCH', 'FPG_P_PATCH', 'PLANT_NDEMAND_POT',
                'PLANT_PDEMAND_POT', 'FROOT_NDEMAND_POT','FROOT_PDEMAND_POT',
                'FUNGI_NDEMAND_POT','FUNGI_PDEMAND_POT',
                'SMINN_TO_NPOOL','SMINP_TO_PPOOL',
                'FUNGI_SOM_TO_NPOOL','FUNGI_SOM_TO_PPOOL', 'CPOOL_TO_FUNGI']
chambers_ordered = {
    'amb': [6, 20, 13, 8, 17], 
    'elev': [19, 11, 4, 16, 10]
}

# Function to perform post-processing for one ensemble member
def postproc(ind):
    pft_list = np.array([2, 3, 11])
    hol_add = 17
    tvec = pd.date_range("2015-01-01", "2021-12-31", freq="1D")
    tvec = tvec[(tvec.month != 2) | (tvec.day != 29)]
    rid = subset_ensemble[ind]

    temp = pd.DataFrame(
        np.nan, index = tvec,
        columns = pd.MultiIndex.from_product([
            chamber_list_complete, var_list_pft, pft_list
        ], names = ['plot','variable','pft'])
    )

    for plot in chamber_list_complete:
        path_data = os.path.join(os.environ["E3SM_ROOT"], "output", "UQ", 
                                    f"{PREFIX}_US-SPR_ICB20TRCNPRDCTCBC", f"g{rid:05g}",
                                    chamber_list_complete_dict[f'P{plot:02d}'])
        for year in range(2015, 2022):
            # print(plot, year)
            filename = os.path.join(path_data, f'{PREFIX}_US-SPR_ICB20TRCNPRDCTCBC.elm.h2.' + 
                                    f'{year}-01-01-00000.nc')
            if not os.path.exists(filename):
                return None
            nc = Dataset(filename)
            for var in var_list_pft:
                # DEBUG: print(plot, year, var)
                # average hummock & hollow
                temp.loc[temp.index.year == year, (plot, var)] = \
                    nc[var][:, pft_list].data * 0.64 + nc[var][:, pft_list+hol_add].data * 0.36
            nc.close()

    # average hummock and hollow, resample to growing season's all-year average
    # unit: g m-2 day-1
    filt = (tvec.month >= 5) & (tvec.month <= 10)
    temp = temp.loc[filt, :].mean().unstack().unstack( \
        ).loc[chambers_ordered['amb'] + chambers_ordered['elev'], :]

    nu_uptake_fraction = pd.DataFrame(np.nan,
        index = chambers_ordered['amb'] + chambers_ordered['elev'],
        columns = pd.MultiIndex.from_product([['N','P'],
                                              ['Spruce','Tamarack','Shrub'], 
                                              ['min froot','min fungi','org fungi']]))
    for nu in ['N','P']:
        for i, (pft, name) in enumerate(zip([2,3,11],['Spruce','Tamarack','Shrub'])):
            froot_min = temp.loc[:, (pft, f'SMIN{nu}_TO_{nu}POOL')].values / \
                temp.loc[:, (pft, f'PLANT_{nu}DEMAND_POT')].values * \
                temp.loc[:, (pft, f'FROOT_{nu}DEMAND_POT')].values
            fungi_min = temp.loc[:, (pft, f'SMIN{nu}_TO_{nu}POOL')].values / \
                temp.loc[:, (pft, f'PLANT_{nu}DEMAND_POT')].values * \
                temp.loc[:, (pft, f'FUNGI_{nu}DEMAND_POT')].values
            fungi_som = temp.loc[:, (pft, f'FUNGI_SOM_TO_{nu}POOL')].values
            retemp = np.vstack([froot_min, fungi_min, fungi_som]).T * 86400
            nu_uptake_fraction.loc[:, (nu,name)] = retemp

    # ('amb', 'elev') x (mean, mean_std, slope, slope_std)
    # ambient -> elevated CO2
    nu_uptake_fraction_summary = pd.DataFrame(np.nan,
        index = ['mean', 'slope'],
        columns = pd.MultiIndex.from_product([['amb', 'elev'], ['N','P'],
                                              ['Spruce','Tamarack','Shrub'],
                                              ['min froot','min fungi','org fungi']]))
    for i, co2 in enumerate(['amb', 'elev']):
        temp = nu_uptake_fraction.loc[chambers_ordered[co2], :]

        # average of all the chambers
        for col in nu_uptake_fraction.columns:
            nu_uptake_fraction_summary.loc['mean', (co2, *col)] = np.mean(temp[col].values)
            res = linregress(np.arange(len(chambers_ordered[co2])), temp[col].values)
            nu_uptake_fraction_summary.loc['slope', (co2, *col)] = res.slope

    collection = nu_uptake_fraction_summary.values
    return collection

# Debug
#collection = postproc(1)

"""
rank = 0
b = 49 
collection = np.full([BLOCK, 2, 36], np.nan)
for ind in range(BLOCK):
    print(f'Start index to ensemble = {ind + b*BLOCK}, block = {b}, rank = {rank}')
    sys.stdout.flush()
    result = postproc(ind + b*BLOCK)
    if not result is None:
        collection[ind, :, :] = result
    else:
        ierr = ind
collection.dump(
    os.path.join(
        path_out, 'extract', PREFIX, f'ensemble_fungifrac_part{b:03g}.bin'
    )
)
"""

if rank == 0:
    n_done = 0

    # send individual jobs
    for b in range(niter):
        # note: the processors are labeled from 0 to size-1
        if b < (size-1):
            comm.send(0, dest=b+1, tag=1)
            comm.send(b, dest=b+1, tag=2)
        else:
            # assign rest of jobs on demand
            process = comm.recv(source=MPI.ANY_SOURCE, tag=3)
            ierr = comm.recv(source=process, tag=4)
            b_id = comm.recv(source=process, tag=5)
            if ierr > -1:
                print(f'Rank = {process}, block = {b_id}, ensemble member = ' + \
                    f'{ierr + b_id*BLOCK} cannot be processed.')
            n_done = n_done + 1
            comm.send(0, dest=process, tag=1)
            comm.send(b, dest=process, tag=2)

    # receive remaining messages and finalize
    while n_done < niter:
        process = comm.recv(source=MPI.ANY_SOURCE, tag=3)
        ierr = comm.recv(source=process, tag=4)
        b_id = comm.recv(source=process, tag=5)
        if ierr > -1:
            print(f'Rank = {process}, block = {b_id}, ensemble member = ' + \
                  f'{ierr + b_id*BLOCK} cannot be processed.')
        n_done = n_done + 1
        comm.send(-1, dest=process, tag=1)
        comm.send(-1, dest=process, tag=2)

# --------------------- Slave process (get data and calculate) --------------
else:
    status = 0
    while status == 0:
        status = comm.recv(source=0, tag=1)
        b = comm.recv(source=0, tag = 2)
        ierr = -1

        if status == 0:
            collection = np.full([BLOCK, 2, 36], np.nan)
            for ind in range(BLOCK):
                print(f'Start ensemble = {ind + b*BLOCK + 1}, block = {b}, rank = {rank}')
                sys.stdout.flush()
                result = postproc(ind + b*BLOCK)
                if not result is None:
                    collection[ind, :, :] = result
                else:
                    ierr = ind

            collection.dump(
                os.path.join(
                    path_out, 'extract', PREFIX, f'ensemble_fungifrac_part{b:03g}.bin'
                )
            )
            comm.send(rank, dest=0, tag=3)
            comm.send(ierr, dest=0, tag=4)
            comm.send(b, dest=0, tag=5)

MPI.Finalize()
