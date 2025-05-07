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

# Down-select the runs because it is impossible to finish processing
# all 4000 runs in reasonable amount of time
PREFIX = 'UQ_20240107'
if PREFIX == 'UQ_20240101':
    subset_ensemble = [2656,147,3470,124,3454,1223,2977,1787,204,3259,877,1322,
                       649,1846,3001,2390,226,3545,3767,3380]
elif PREFIX == 'UQ_20240102':
    subset_ensemble = [3941,3377,1899,1876,744,3690,2348,3681,1810,3641,776,1257,1089,1921,
                       1929,1933,3263,1967,1693,3305]
elif PREFIX == 'UQ_20240107':
    subset_ensemble = [276,2578,3203,2263,2366,1002,174,2583,605,1458,3420,873,1485,2786,
                       140,1931,2030,1741,1535,3384,103,2541,2702,936,1881,1876,722,715,
                       64,653,1081,3119,1847,249,35,2316,49,1917,61,2100]
elif PREFIX == 'UQ_20240107':
    subset_ensemble = [1870,1944,1002,1330,1071,448,654,3833,161,2601,2585,867,2764,1079,
                       3629,1372,2364,3116,2675,1838,2421,1172,3609,3204,55,2204,1221,
                       34,1480,2540,3558,1542,3726,2234,2422,2971,2165,840,3490,2495]
N = len(subset_ensemble)

# Assume folder already exists
#time.sleep(0.02*rank) # ensure the mkdir doesn't conflict with each other
#if not os.path.exists(os.path.join(path_out, 'extract', PREFIX)):
#    os.mkdir(os.path.join(path_out, 'extract', PREFIX))

# number of ensembles to save in each bin file
# this avoids having difficulty in dumping file
BLOCK = 10
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
