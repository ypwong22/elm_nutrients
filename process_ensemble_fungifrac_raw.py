# Use netCDF4 because xarray is too slow
# Save the whole matrix instead of just mean and slope
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
PREFIX = 'UQ_20240112'
if PREFIX == 'UQ_20240107':
    subset_ensemble = [249,2702,61,936,174,1535,1485,1002,3119,1741,2263,1081,1590,103,
                       1918,2786,2541,1532,559,1685,902,1026,1392,1881,1377,316,722,1069,
                       1957,131,1198,1027,616,2120,1963,3384,3510,1847,653,3203]
elif PREFIX == 'UQ_20240112':
    subset_ensemble = [1944,2204,1429,2540,2764,867,654,1870,1079,1362,2922,1330,448,2364,
                       2601,2165,2421,3609,3554,34,1172,1480,161,1980,2259,840,3116,1071,
                       3744,1372,2495,2096,3204,55,3726,3297,2971,1838,1002,3650]
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
    'amb': [7, 6, 20, 13, 8, 17], 
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

    # (chamber) x (mean average)
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

    collection = nu_uptake_fraction.values
    return collection

## Debug
##collection = postproc(1)


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
            collection = np.full([BLOCK, 11, 18], np.nan)
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
                    path_out, 'extract', PREFIX, f'ensemble_fungifrac_raw_part{b:03g}.bin'
                )
            )
            comm.send(rank, dest=0, tag=3)
            comm.send(ierr, dest=0, tag=4)
            comm.send(b, dest=0, tag=5)

MPI.Finalize()
