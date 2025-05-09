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
from utils.analysis import vert_interp, get_dissolved_nutrients

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

workdir = os.getcwd()

# Down-select the runs because it is impossible to finish processing
# all 4000 runs in reasonable amount of time
PREFIX = 'UQ_20240107'
if PREFIX == 'UQ_20231118':
    subset_ensemble =  [1903,823,2598,2028,3796,3165,2526,2691,3453,1534,3108,2472,231,1196,
                        907,3620,1181,262,2848,3613,3067,2562,1125,1757,1679,1682,2154,98,
                        1103,3003,694,3661,652,709,3711,255,1660,3767,3564,1130]
elif PREFIX == 'UQ_20240107':
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

chambers_ordered = {
    'amb': [6, 20, 13, 8, 17], 
    'elev': [19, 11, 4, 16, 10]
}

# Read the observed nutrient level because I need to get the dates
DEPTH = 0.3
obs_data = get_dissolved_nutrients(DEPTH)

# Function to perform post-processing for one ensemble member
def postproc(ind):
    tvec = pd.date_range("2015-01-01", "2021-12-31", freq="1D")
    tvec = tvec[(tvec.month != 2) | (tvec.day != 29)]
    rid = subset_ensemble[ind]

    LEVGRND = np.array([0.007100635, 0.027925, 0.06225858, 0.1188651, 0.2121934,
                        0.3660658, 0.6197585, 1.038027, 1.727635, 2.864607, 4.739157,
                        7.829766, 12.92532, 21.32647, 35.17762])

    var_list = ['SMINN_vr','SMIN_NH4_vr','SMIN_NO3_vr','SOLUTIONP_vr']
    obs_list = ['NH4','NO3','SRP']


    temp = pd.DataFrame(
        np.nan, index = tvec,
        columns = pd.MultiIndex.from_product([chamber_list_complete, var_list],
                                             names = ['plot','variable'])
    )

    for plot in chamber_list_complete:
        path_data = os.path.join(os.environ["E3SM_ROOT"], "output", "UQ", 
                                    f"{PREFIX}_US-SPR_ICB20TRCNPRDCTCBC", f"g{rid:05g}",
                                    chamber_list_complete_dict[f'P{plot:02d}'])
        for year in range(2015, 2022):
            # print(plot, year)
            filename = os.path.join(path_data, f'{PREFIX}_US-SPR_ICB20TRCNPRDCTCBC.elm.h1.' + 
                                    f'{year}-01-01-00000.nc')
            if not os.path.exists(filename):
                return None
            nc = Dataset(filename)
            for var in var_list:

                target_nodes = [DEPTH]
                input_nodes = LEVGRND
                target_single_level = True

                # hollow only
                col = 1
                input_data = nc[var][:,:,col] / nc['H2OSOI'][:,:,col]
                result = vert_interp(target_nodes, input_nodes, input_data, 
                                     target_single_level)
                temp.loc[temp.index.year == year, (plot, var)] = result[:, 0]
            nc.close()

    # For this simulation, subset to the observed dates, calculate temporal trend and average values
    collect_mean = pd.DataFrame(np.nan,
                                index = chambers_ordered['amb'] + chambers_ordered['elev'],
                                columns = ['NH4','NO3','SRP'])
    collect_trend = pd.DataFrame(np.nan,
                                 index = chambers_ordered['amb'] + chambers_ordered['elev'],
                                 columns = ['NH4','NO3','SRP'])
    collect_trend_p = collect_trend.copy()

    for sim_name, obs_name in zip(var_list, obs_list):

        # get the days since begin from observed data
        obs_data_var = obs_data[['PLOT', 'DATE', obs_name]
                                ].set_index(['DATE', 'PLOT']).iloc[:,0].dropna().unstack()
        day_since_begin = (obs_data_var.index - obs_data_var.index[0]).days

        for col in chambers_ordered['amb'] + chambers_ordered['elev']:

            filt = ~np.isnan(obs_data_var[f'{col:02d}'].values)

            # simulated data on the same dates as observed
            sims_data_matched = pd.DataFrame(np.nan, index = obs_data_var.index, 
                                             columns = obs_data_var.columns)
            for plot in chambers_ordered['amb'] + chambers_ordered['elev']:
                sims_data_matched.loc[:, plot] = temp.loc[obs_data_var.index, (plot,sim_name)]

            for col in chambers_ordered['amb'] + chambers_ordered['elev']:
                filt = ~np.isnan(obs_data_var[f'{col:02d}'].values)

                collect_mean.loc[col,obs_name] = np.mean(sims_data_matched[col].values[filt])

                res = linregress(day_since_begin[filt], sims_data_matched[col].values[filt])
                collect_trend.loc[col, obs_name] = res.slope
                collect_trend_p.loc[col, obs_name] = res.pvalue

    collection = np.stack([collect_mean.values, collect_trend.values, 
                                 collect_trend_p.values], axis=2)

    return collection

# Debug
#collection = postproc(1)

"""
rank = 0
b = 2
collection = np.full([BLOCK, len(chambers_ordered['amb']+chambers_ordered['elev']) , 3, 3], np.nan)
for ind in range(BLOCK):
    print(f'Start index to ensemble = {ind + b*BLOCK}, block = {b}, rank = {rank}')
    sys.stdout.flush()
    result = postproc(ind + b*BLOCK)
    if not result is None:
        collection[ind, :, :, :] = result
    else:
        ierr = ind
collection.dump(
    os.path.join(
        path_out, 'extract', PREFIX, f'ensemble_dissolved_nutrients_part{b:03g}.bin'
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
            collection = np.full([BLOCK, len(chambers_ordered['amb']+chambers_ordered['elev']),
                                  3, 3], np.nan)
            for ind in range(BLOCK):
                print(f'Start ensemble = {ind + b*BLOCK + 1}, block = {b}, rank = {rank}')
                sys.stdout.flush()
                result = postproc(ind + b*BLOCK)
                if not result is None:
                    collection[ind, :, :, :] = result
                else:
                    ierr = ind

            collection.dump(
                os.path.join(
                    path_out, 'extract', PREFIX, f'ensemble_dissolved_nutrients_part{b:03g}.bin'
                )
            )
            comm.send(rank, dest=0, tag=3)
            comm.send(ierr, dest=0, tag=4)
            comm.send(b, dest=0, tag=5)

MPI.Finalize()
