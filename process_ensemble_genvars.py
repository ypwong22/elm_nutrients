#from mpi4py import MPI
import itertools as it
import xarray as xr
import pandas as pd
import numpy as np
import os
from scipy.stats import t,linregress
from glob import glob
from utils.constants import *
from utils.analysis import *
from utils.paths import *

#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()
#size = comm.Get_size()

workdir = os.getcwd()

N = 4000

PREFIX = "UQ_20240311_2"

# Assume folder already exists
#time.sleep(0.02*rank) # ensure the mkdir doesn't conflict with each other
#if not os.path.exists(os.path.join(path_out, 'extract', PREFIX)):
#    os.mkdir(os.path.join(path_out, 'extract', PREFIX))

# number of ensembles to save in each bin file
# this avoids having difficulty in dumping file
BLOCK = 200
if np.mod(N, BLOCK) != 0:
    raise Exception("N must be a multiply of BLOCK")
niter = int(N / BLOCK)

var_list = {}
var_list['col'] = []
var_list['const'] = []
var_list['pft'] = []
var_list['pft'].extend(['FPG_PATCH', 'FPG_P_PATCH',
                        'PLANT_NDEMAND_POT','PLANT_PDEMAND_POT',
                        'FROOT_NDEMAND_POT','FROOT_PDEMAND_POT',
                        'FUNGI_NDEMAND_POT','FUNGI_PDEMAND_POT',
                        'SMINN_TO_NPOOL','SMINP_TO_PPOOL',
                        'FUNGI_SOM_TO_NPOOL','FUNGI_SOM_TO_PPOOL',
                        'CPOOL_TO_FUNGI'])
chambers_ordered = {
    'amb': ['6', '20', '13', '8', '17'], 
    'elev': ['19', '11', '4', '16', '10']
}

# Function to perform post-processing for one ensemble member
def postproc(thisjob, collection):
    # Debug
    #thisjob = 2462
    #collection = np.empty([2, 36])
    ierr = 0

    try:
        collection_ts, _ = extract_sims(PREFIX, var_list, thisjob)
    except:
        ierr = 1
        return ierr

    # average hummock and hollow, resample to growing season
    temp_data = collection_ts.groupby(['plot','variable','pft'], axis = 1
        ).apply(lambda df: df.iloc[:,0] * 0.64 + df.iloc[:,1] * 0.36)
    filt = (temp_data.index.month >= 5) & (temp_data.index.month <= 10)
    temp_data = temp_data.loc[filt, :]
    temp_data.columns = pd.MultiIndex.from_tuples([(str(c[0]), str(c[1]), str(c[2])) for c in temp_data.columns])
    # average over all the years
    temp = temp_data.mean(axis = 0).unstack().unstack( \
        ).loc[chambers_ordered['amb'] + chambers_ordered['elev'], :]

    nu_uptake_fraction = pd.DataFrame(np.nan,
        index = chambers_ordered['amb'] + chambers_ordered['elev'],
        columns = pd.MultiIndex.from_product([['N','P'],
                                              ['Spruce','Tamarack','Shrub'], 
                                              ['min froot','min fungi','org fungi']]))
    for nu in ['N','P']:
        for i, (pft, name) in enumerate(zip(['2','3','11'],['Spruce','Tamarack','Shrub'])):
            froot_min = temp.loc[:, (pft, f'SMIN{nu}_TO_{nu}POOL')] / \
                temp.loc[:, (pft, f'PLANT_{nu}DEMAND_POT')] * \
                temp.loc[:, (pft, f'FROOT_{nu}DEMAND_POT')]
            fungi_min = temp.loc[:, (pft, f'SMIN{nu}_TO_{nu}POOL')] / \
                temp.loc[:, (pft, f'PLANT_{nu}DEMAND_POT')] * \
                temp.loc[:, (pft, f'FUNGI_{nu}DEMAND_POT')]
            fungi_som = temp.loc[:, (pft, f'FUNGI_SOM_TO_{nu}POOL')]
            retemp = pd.concat([froot_min, fungi_min, fungi_som], axis = 1) * 86400
            nu_uptake_fraction.loc[retemp.index, (nu,name)] = retemp.values

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

    collection[:, :] = nu_uptake_fraction_summary.values

    # Debug
    return collection
    #return ierr

# Debug
collection = postproc(None, None)

for b in range(niter):
    print("rank = ", rank, "b = ", b, flush = True)

    if rank == 0:
        # process a whole block
        start = b * BLOCK
        collection_all = np.full([BLOCK, 2, 36], np.nan)

        # --------------------------Collect and save data---------------------
        n_done = 0

        # send first np-1 jobs where np is number of processes
        for n_job in range(1, size):
            comm.send(n_job, dest=n_job, tag=1)
            comm.send(0, dest=n_job, tag=2)
            comm.send(start, dest=n_job, tag=3)

        # assign rest of jobs on demand
        for n_job in range(size, BLOCK + 1):
            process = comm.recv(source=MPI.ANY_SOURCE, tag=4)
            thisjob = comm.recv(source=process, tag=5)
            collection = comm.recv(source=process, tag=6)
            collection_all[thisjob - 1, :, :] = collection
            n_done = n_done + 1
            comm.send(n_job, dest=process, tag=1)
            comm.send(0, dest=process, tag=2)
            comm.send(start, dest=process, tag=3)

        # receive remaining messages and finalize
        while n_done < BLOCK:
            process = comm.recv(source=MPI.ANY_SOURCE, tag=4)
            thisjob = comm.recv(source=process, tag=5)
            collection = comm.recv(source=process, tag=6)
            collection_all[thisjob - 1, :, :] = collection
            n_done = n_done + 1
            comm.send(-1, dest=process, tag=1)
            comm.send(-1, dest=process, tag=2)
            comm.send(-1, dest=process, tag=3)

        collection_all.dump(
            os.path.join(
                path_out, "extract", PREFIX, f"ensemble_fungifrac_part{b:03g}.bin"
            )
        )

    # --------------------- Slave process (get data and calculate) --------------
    else:
        # process individual members inside the block
        collection = np.full([2, 36], np.nan)

        status = 0
        while status == 0:
            myjob = comm.recv(source=0, tag=1)
            status = comm.recv(source=0, tag=2)
            start = comm.recv(source=0, tag=3)

            if status == 0:
                ierr = postproc(start + myjob, collection)
                comm.send(rank, dest=0, tag=4)
                comm.send(myjob, dest=0, tag=5)
                comm.send(collection, dest=0, tag=6)
        #print(collection[:, 1, 2])

MPI.Finalize()
