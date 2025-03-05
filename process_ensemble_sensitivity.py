#!/usr/bin/env python3
import sys, os, time
import numpy as np
from scipy.stats import linregress, t
from mpi4py import MPI
from utils.paths import *
from utils.constants import chamber_list_names_complete
import pandas as pd
import time
from utils.analysis import get_sim_carbonfluxes

# delete when debug
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

workdir = os.getcwd()

# number of simulations
N = 4000
#N = 2000
#N = 3125
#N = 1850

PREFIX = "UQ_20240102"
#time.sleep(0.02*rank) # ensure the mkdir doesn't conflict with each other
if not os.path.exists(os.path.join(path_out, 'extract', PREFIX)):
    os.mkdir(os.path.join(path_out, 'extract', PREFIX))

# number of ensembles to save in each bin file
# this avoids having difficulty in dumping file
BLOCK = 200
#BLOCK = 125
#BLOCK = 99
#BLOCK = 50
if np.mod(N, BLOCK) != 0:
    raise Exception("N must be a multiply of BLOCK")

RUNROOT = os.path.join(os.environ["E3SM_ROOT"], "output")
VAR_LIST = ['Tair', 'AGBiomass_Spruce', 'AGBiomass_Tamarack', 'AGBiomass_Shrub',
            'AGNPPtoBiomass_Spruce', 'AGNPPtoBiomass_Tamarack', 'AGNPPtoBiomass_Shrub',
            'AGNPP_Spruce', 'AGNPP_Tamarack', 'AGNPP_Shrub', 'NPP_moss',
            'BGNPP_TreeShrub', 'BGtoAG_TreeShrub', 'NPP', 'HR', 'NEE']
YEAR_LIST = range(2015, 2022)  # skip 2015 and 2021 because no observation/no model data

MOSSFRAC = pd.read_excel("Sphagnum_fraction.xlsx", index_col=0, skiprows=1, engine="openpyxl"
                         ).drop(["plot", "Temp", "CO2"], axis=1)
MOSSFRAC[2015] = MOSSFRAC[2016]
# MOSSFRAC = MOSSFRAC.drop(2021, axis=1)

niter = int(N / BLOCK)


# Function to perform post-processing for one ensemble member
def postproc(thisjob, collection):
    ## Debug
    ##thisjob = 3852
    ##collection = np.empty([len(VAR_LIST), 2, 4])

    ierr = 0

    casename = f"{PREFIX}_US-SPR_ICB20TRCNPRDCTCBC"
    baserundir = os.path.join(RUNROOT, "UQ", casename, f"g{thisjob:05g}")
    # print(baserundir)    
    # print(thisjob)

    try:
         values = get_sim_carbonfluxes(YEAR_LIST, baserundir, False, extra_col_vars=['TOTSOMC'])
    except:
        print(f'unsuccessful reading of {baserundir}')
        ierr = 1
        return ierr

    values = values.loc['average', :]

    # ('amb', 'elev') x (mean, mean_std, slope, slope_std)
    # ambient -> elevated CO2
    chamber_list = {'amb': ['P07', 'P06', 'P20', 'P13', 'P08', 'P17'], 
                    'elev': ['P19', 'P11', 'P04', 'P16', 'P10']}
    for i, co2 in enumerate(['amb', 'elev']):
        temp = values.index.get_level_values(0).isin(chamber_list[co2])

        # average of all the chambers
        for k, var in enumerate(VAR_LIST):
            collection[k, i, 0] = np.mean(values.loc[temp, var].values)
            collection[k, i, 1] = np.std(values.loc[temp, var].values)

        ts = abs(t.ppf(0.05, len(temp) * len(YEAR_LIST) - 2))
        for k, var in enumerate(VAR_LIST):
            res = linregress(values.loc[temp, 'Tair'].values, values.loc[temp, var].values)
            collection[k, i, 2] = res.slope
            collection[k, i, 3] = ts * res.stderr

    # Debug
    #return collection
    return ierr

# Debug
#collection = postproc(None, None)


for b in range(niter):
    print("rank = ", rank, "b = ", b, flush = True)

    if rank == 0:
        # process a whole block
        start = b * BLOCK
        collection_all = np.full([BLOCK, len(VAR_LIST), 2, 4], np.nan)

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
            collection_all[thisjob - 1, :, :, :] = collection
            n_done = n_done + 1
            comm.send(n_job, dest=process, tag=1)
            comm.send(0, dest=process, tag=2)
            comm.send(start, dest=process, tag=3)

        # receive remaining messages and finalize
        while n_done < BLOCK:
            process = comm.recv(source=MPI.ANY_SOURCE, tag=4)
            thisjob = comm.recv(source=process, tag=5)
            collection = comm.recv(source=process, tag=6)
            collection_all[thisjob - 1, :, :, :] = collection
            n_done = n_done + 1
            comm.send(-1, dest=process, tag=1)
            comm.send(-1, dest=process, tag=2)
            comm.send(-1, dest=process, tag=3)

        collection_all.dump(
            os.path.join(
                path_out, "extract", PREFIX, f"ensemble_collection_part{b:03g}.bin"
            )
        )

    # --------------------- Slave process (get data and calculate) --------------
    else:
        # process individual members inside the block
        collection = np.full([len(VAR_LIST), 2, 4], np.nan)

        status = 0
        while status == 0:
            myjob = comm.recv(source=0, tag=1)
            status = comm.recv(source=0, tag=2)
            start = comm.recv(source=0, tag=3)

            if status == 0:
                collection[:, :, :] = np.nan
                ierr = postproc(start + myjob, collection)
                comm.send(rank, dest=0, tag=4)
                comm.send(myjob, dest=0, tag=5)
                comm.send(collection, dest=0, tag=6)
        #print(collection[:, 1, 2])

MPI.Finalize()
