#!/usr/bin/env python
# /software/dev_tools/swtree/cs400_centos7.2_pe2016-08/openmpi/3.0.0/centos7.2_gnu5.3.0/bin/mpirun -np 64 python process_ensemble_sensitivity.py
import sys, os, time
import numpy as np
from netCDF4 import Dataset
from scipy.stats import linregress, t
from mpi4py import MPI
from utils.paths import *
from utils.constants import chamber_list_names_complete


N = 1000  # number of simulations
PREFIX = "UQ_default2_optimized"
RUNROOT = os.path.join(os.environ["PROJDIR"], "E3SM", "output")
VAR_LIST = ["ANPPtree", "ANPPshrub", "NPPmoss", "BGNPP", "NPP", "HR", "NEE"]
YEAR_LIST = range(2016, 2022)


# Define function to perform ensemble member post-processing
def postproc(thisjob, collection):
    ierr = 0

    casename = f"{PREFIX}_US-SPR_ICB20TRCNPRDCTCBC"
    baserundir = os.path.join(RUNROOT, "UQ", casename, f"g{thisjob:05g}")
    print(thisjob)

    values = np.empty([len(chamber_list_names_complete), len(YEAR_LIST), len(VAR_LIST)])
    t2m = np.empty([len(chamber_list_names_complete), len(YEAR_LIST)])

    for i, folder in enumerate(chamber_list_names_complete):
        for j, year in enumerate(YEAR_LIST):
            hr = Dataset(
                os.path.join(
                    baserundir,
                    folder,
                    f"{casename}.clm2.h1.{year}-01-01-00000.nc",
                )
            )
            hr2 = Dataset(
                os.path.join(
                    baserundir,
                    folder,
                    f"{casename}.clm2.h0.{year}-01-01-00000.nc",
                )
            )

            for k, varname in enumerate(VAR_LIST):
                try:
                    if varname in ["NPPmoss", "ANPPtree", "ANPPshrub", "BGNPP"]:
                        if varname == "NPPmoss":
                            var = hr["NPP"][:, :].mean(axis=0)
                            values[i, j, k] = (
                                var[12] * 0.64 + var[12 + 17] * 0.36
                            ) * 0.25
                        elif varname == "ANPPtree":
                            var = hr["ANPP"][:, :].mean(axis=0)
                            values[i, j, k] = (
                                var[2] * 0.64 + var[2 + 17] * 0.36
                            ) * 0.36 + (var[3] * 0.64 + var[3 + 17] * 0.36) * 0.14
                        elif varname == "ANPPshrub":
                            var = hr["ANPP"][:, :].mean(axis=0)
                            values[i, j, k] = (
                                var[11] * 0.64 + var[11 + 17] * 0.36
                            ) * 0.25
                        elif varname == "BGNPP":
                            # tree and shrub only
                            var = hr["BGNPP"][:, :].mean(axis=0)
                            values[i, j, k] = (
                                (var[2] * 0.64 + var[2 + 17] * 0.36) * 0.36
                                + (var[3] * 0.64 + var[3 + 17] * 0.36) * 0.14
                                + (var[11] * 0.64 + var[11 + 17] * 0.36) * 0.25
                            )
                    else:
                        var = hr2[varname][:, :].mean(axis=0) - 273.15
                        values[i, j, k] = var[0] * 0.64 + var[1] * 0.36
                except:
                    ierr = 1

            var = hr2["TBOT"][:, :].mean(axis=0) - 273.15
            t2m[i, j] = var[0] * 0.64 + var[1] * 0.36

            hr.close()
            hr2.close()

    # ('amb', 'elev') x (mean, mean_std, slope, slope_std)
    # collection = np.empty([len(VAR_LIST), 2, 4])
    # ambient -> elevated CO2
    for i, temp in enumerate(
        [
            np.insert(np.arange(1, len(chamber_list_names_complete), 2), 0, 0),
            np.arange(2, len(chamber_list_names_complete), 2),
        ]
    ):
        # Note this is the ambient or T0.00CO2 chamber
        for k, _ in enumerate(VAR_LIST):
            collection[k, i, 0] = values[temp[0], :, k]
            collection[k, i, 1] = np.std(values[temp[0], :, k])

        ts = abs(t.ppf(0.05, len(temp) * len(YEAR_LIST) - 2))
        for k, _ in enumerate(VAR_LIST):
            res = linregress(t2m[temp, :].reshape(-1), values[temp, :, k].reshape(-1))
            collection[k, i, 2] = res.slope
            collection[k, i, 3] = ts * res.stderr

    return ierr


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

print(rank)

workdir = os.getcwd()

niter = 1

if rank == 0:
    collection_all = np.empty([N, len(VAR_LIST), 2, 4], float)

    # --------------------------Collect and save data---------------------
    for thisiter in range(0, niter):
        n_done = 0

        # send first np-1 jobs where np is number of processes
        for n_job in range(1, size):
            comm.send(n_job, dest=n_job, tag=1)
            comm.send(0, dest=n_job, tag=2)

        # assign rest of jobs on demand
        for n_job in range(size, N + 1):
            process = comm.recv(source=MPI.ANY_SOURCE, tag=3)
            thisjob = comm.recv(source=process, tag=4)
            collection = comm.recv(source=process, tag=5)
            collection_all[thisjob - 1, :, :, :] = collection
            n_done = n_done + 1
            comm.send(n_job, dest=process, tag=1)
            comm.send(0, dest=process, tag=2)

        # receive remaining messages and finalize
        while n_done < N:
            process = comm.recv(source=MPI.ANY_SOURCE, tag=3)
            thisjob = comm.recv(source=process, tag=4)
            collection = comm.recv(source=process, tag=5)
            collection_all[thisjob - 1, :, :, :] = collection
            n_done = n_done + 1
            comm.send(-1, dest=process, tag=1)
            comm.send(-1, dest=process, tag=2)

    collection_all.dump(
        os.path.join(path_out, "extract", PREFIX, f"ensemble_collection.bin"),
    )

    MPI.Finalize()

# --------------------- Slave process (get data and calculate) --------------
else:
    collection = np.empty([len(VAR_LIST), 2, 4])

    for thisiter in range(0, niter):
        status = 0
        while status == 0:
            myjob = comm.recv(source=0, tag=1)
            status = comm.recv(source=0, tag=2)

            if status == 0:
                ierr = postproc(myjob, collection)
                comm.send(rank, dest=0, tag=3)
                comm.send(myjob, dest=0, tag=4)
                comm.send(collection, dest=0, tag=5)

    MPI.Finalize()
