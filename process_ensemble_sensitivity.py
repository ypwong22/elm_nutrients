#!/usr/bin/env python
# /software/dev_tools/swtree/cs400_centos7.2_pe2016-08/openmpi/3.0.0/centos7.2_gnu5.3.0/bin/mpirun -np 64 python process_ensemble_sensitivity.py
import sys, os, time
import numpy as np
from netCDF4 import Dataset
from scipy.stats import linregress, t
from mpi4py import MPI
from utils.paths import *
from utils.constants import chamber_list_names_complete
import pandas as pd


# number of simulations
# N = 1000
# PREFIX = "UQ_default2_optimized"

# jobid = 4113483, burst
# N = 4000
# PREFIX = "UQ_default2"

# jobid = 4113489
N = 2000
PREFIX = "UQ_root"

# number of ensembles to save in each bin file
# this avoids having difficulty in dumping file
BLOCK = 200
if np.mod(N, BLOCK) != 0:
    raise Exception("N must be a multiply of BLOCK")


RUNROOT = os.path.join(os.environ["PROJDIR"], "E3SM", "output")
VAR_LIST = ["ANPPtree", "ANPPshrub", "NPPmoss", "BGNPP", "NPP", "HR", "NEE"]
YEAR_LIST = range(2016, 2021)  # skip 2015 and 2021 because no observation/no model data


MOSSFRAC = pd.read_excel(
    os.path.join(path_input, "Sphagnum_fraction.xlsx"),
    index_col=0,
    skiprows=1,
    engine="openpyxl",
).drop(["plot", "Temp", "CO2"], axis=1)
MOSSFRAC[2015] = MOSSFRAC[2016]
# MOSSFRAC = MOSSFRAC.drop(2021, axis=1)


# Define function to perform ensemble member post-processing
def postproc(thisjob, collection):
    ierr = 0

    casename = f"{PREFIX}_US-SPR_ICB20TRCNPRDCTCBC"
    baserundir = os.path.join(RUNROOT, "UQ", casename, f"g{thisjob:05g}")
    print(thisjob)

    values = np.empty([len(chamber_list_names_complete), len(YEAR_LIST), len(VAR_LIST)])
    t2m = np.empty([len(chamber_list_names_complete), len(YEAR_LIST)])

    for i, folder in enumerate(chamber_list_names_complete):
        # DEBUG
        ##casename = f'20230526_plot{chamber_list_complete[i]:02g}_US-SPR_ICB20TRCNPRDCTCBC'
        ##baserundir = os.path.join(RUNROOT, casename, "run")

        for j, year in enumerate(YEAR_LIST):
            fname_pft = os.path.join(
                baserundir,
                folder,
                f"{casename}.clm2.h1.{year}-01-01-00000.nc",
            )
            fname_col = os.path.join(
                baserundir,
                folder,
                f"{casename}.clm2.h0.{year}-01-01-00000.nc",
            )

            # DEBUG
            ##fname_pft = os.path.join(
            ##    baserundir,
            ##    f"{casename}.clm2.h2.{year}-01-01-00000.nc",
            ##)
            ##fname_col = os.path.join(
            ##    baserundir,
            ##    f"{casename}.clm2.h1.{year}-01-01-00000.nc",
            ##)
            ##print(fname_pft, fname_col)

            # note some parameter combinations can fail
            if not (os.path.exists(fname_pft) or os.path.exists(fname_col)):
                ierr = -1
                return ierr

            hr = Dataset(fname_pft)
            hr2 = Dataset(fname_col)

            for k, varname in enumerate(VAR_LIST):
                if varname in ["NPPmoss", "ANPPtree", "ANPPshrub", "BGNPP"]:
                    if varname == "NPPmoss":
                        var = (
                            hr["NPP"][:, :].mean(axis=0)
                            * MOSSFRAC.loc[folder, year]
                            / 100
                        )
                        values[i, j, k] = var[12] * 0.64 + var[12 + 17] * 0.36
                    elif varname == "ANPPtree":
                        var = hr["AGNPP"][:, :].mean(axis=0)
                        values[i, j, k] = (
                            var[2] * 0.64 + var[2 + 17] * 0.36
                        ) * 0.36 + (var[3] * 0.64 + var[3 + 17] * 0.36) * 0.14
                    elif varname == "ANPPshrub":
                        var = hr["AGNPP"][:, :].mean(axis=0)
                        values[i, j, k] = (var[11] * 0.64 + var[11 + 17] * 0.36) * 0.25
                    elif varname == "BGNPP":
                        # tree and shrub only
                        var = hr["BGNPP"][:, :].mean(axis=0)
                        values[i, j, k] = (
                            (var[2] * 0.64 + var[2 + 17] * 0.36) * 0.36
                            + (var[3] * 0.64 + var[3 + 17] * 0.36) * 0.14
                            + (var[11] * 0.64 + var[11 + 17] * 0.36) * 0.25
                        )
                else:
                    var = hr2[varname][:, :].mean(axis=0)
                    values[i, j, k] = var[0] * 0.64 + var[1] * 0.36

            var = hr2["TBOT"][:, :].mean(axis=0) - 273.15
            t2m[i, j] = var[0] * 0.64 + var[1] * 0.36

            hr.close()
            hr2.close()

    # ('amb', 'elev') x (mean, mean_std, slope, slope_std)
    # collection = np.empty([len(VAR_LIST), 2, 4])
    # ambient -> elevated CO2
    for i, temp in enumerate(
        [
            np.arange(1, len(chamber_list_names_complete), 2),
            np.arange(2, len(chamber_list_names_complete), 2),
        ]
    ):
        # Note this is the ambient or T0.00CO2 chamber
        for k, _ in enumerate(VAR_LIST):
            collection[k, i, 0] = np.mean(values[temp[0], :, k])
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

workdir = os.getcwd()

niter = int(N / BLOCK)

for b in range(niter):
    print("rank = ", rank, "b = ", b)

    if rank == 0:
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
        collection = np.full([len(VAR_LIST), 2, 4], np.nan)

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

MPI.Finalize()
