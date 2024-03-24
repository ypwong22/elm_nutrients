#!/bin/csh
#SBATCH --time=24:0:00
#SBATCH -J ens_UQ_20240311_US-SPR_ICB20TRCNPRDCTCBC
#SBATCH --nodes=2
#SBATCH -A CLI185
#SBATCH -p batch
#SBATCH --ntasks-per-node 128

setenv LD_LIBRARY_PATH /sw/baseline/spack-envs/base/opt/linux-rhel8-zen3/gcc-12.2.0/netlib-lapack-3.11.0-lpwyqsehj7wuz2i45umfhwa5ymv2dz5b/lib64:/sw/baseline/spack-envs/base/opt/linux-rhel8-zen3/gcc-12.2.0/openmpi-4.0.4-bxes2wvty3q7v55qep7hiuud6rocd4bl/lib:/sw/baseline/gcc/12.2.0/lib64:/ccsopen/home/zdr/opt/lib:${HOME}/.conda/envs/olmt/lib
cd ${HOME}/Git/phenology_elm
srun -n 128 ${HOME}/.conda/envs/olmt/bin/python process_ensemble_sensitivity.py
