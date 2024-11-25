#!/bin/csh
#SBATCH --time=6:0:00
#SBATCH -J ens_UQ_20240323_OAT
#SBATCH --nodes=1
#SBATCH -A CLI185
#SBATCH -p batch
#SBATCH --ntasks-per-node 100

#setenv LD_LIBRARY_PATH /sw/baseline/spack-envs/base/opt/linux-rhel8-zen3/gcc-12.2.0/netlib-lapack-3.11.0-lpwyqsehj7wuz2i45umfhwa5ymv2dz5b/lib64:/sw/baseline/spack-envs/base/opt/linux-rhel8-zen3/gcc-12.2.0/openmpi-4.0.4-bxes2wvty3q7v55qep7hiuud6rocd4bl/lib:/sw/baseline/gcc/12.2.0/lib64:/ccsopen/home/zdr/opt/lib:${HOME}/.conda/envs/olmt/lib
cd ${HOME}/Git/phenology_elm
srun -n 100 ${HOME}/.conda/envs/olmt/bin/python -u process_ensemble_sensitivity.py
