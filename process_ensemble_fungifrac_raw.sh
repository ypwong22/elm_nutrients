#!/bin/csh
#SBATCH --time=24:00:00
#SBATCH -J fungifrac
#SBATCH --nodes=1
#SBATCH -A CLI185
#SBATCH -p batch_ccsi
#SBATCH --ntasks-per-node 3

#setenv LD_LIBRARY_PATH /sw/baseline/spack-envs/base/opt/linux-rhel8-zen3/gcc-12.2.0/netlib-lapack-3.11.0-lpwyqsehj7wuz2i45umfhwa5ymv2dz5b/lib64:/sw/baseline/spack-envs/base/opt/linux-rhel8-zen3/gcc-12.2.0/openmpi-4.0.4-bxes2wvty3q7v55qep7hiuud6rocd4bl/lib:/sw/baseline/gcc/12.2.0/lib64:/ccsopen/home/zdr/opt/lib:${HOME}/.conda/envs/olmt/lib

# From load-balancing perspective, 128 cores seems inefficient. Fewer cores are better
cd ${HOME}/Git/phenology_elm
srun -n 3 ${HOME}/.conda/envs/olmt/bin/python -u process_ensemble_fungifrac_raw.py # u flag to print un-buffered output
