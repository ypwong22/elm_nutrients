#!/bin/csh -f
#SBATCH --time=06:00:00
#SBATCH -J process_ensemble_sensitivity
#SBATCH --nodes=2
#SBATCH -p batch
#SBATCH --mem=64G
#SBATCH --ntasks-per-node 32
#SBATCH --export=ALL
#SBATCH -A ccsi

# Before running this script, in the external environment
# module load PE-gnu
# module load mpich
# conda activate olmt

echo `which python`
cd ~/Git/phenology_elm
/software/dev_tools/swtree/cs400_centos7.2_pe2016-08/openmpi/3.0.0/centos7.2_gnu5.3.0/bin/mpirun -np 32 python process_ensemble_sensitivity.py
