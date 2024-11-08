#!/bin/csh
#SBATCH --time=01:0:00
#SBATCH -J extract_ts_genvars
#SBATCH --nodes=1
#SBATCH -A CLI185
#SBATCH -p batch_ccsi
#SBATCH --ntasks-per-node 1

cd ${HOME}/Git/phenology_elm
${HOME}/.conda/envs/olmt/bin/python extract_ts_genvars.py
