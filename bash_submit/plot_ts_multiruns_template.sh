#PBS -S /bin/bash
#PBS -A ACF-UTK0011
#PBS -l nodes=1:ppn=16,walltime=10:00:00

cd ~/Git/phenology_elm

conda deactivate
conda activate myenv
python plot_ts_multiruns.py --level REPLACE
