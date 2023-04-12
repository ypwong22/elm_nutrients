import os
import numpy as np

f = open(os.environ['HOME'] + '/models/OLMT/UQ_output/UQ_default_US-SPR_ICB20TRCNPRDCTCBC/MCMC_output/MCMC_chain.txt', 'r')
allines = f.read().split('\n')
f.close()

np.random.seed(125)
array = np.random.choice(range(len(allines)), 4000, replace = False)
sublines = [allines[m] for m in array]
f = open(os.environ['HOME'] + '/Git/phenology_elm/temp/MCMC_subset.txt', 'w')
for line in sublines:
    f.write(line + '\n')
f.close()