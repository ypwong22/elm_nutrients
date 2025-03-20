""" Check the tree transpiration for Jeff using the default model """
import os
from glob import glob
import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from utils.constants import chamber_list_complete, chamber_list_names_complete


path_in = os.path.join(os.environ['E3SM_ROOT'], 'output', 
                       '20231116_US-SPR_ICB20TRCNPRDCTCBC', 'spruce_treatments')
path_out = os.path.join(os.environ['PROJDIR'], 'ELM_Phenology', 'output', 'extract', '20231116')

pftfracs = [0.36, 0.14, 0.25, 0.25]

transpiration = {}
for plot in chamber_list_complete:
    flist = sorted(glob(os.path.join(path_in, f'plot{plot:02d}_US-SPR_ICB20TRCNPRDCTCBC', 'run', 
                                     '*.h2.*.nc')))[:-1]
    hr = xr.open_mfdataset(flist)
    for pft, name in zip([2,3,11,12], ['Spruce','Tamarack','Shrub','Moss']):
        temp = hr['QVEGT'][:, pft].resample(time = '1YS').mean()
        tvec = [pd.to_datetime(f'{t.year}-{t.month:02d}-01') for t in temp['time'].values]
        transpiration[(plot,name)] = pd.Series(temp.values * 24 * 3600, index = tvec)
    hr.close()

transpiration = pd.DataFrame(transpiration)


fig, axes = plt.subplots(4, 3, figsize = (12, 12), sharex = True, sharey = True)
fig.subplots_adjust(wspace = 0.05)
for i,plot in enumerate(chamber_list_complete):
    ax = axes.flat[i]

    cum = 0

    for j, pft in enumerate(['Spruce','Tamarack','Shrub','Moss']):
        ax.plot(transpiration.index, transpiration[(plot,pft)], '-o', label = pft)

        cum = cum + transpiration[(plot,pft)] * pftfracs[j]
    
    ax.plot(transpiration.index, cum, '-ok', lw = 2, label = 'Average')
    ax.set_title(chamber_list_names_complete[i])

    if np.mod(i, 3) == 0:
        ax.set_ylabel('Transpiration (mm/day)')

ax.legend(loc = (1.2, 0.5))
axes.flat[-1].axis('off')
fig.savefig(os.path.join(path_out, 'QVEGT.png'), dpi = 600., bbox_inches = 'tight')
plt.close(fig)