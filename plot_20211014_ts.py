import itertools as it
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import numpy as np
import os
from glob import glob
from datetime import datetime, timedelta
from scipy.signal import savgol_filter

level = 'plot10' # 'ambient'
date  = '20211201' # '20211201' (modified), '20211008' (original)
path_data = os.path.join(os.environ['SCRATCHDIR'], 'E3SM', 'output', 
                         date + '_' + level + '_US-SPR_ICB20TRCNPRDCTCBC','run')
path_out = os.path.join(os.environ['PROJDIR'], 'Phenology_ELM', 'output_elm')


###################################################################################################
# Plot the time series of vegetation variables
###################################################################################################
var_list = ['TLAI', 'GPP', 'QVEGE', 'QVEGT', 'LEAFC_STORAGE']

# -------------------------------------------------------------------------------------------------
# PFT level
# -------------------------------------------------------------------------------------------------
pft_list = {'needleleaf evergreen boreal tree': 2, 
            'needleleaf deciduous boreal tree': 3, 
            'broadleaf deciduous boreal shrub': 11}
clist = ['#1b9e77', '#d95f02', '#7570b3']

flist = glob(os.path.join(path_data, '*h2*.nc'))

fig, axes = plt.subplots(len(var_list), 1, figsize = (12, 12), sharex = True)
for ii in range(len(var_list)):
    hr = xr.open_mfdataset(flist, decode_times = True)

    ax = axes.flat[ii]

    count = 0
    nlist = []
    for jj, kk in pft_list.items():
        tvec = [datetime.strptime(x, '%Y-%m-%d %H') for x in \
                hr['time'].indexes['time'].strftime('%Y-%m-%d %H')]
        retype = pd.DataFrame(hr[var_list[ii]][:, kk].values, index = tvec)
        if var_list[ii] != 'TLAI':
            # s to day, add up to month, divide by years
            retype = retype * 24 * 3600

        if var_list[ii] == 'GPP':
            retype2 = pd.DataFrame(savgol_filter(hr[var_list[ii]][:, kk].values, 81, 3),
                                   index = tvec)
            retype2 = retype2 * 24 * 3600

        if kk == 2:
            if ii == 0:
                collect = pd.DataFrame(np.nan, index = retype.index, columns = var_list)
            collect.loc[:, var_list[ii]] = retype.values

        if var_list[ii] == 'QVEGE':
            retype = retype.iloc[1:]

        if var_list[ii] == 'GPP':
            ax.plot(retype.index, retype.values, '-', color = clist[count],
                    lw = 0.5)
            ax.plot(retype2.index, retype2.values, '-', color = clist[count], lw = 1)
        else:
            ax.plot(retype.index, retype.values, '-', color = clist[count])
        ax.set_title(var_list[ii])
        #ax.set_xlim([0.5, 365.5])
        #ax.set_xticks(range(1,13))

        count += 1
        nlist.append(jj)
    hr.close()
ax.legend(nlist, loc = 'upper right', bbox_to_anchor = (0.5, -.2))
fig.savefig(os.path.join(path_out, date + '_' + level + '_1_ts.png'), dpi = 600., bbox_inches = 'tight')
plt.close()

collect.to_csv(os.path.join(path_out, date + '_' + level + '_1_ts.csv'))
