""" Plot multiple variables on different panels on the same figure. 
    Each simulation is a separate plot. """
import itertools as it
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import numpy as np
import os
from glob import glob
from datetime import datetime, timedelta
from scipy.signal import savgol_filter
from optparse import OptionParser


parser = OptionParser()
parser.add_option('--plot', dest = 'level', default = '', help = 'Write the chamber number.')
parser.add_option('--date', dest = 'date', default = '', help = 'Simulation prefix.')
parser.add_option('--extra', dest = 'extra', default = '', help = 'Extra prefix.')
(options, args) = parser.parse_args()


#level = 'plot10' # 'ambient'
level = 'plot' + options.level

# 20211008: original
# 20211201: modified
# 20220103: modified, air temperature-based pheno
# 20220108: modified, air temperature-based pheno, dynamic root
# 20220407: modified, root phenology; example call
#           python plot_ts_multipanel.py --date 20220407 --extra _stopdate171_switch1 --plot 04
date  = options.date
extra = options.extra
#date = '20220407_stopdate171_switch2'
#extra = ''

path_data = os.path.join(os.environ['SCRATCHDIR'], 'E3SM', 'output', 
                         date + extra + '_' + level + '_US-SPR_ICB20TRCNPRDCTCBC','run')
print(path_data)
path_out = os.path.join(os.environ['PROJDIR'], 'Phenology_ELM', 'output_elm')
###################################################################################################
# Plot the time series of vegetation variables
###################################################################################################
var_list = ['TLAI', 'GPP', 'QVEGT', 'FROOTC']
#var_list = ['FROOTC','PROOT_ONSET_FLAG','PROOT_ONSET_OFFSET','PROOT_GDD2','PROOT_CWD2','PROOT_SUMFRAC']
#var_list = ['LEAFC_STORAGE_TO_XFER','LEAFC_XFER','LEAFC_XFER_TO_LEAFC','FROOTC_STORAGE','FROOTC_XFER'] # 'LITFALL'

# -------------------------------------------------------------------------------------------------
# PFT level
# -------------------------------------------------------------------------------------------------
pft_list = {'needleleaf evergreen boreal tree': 2, 
            'needleleaf deciduous boreal tree': 3, 
            'broadleaf deciduous boreal shrub': 11}
clist = ['#1b9e77', '#d95f02', '#7570b3']

flist = sorted(glob(os.path.join(path_data, '*.h2.*.nc')))

fig, axes = plt.subplots(len(var_list), 1, figsize = (12, 14), sharex = True)
for ii in range(len(var_list)):
    # autocoversion misidentified PROOT_ONSET_OFFSET as a datetime variable because its unit "days"
    hr = xr.open_mfdataset(flist, decode_times = False)
    units, reference_date = hr.time.attrs['units'].split('since')
    tvec = pd.date_range(start = reference_date, end = '2021-01-01', freq='D')
    tvec = tvec[(tvec.month != 2) | (tvec.day != 29)]
    hr['time'] = tvec

    ax = axes.flat[ii]

    count = 0
    nlist = []
    for jj, kk in pft_list.items():
        tvec = [datetime.strptime(x, '%Y-%m-%d %H') for x in \
                hr['time'].indexes['time'].strftime('%Y-%m-%d %H')]
        retype = pd.DataFrame(hr[var_list[ii]][:, kk].values, index = tvec)

        if var_list[ii] in ['GPP', 'QVEGE', 'QVEGT', 'LEAFC_STORAGE', 'FROOTC']:
            # s to day, add up to month, divide by years
            retype = retype * 24 * 3600

        if var_list[ii] == 'GPP':
            retype2 = pd.DataFrame(savgol_filter(hr[var_list[ii]][:, kk].values, 81, 3),
                                   index = tvec)
            retype2 = retype2 * 24 * 3600

        if (ii == 0) & (count == 0):
            collect = pd.DataFrame(np.nan, index = retype.index,
                columns = pd.MultiIndex.from_product([var_list, sorted(pft_list.keys())]))

        collect.loc[:, (var_list[ii],jj)] = retype.values

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
fig.savefig(os.path.join(path_out, date + extra, level + '_ts_bypft.png'), dpi = 600., bbox_inches = 'tight')
plt.close()

collect.to_csv(os.path.join(path_out, date + extra, level + '_ts_bypft.csv'))
