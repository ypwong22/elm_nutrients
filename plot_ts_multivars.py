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


"""
parser = OptionParser()
parser.add_option('--plot', dest = 'level', default = '', help = 'Write the chamber number.')
parser.add_option('--date', dest = 'date', default = '', help = 'Simulation prefix.')
parser.add_option('--extra', dest = 'extra', default = '', help = 'Extra prefix.')
(options, args) = parser.parse_args()

# obsolete
# 20211008: original
# 20211201: modified
# 20220103: modified, air temperature-based pheno
# 20220108: modified, air temperature-based pheno, dynamic root
# 20220407: modified, root phenology; example call
#           python plot_ts_multipanel.py --date 20220407 --extra _stopdate171_switch1 --plot 04
"""

date  = '20221212'


###################################################################################################
# Plot the time series of vegetation variables
###################################################################################################
var_list = ['TLAI', 'GPP', 'QVEGT', 'FROOTC']
unit_list = ['', 'gC m-2 day-1', 'mm day-1', 'gC m-2']
# var_list = ['FROOTC','FROOTC_XFER', 'LEAFC_TO_LITTER', 'FROOTC_TO_LITTER']
# unit_list = ['gC m-2', 'gC m-2', 'gC m-2 day-1', 'gC m-2 day-1']
# var_list = ['LEAFC_STORAGE_TO_XFER','LEAFC_XFER','LEAFC_XFER_TO_LEAFC','FROOTC_STORAGE','FROOTC_XFER'] # 'LITFALL'
# unit_list = ['gC m-2 day-1', 'gC m-2', 'gC m-2 day-1', 'gC m-2', 'gC m-2']
path_out = os.path.join(os.environ['PROJDIR'], 'Phenology_ELM', 'output_elm')
plot_list = [4, 10, 11, 16, 19, 6, 8, 13, 17, 20]


# -------------------------------------------------------------------------------------------------
# PFT level
# -------------------------------------------------------------------------------------------------
pft_list = {'needleleaf evergreen boreal tree': 2, 
            'needleleaf deciduous boreal tree': 3, 
            'broadleaf deciduous boreal shrub': 11}
clist = ['#fcbba1', '#fc9272', '#fb6a4a', '#de2d26', '#a50f15', '#c6dbef', '#9ecae1', '#6baed6', '#3182bd', '#08519c']


fig, axes = plt.subplots(len(var_list), 3, figsize = (12, 14), sharex = True)
for ii, var in enumerate(var_list):
    for jj, pft in enumerate(list(pft_list.keys())):

        h = [None] * 10
        for kk, plot in enumerate(plot_list):
            path_data = os.path.join(os.environ['PROJDIR'], 'E3SM', 'output', f'{date}_plot{plot:02d}_US-SPR_ICB20TRCNPRDCTCBC', 'run')
            flist = sorted(glob(os.path.join(path_data, '*.h2.*.nc')))

            # osbolete # autocoversion misidentified PROOT_ONSET_OFFSET as a datetime variable because its unit "days"
            hr = xr.open_mfdataset(flist, decode_times = False)
            units, reference_date = hr.time.attrs['units'].split('since')
            tvec = pd.date_range(start = reference_date, end = '2021-01-01', freq='D')
            tvec = tvec[(tvec.month != 2) | (tvec.day != 29)]
            hr['time'] = tvec

            ax = axes[ii, jj]

            tvec = [datetime.strptime(x, '%Y-%m-%d %H') for x in hr['time'].indexes['time'].strftime('%Y-%m-%d %H')]
            retype = pd.DataFrame(hr[var][:, pft_list[pft]].values * 0.64 + hr[var][:, pft_list[pft] + 17].values * 0.64, index = tvec)

            hr.close()

            if var in ['GPP', 'QVEGE', 'QVEGT', 'LEAFC_TO_LITTER', 'FROOTC_TO_LITTER', 'LEAFC_STORAGE_TO_XFER', 'LEAFC_XFER_TO_LEAFC']:
                # s to day
                retype = retype * 24 * 3600

            if var == 'GPP':
                retype2 = pd.DataFrame(savgol_filter(retype.values, 81, 3, mode = 'constant', cval = np.nan), index = tvec)

            if var == 'QVEGE':
                retype = retype.iloc[1:]

            if var_list[ii] == 'GPP':
                htemp, = ax.plot(retype.index, retype.values, '-', color = clist[kk], lw = 0.5)
                ax.plot(retype2.index, retype2.values, '-', color = clist[kk], lw = 1)
            else:
                htemp, = ax.plot(retype.index, retype.values, '-', color = clist[kk])
            h[kk] = htemp
            
            if jj == 0:
                ax.set_ylabel(var_list[ii] + ' ' + unit_list[ii])
            if ii == 0:
                ax.set_title(pft)

            #ax.set_xlim([0.5, 365.5])
            #ax.set_xticks(range(1,13))
ax.legend(h, [f'plot{plot:02d}' for plot in plot_list], loc = 'upper right', bbox_to_anchor = (0.5, -.2), ncol = 3)
fig.savefig(os.path.join(path_out, date, 'ts_bypft.png'), dpi = 600., bbox_inches = 'tight')
plt.close()