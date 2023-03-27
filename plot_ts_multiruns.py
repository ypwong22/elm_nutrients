""" Compare the time series of multiple runs on the same graph. Plot level. """
import itertools as it
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import numpy as np
import os
from glob import glob
from datetime import datetime, timedelta
from scipy.signal import savgol_filter
from scipy.stats import pearsonr
from optparse import OptionParser
from utils.constants import *
from utils.paths import *


var_list    = ['LITFALL', 'GROSS_NMIN', 'GROSS_PMIN']
unit_list   = ['gC m-2 day-1', 'gN m-2 day-1', 'gP m-2 day-1']
sims_prefix = ['20221212', '20221231', '20230212']
sims_names  = ['Default', 'Evergreen', 'EvgrRoot']
clist = ['#fcbba1', '#c6dbef', '#fc9272', '#9ecae1', '#fb6a4a', '#6baed6', '#de2d26', '#3182bd', '#a50f15', '#08519c']
pft_list = [2, 3, 11]


###########################################################################################
# Time series
###########################################################################################
fig, axes = plt.subplots(len(var_list), len(sims_prefix), figsize = (12, 3 * len(var_list)), sharex = True, sharey = False)
fig.subplots_adjust(hspace = 0.05, wspace = 0.05)
for i, varname in enumerate(var_list):
    ymin =  999999
    ymax = -999999

    for j, prefix in enumerate(sims_prefix):
        ax = axes[i, j]

        h = [None] * len(chamber_list)
        for k, plot in enumerate(chamber_list):
            path_data = os.path.join(os.environ['PROJDIR'], 'E3SM', 'output', f'{prefix}_plot{plot:02d}_US-SPR_ICB20TRCNPRDCTCBC', 'run')
            flist = sorted(glob(os.path.join(path_data, '*.h0.*.nc')))

            hr = xr.open_mfdataset(flist, decode_times = False)
            units, reference_date = hr.time.attrs['units'].split('since')
            tvec = pd.date_range(start = reference_date, end = '2020-12-31', freq='MS')
            hr['time'] = tvec
            retype = pd.DataFrame(hr[varname][:, 0].values * 0.64 + hr[varname][:,  1].values * 0.36, index = tvec)
            if ('/s' in hr[varname].attrs['units']) and ('day-1' in unit_list[i]):
                # s to year
                retype = retype * 24 * 3600 * 365
            hr.close()

            h[k], = ax.plot(retype.index, retype.values, color = clist[k])

        if j == 0:
            ax.set_ylabel(f'{varname} {unit_list[i]}')
        else:
            ax.set_yticklabels([])
        if i == 0:
            ax.set_title(sims_names[j])

        ymin_, ymax_ = ax.get_ylim()
        ymin = min(ymin, ymin_)
        ymax = max(ymax, ymax_)

    for j in range(len(sims_prefix)):
        axes[i, j].set_ylim([ymin, ymax])

ax.legend(h, chamber_list_names, loc = [-2.2, -0.4], ncol = 5)
fig.savefig(os.path.join(path_out_elm, f'ts_multiruns_soil.png'), dpi = 600.,  bbox_inches = 'tight')
plt.close(fig)


###########################################################################################
# Mean annual cycle
###########################################################################################
fig, axes = plt.subplots(len(var_list), len(sims_prefix), figsize = (3 * len(var_list), 12), sharex = True, sharey = False)
fig.subplots_adjust(hspace = 0.05, wspace = 0.05)
for i, varname in enumerate(var_list):
    ymin =  999999
    ymax = -999999

    for j, prefix in enumerate(sims_prefix):
        ax = axes[i, j]

        h = [None] * len(chamber_list)
        for k, plot in enumerate(chamber_list):
            path_data = os.path.join(os.environ['PROJDIR'], 'E3SM', 'output', f'{prefix}_plot{plot:02d}_US-SPR_ICB20TRCNPRDCTCBC', 'run')
            flist = sorted(glob(os.path.join(path_data, '*.h0.*.nc')))

            # osbolete # autocoversion misidentified PROOT_ONSET_OFFSET as a datetime variable because its unit "days"
            hr = xr.open_mfdataset(flist, decode_times = False)
            units, reference_date = hr.time.attrs['units'].split('since')
            tvec = pd.date_range(start = reference_date, end = '2020-12-31', freq='MS')
            hr['time'] = tvec
            retype = pd.DataFrame(hr[varname][:, 0].values * 0.64 + hr[varname][:,  1].values * 0.36, index = tvec)
            hr.close()

            if ('/s' in hr[varname].attrs['units']) and ('day-1' in unit_list[i]):
                # s to year
                retype = retype * 24 * 3600 * 365

            retype = retype.groupby(retype.index.month).mean()

            h[k], = ax.plot(retype.index, retype.values, color = clist[k])

        if j == 0:
            ax.set_ylabel(f'{varname} {unit_list[i]}')
        else:
            ax.set_yticklabels([])
        if i == 0:
            ax.set_title(sims_names[j])

        if i == (len(var_list) - 1):
            ax.set_xlabel('Month')

        ymin_, ymax_ = ax.get_ylim()
        ymin = min(ymin, ymin_)
        ymax = max(ymax, ymax_)

    for j in range(len(sims_prefix)):
        axes[i, j].set_ylim([ymin, ymax])

ax.legend(h, chamber_list_names, loc = [-2.2, -0.4], ncol = 5)
fig.savefig(os.path.join(path_out_elm, f'ts_multiruns_seasonal_soil.png'), dpi = 600.,  bbox_inches = 'tight')
plt.close(fig)
