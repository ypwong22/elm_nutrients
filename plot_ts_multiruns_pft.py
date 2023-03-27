""" Compare the time series of multiple runs on the same graph. """
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


var_list    = ['GPP', 'LEAFC', 'LEAFC_TO_LITTER', 'FROOTC', 'FROOTC_TO_LITTER']
unit_list   = ['gC m-2 year-1', 'gC m-2', 'gC m-2 year-1', 'gC m-2', 'gC m-2 year-1']
sims_prefix = ['20221212', '20221231', '20230101']
sims_names  = ['Default', 'Evergreen', 'EvgrRoot']
clist = ['#fcbba1', '#c6dbef', '#fc9272', '#9ecae1', '#fb6a4a', '#6baed6', '#de2d26', '#3182bd', '#a50f15', '#08519c']
pft_list = [2, 3, 11]


###########################################################################################
# Time series
###########################################################################################
for p, pft in enumerate(pft_list):

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
                flist = sorted(glob(os.path.join(path_data, '*.h2.*.nc')))

                # osbolete # autocoversion misidentified PROOT_ONSET_OFFSET as a datetime variable because its unit "days"
                hr = xr.open_mfdataset(flist, decode_times = False)
                units, reference_date = hr.time.attrs['units'].split('since')
                tvec = pd.date_range(start = reference_date, end = '2021-01-01', freq='D')
                tvec = tvec[(tvec.month != 2) | (tvec.day != 29)]
                hr['time'] = tvec
                tvec = [datetime.strptime(x, '%Y-%m-%d %H') for x in hr['time'].indexes['time'].strftime('%Y-%m-%d %H')]
                retype = pd.DataFrame(hr[varname][:, pft].values * 0.64 + hr[varname][:,  pft + 17].values * 0.36, index = tvec)
                hr.close()

                if varname in ['GPP', 'QVEGE', 'QVEGT', 'LEAFC_TO_LITTER', 'FROOTC_TO_LITTER', 'LEAFC_STORAGE_TO_XFER', 'LEAFC_XFER_TO_LEAFC']:
                    # s to year
                    retype = retype * 24 * 3600 * 365
                if varname == 'GPP':
                    retype2 = pd.DataFrame(savgol_filter(retype.values, 81, 3, mode = 'constant', cval = np.nan), index = tvec)
                if varname == 'QVEGE':
                    retype = retype.iloc[1:]

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
    fig.savefig(os.path.join(path_out_elm, f'ts_multiruns_{pft}.png'), dpi = 600.,  bbox_inches = 'tight')
    plt.close(fig)


###########################################################################################
# Mean annual cycle
###########################################################################################
for p, pft in enumerate(pft_list):

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
                flist = sorted(glob(os.path.join(path_data, '*.h2.*.nc')))

                # osbolete # autocoversion misidentified PROOT_ONSET_OFFSET as a datetime variable because its unit "days"
                hr = xr.open_mfdataset(flist, decode_times = False)
                units, reference_date = hr.time.attrs['units'].split('since')
                tvec = pd.date_range(start = reference_date, end = '2021-01-01', freq='D')
                tvec = tvec[(tvec.month != 2) | (tvec.day != 29)]
                hr['time'] = tvec
                tvec = [datetime.strptime(x, '%Y-%m-%d %H') for x in hr['time'].indexes['time'].strftime('%Y-%m-%d %H')]
                retype = pd.DataFrame(hr[varname][:, pft].values * 0.64 + hr[varname][:,  pft + 17].values * 0.36, index = tvec)
                hr.close()

                if varname in ['GPP', 'QVEGE', 'QVEGT', 'LEAFC_TO_LITTER', 'FROOTC_TO_LITTER', 'LEAFC_STORAGE_TO_XFER', 'LEAFC_XFER_TO_LEAFC']:
                    # s to year
                    retype = retype * 24 * 3600 * 365
                if varname == 'GPP':
                    retype2 = pd.DataFrame(savgol_filter(retype.values, 81, 3, mode = 'constant', cval = np.nan), index = tvec)
                if varname == 'QVEGE':
                    retype = retype.iloc[1:]

                dayofyear = []
                for yy in np.unique(retype.index.year[:-1]):
                    dayofyear = dayofyear + list(np.arange(1, 366))
                dayofyear = dayofyear + [1] # restart file
                dayofyear = np.array(dayofyear)
                retype = retype.groupby(dayofyear).mean()

                h[k], = ax.plot(retype.index, retype.values, color = clist[k])

            if j == 0:
                ax.set_ylabel(f'{varname} {unit_list[i]}')
            else:
                ax.set_yticklabels([])
            if i == 0:
                ax.set_title(sims_names[j])

            if i == (len(var_list) - 1):
                ax.set_xlabel('Day of year')

            ymin_, ymax_ = ax.get_ylim()
            ymin = min(ymin, ymin_)
            ymax = max(ymax, ymax_)

        for j in range(len(sims_prefix)):
            axes[i, j].set_ylim([ymin, ymax])

    ax.legend(h, chamber_list_names, loc = [-2.2, -0.4], ncol = 5)
    fig.savefig(os.path.join(path_out_elm, f'ts_multiruns_seasonal_{pft}.png'), dpi = 600.,  bbox_inches = 'tight')
    plt.close(fig)


"""
# Sensitivity of GPP to the variables
for prefix, model in zip(sims_prefix, sims_names):
    for loc, suffix in it.product(['hummock', 'hollow'], ['_pima', '_lala', '_shrub']):
        fig, axes = plt.subplots(len(chamber_list), 4, figsize = (3 * len(chamber_list), 16))
        for j,ch in enumerate(chamber_list):
            data_list = pd.read_csv(os.path.join(path_out_elm, prefix, f'plot{ch}_ts_extract.csv'),
                                    index_col = 0, parse_dates = True, header = [0,1])
            for i, (varname,unit) in enumerate([('FSDS','W m-2'), 
                                                ('TLAI',''),
                                                ('FROOTC','gC m-2')]):
                ax = axes[j,i]
                if varname == 'FSDS':
                    x = data_list[(loc, varname)].values
                else:
                    x = data_list[(loc, varname + suffix)].values
                y = data_list[(loc, 'GPP' + suffix)].values
                ax.plot(x, y, 'o', markersize = 3)
                r, pval = pearsonr(x, y)
                ax.text(0.05, 0.85, '{:.4f}, p = {:.3f}'.format(r, pval), 
                        transform = ax.transAxes)

                if j == 0:
                    ax.set_title(varname)
                elif j == (len(chamber_list)-1):
                    ax.set_xlabel(unit)
                if i == 0:
                    ax.set_ylabel('GPP (gC m-2)')
                if i == 2:
                    ax.text(-0.5, 0.5, 
                            f'P{ch} ({chamber_levels[ch][0]}$^o$C, {chamber_levels[ch][1]})', rotation = 90, transform = ax.transAxes)

        for j,ch in enumerate(chamber_list):
            data_list = pd.read_csv(os.path.join(path_out_elm, prefix,
                                                 f'plot{ch}_ts_extract.csv'),
                                    index_col = 0, parse_dates = True, header = [0,1])

            ax = axes[j,3]

            x = data_list[(loc, 'TLAI' + suffix)].values
            # y = data_list[(loc, 'FSDS')].values
            y = data_list[(loc, 'FROOTC' + suffix)].values
            z = data_list[(loc, 'GPP' + suffix)].values
            cf = ax.scatter(x, y, c = z, s = 3, cmap = 'Spectral')
            plt.colorbar(cf, ax = ax)
            ax.set_title('GPP (gC m-2)')
            ax.set_xlabel('TLAI')
            #ax.set_ylabel('FSDS (W m-2)')
            ax.set_ylabel('FROOTC (gC m-2)')

        fig.savefig(os.path.join(path_out_elm, 'ts_multiruns', model + '_' + loc + \
                                 '_photosyntheis' + suffix + '.png'), dpi = 600., 
                    bbox_inches = 'tight')
        plt.close(fig)
"""