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


date   = '20221212' # '20230214', '20221231', '20221212'
varset = 'vegc' # ['vegc', 'root', 'xfer', 'soil', 'diag']


###################################################################################################
# Plot the time series of vegetation variables
###################################################################################################
if varset == 'vegc':
    var_list = ['TLAI', 'GPP', 'QVEGT', 'FROOTC']
    unit_list = ['', 'gC m-2 day-1', 'mm day-1', 'gC m-2']
elif varset == 'root':
    var_list = ['FROOTC', 'FROOTC_XFER', 'FROOTC_STORAGE_TO_XFER', 'FROOTC_TO_LITTER']
    unit_list = ['gC m-2', 'gC m-2', 'gC m-2 day-1', 'gC m-2 day-1']
elif varset == 'xfer':
    var_list = ['LEAFC_STORAGE_TO_XFER','LEAFC_XFER','LEAFC_XFER_TO_LEAFC','FROOTC_STORAGE','FROOTC_XFER']
    unit_list = ['gC m-2 day-1', 'gC m--2', 'gC m-2 day-1', 'gC m-2', 'gC m-2']
elif varset == 'soil':
    var_list = ['HR', 'H2OSOI', 'H2OSFC', 'ZWT', 'GROSS_NMIN', 'GROSS_PMIN']
    unit_list = ['gC m-2 day-1', 'mm3/mm3', 'mm', 'm', 'gN m-2 day-1', 'gP m-2 day-1']
elif varset == 'diag':
    # diagnostic variables for the new root scheme

    # Soil growth rates
    # - downregm
    # - h2osfc
    # - t_soisno(c,3) < 0 or not

    # Soil litterfall rates

    # - t_soisno(c,j) averaged
    # - smp_l(c,j) max
    # - h2osoi_liqvol(c,j) averaged

    var_list = ['ONSET_RATE_FROOT', 'COMPS_RATE_FROOT', 'F_NMIN', 'F_H2OSFC', 'BGLFR_FROOT', 'T_D', 'W_D']
    unit_list = ['s-1', 's-1', '', '', 's-1', '', '']
path_out = os.path.join(os.environ['PROJDIR'], 'Phenology_ELM', 'output_elm')
plot_list = [4, 10, 11, 16, 19, 6, 8, 13, 17, 20]


# -------------------------------------------------------------------------------------------------
# PFT level
# -------------------------------------------------------------------------------------------------
pft_list = {'needleleaf evergreen boreal tree': 2, 
            'needleleaf deciduous boreal tree': 3, 
            'broadleaf deciduous boreal shrub': 11,
            'sphagnum moss': 12}
clist = ['#fcbba1', '#fc9272', '#fb6a4a', '#de2d26', '#a50f15', '#c6dbef', '#9ecae1', '#6baed6', '#3182bd', '#08519c']


if varset == 'diag':
    for pft in list(pft_list.keys()):
        pft_id = pft_list[pft]

        fig, axes = plt.subplots(4, 2, figsize = (12, 12), sharex = True)
        for ii, var in enumerate(var_list):
            print(pft, var)

            ax = axes.flat[ii]
            h = [None] * len(plot_list)
            for kk, plot in enumerate(plot_list):
                path_data = os.path.join(os.environ['PROJDIR'], 'E3SM', 'output', f'{date}_plot{plot:02d}_US-SPR_ICB20TRCNPRDCTCBC', 'run')

                if var in ['F_H2OSFC', 'T_D', 'W_D']:
                    flist = sorted(glob(os.path.join(path_data, '*.h0.*.nc')))
                    hr = xr.open_mfdataset(flist, decode_times = False)
                else:
                    flist = sorted(glob(os.path.join(path_data, '*.h2.*.nc')))[:-1]
                    hr = xr.open_mfdataset(flist, decode_times = True)
                    hr = hr.resample(time = '1MS').mean().load()
                tvec = pd.date_range(start = '2015-01-01', end = '2020-12-31', freq='MS')
                hr['time'] = tvec

                if var in ['ONSET_RATE_FROOT', 'COMPS_RATE_FROOT', 'BGLFR_FROOT']:
                    retype = pd.Series(hr[var][:, pft_id].values * 0.64 + hr[var][:, pft_id + 17].values * 0.36, index = hr['time'].to_index())
                elif var == 'F_NMIN':
                    temp = 2 - (1 - (hr['DOWNREG'][:, pft_id].values * 0.64 + hr['DOWNREG'][:, pft_id + 17].values * 0.36))
                    retype = pd.Series(temp, index = hr['time'].to_index())
                elif var == 'F_H2OSFC':
                    col1 = np.clip(1 - hr['H2OSFC'][:, 0].values / 100, a_min = 0.1, a_max = 1.)
                    col2 = np.clip(1 - hr['H2OSFC'][:, 1].values / 100, a_min = 0.1, a_max = 1.)
                    retype = pd.Series(col1 * 0.64 + col2 * 0.36, index = hr['time'].to_index())
                elif var == 'T_D':
                    hr2 = xr.open_mfdataset(sorted(glob(os.path.join(path_data, '*.h2.*.nc')))[:-1], decode_times = True)
                    hr2 = hr2.resample(time = '1MS').mean().load()
                    hr2['time'] = pd.date_range(start = '2015-01-01', end = '2020-12-31', freq='MS')
                    rootfr = hr2['ROOTFR'][:, :, [pft_id, pft_id + 17]].values
                    hr2.close()

                    t_soi_avg = np.zeros([len(tvec), 2])
                    for gr in range(2):
                        t_soi_avg[:, gr] = np.sum(hr['TSOI'].values[:, :, gr] * rootfr[:, :, gr], axis = 1)
                    td = np.clip(np.power(t_soi_avg - 273.15 - 10, 2) / 4 * 0.00175 + 0.1, a_min = None, a_max = 0.5)
                    retype = pd.Series(td[:,0] * 0.64 + td[:,1] * 0.36, index = hr['time'].to_index())
                elif var == 'W_D':
                    hr2 = xr.open_mfdataset(sorted(glob(os.path.join(path_data, '*.h2.*.nc')))[:-1], decode_times = True)
                    hr2 = hr2.resample(time = '1MS').mean().load()
                    hr2['time'] = pd.date_range(start = '2015-01-01', end = '2020-12-31', freq='MS')
                    rootfr = hr2['ROOTFR'][:, :, [pft_id, pft_id + 17]].values
                    hr2.close()

                    hr2 = xr.open_mfdataset(sorted(glob(os.path.join(path_data, '*.h1.*.nc')))[:-1], decode_times = True)
                    smp = hr2['SMP'].resample(time = '1MS').mean().load()
                    smp['time'] = pd.date_range(start = '2015-01-01', end = '2020-12-31', freq='MS')
                    hr2.close()

                    # exclude the non-root zone
                    filt = np.mean(np.mean(rootfr, axis = 0), axis = 1) > 1e-6
                    rootfr = rootfr[:, filt, :]
                    smp = smp[:, filt, :]

                    h2osoi_liq_avg = np.zeros([len(tvec), 2])
                    psi_max = np.zeros([len(tvec), 2]) - 9999
                    psi_t6 = np.zeros([len(tvec), 2])
                    for gr in range(2):
                        h2osoi_liq_avg[:, gr] = np.sum(hr['SOILLIQ'].values[:, filt, gr] / \
                                                       hr['DZSOI'].values[filt, gr].reshape(1, -1) / 1000. \
                                                       * rootfr[:, :, gr], axis = 1)

                        psi_max[:, gr] = np.max(smp.values[:, :, gr], axis = 1)

                        psi_t6[:, gr] = - hr['SUCSAT'][2, gr] * np.power(np.clip(0.6 / hr['WATSAT'].values[2,gr], a_min = 0.01, a_max = None), - hr['BSW'].values[2,gr])

                    wd_le6 = 0.5 + np.arctan(0.05 * np.pi * (np.abs(psi_max) - 35)) / np.pi
                    wd6 = 0.5 + np.arctan(0.05 * np.pi * (np.abs(psi_t6) - 35)) / np.pi
                    wd_gt6 = np.clip(1 - np.log(0.6 / (1.6 - h2osoi_liq_avg)), a_min = None, a_max = 0.9) * (1 - wd6) + wd6

                    wd = np.where(h2osoi_liq_avg <= 0.6, wd_le6, wd_gt6)
                    retype = pd.Series(wd[:,0] * 0.64 + wd[:, 1] * 0.36, index = hr['time'].to_index())

                hr.close()

                htemp, = ax.plot(retype.index, retype.values, '-', lw = 0.5, color = clist[kk])
                h[kk] = htemp

            ax.set_title(var_list[ii] + ' ' + unit_list[ii])

            #ax.set_xlim([0.5, 365.5])
            #ax.set_xticks(range(1,13))
        ax.legend(h, [f'plot{plot:02d}' for plot in plot_list], loc = 'upper right', bbox_to_anchor = (0.5, -.2), ncol = 3)
        fig.savefig(os.path.join(path_out, date, f'ts_bypft_{varset}_{pft_id}.png'), dpi = 600., bbox_inches = 'tight')
        plt.close()

elif varset == 'soil':
    ncols = 2
    nrows = int(np.ceil(len(var_list) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize = (ncols * 3, nrows * 3), sharex = True)
    for ii, var in enumerate(var_list):
        ax = axes.flat[ii]
        h = [None] * len(plot_list)
        for kk, plot in enumerate(plot_list):
            path_data = os.path.join(os.environ['PROJDIR'], 'E3SM', 'output', f'{date}_plot{plot:02d}_US-SPR_ICB20TRCNPRDCTCBC', 'run')
            flist = sorted(glob(os.path.join(path_data, '*.h0.*.nc')))

            # osbolete # autocoversion misidentified PROOT_ONSET_OFFSET as a datetime variable because its unit "days"
            hr = xr.open_mfdataset(flist, decode_times = False)
            units, reference_date = hr.time.attrs['units'].split('since')
            tvec = pd.date_range(start = reference_date, end = '2020-12-31', freq='MS')
            hr['time'] = tvec
            if var == 'H2OSOI':
                # 3rd level (Fortran index starts from 1)
                retype = pd.Series(hr[var][:, 2, 0].values * 0.64 + hr[var][:, 2, 1].values * 0.36, index = hr['time'].to_index())
            else:
                retype = pd.Series(hr[var][:, 0].values * 0.64 + hr[var][:, 1].values * 0.36, index = hr['time'].to_index())
            if ('/s' in hr[var].attrs['units']) and ('day-1' in unit_list[ii]):
                # s to day
                retype = retype * 24 * 3600
            hr.close()

            ax = axes.flat[ii]
            htemp, = ax.plot(retype.index, retype.values, '-', color = clist[kk], lw = 0.5)
            h[kk] = htemp

            ax.set_title(var_list[ii] + ' ' + unit_list[ii])
    ax.legend(h, [f'plot{plot:02d}' for plot in plot_list], loc = 'upper right', bbox_to_anchor = (0.5, -.2), ncol = 3)
    fig.savefig(os.path.join(path_out, date, f'ts_bycol_{varset}.png'), dpi = 600., bbox_inches = 'tight')
    plt.close()

else:
    fig, axes = plt.subplots(len(var_list), 4, figsize = (14, 14), sharex = True)
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
                tvec = [datetime.strptime(x, '%Y-%m-%d %H') for x in hr['time'].indexes['time'].strftime('%Y-%m-%d %H')]
                retype = pd.DataFrame(hr[var][:, pft_list[pft]].values * 0.64 + hr[var][:, pft_list[pft] + 17].values * 0.36, index = tvec)
                if ('/s' in hr[var].attrs['units']) and ('day-1' in unit_list[ii]):
                    # s to day
                    retype = retype * 24 * 3600
                hr.close()

                if var == 'GPP':
                    retype2 = pd.DataFrame(savgol_filter(retype.values, 81, 3, mode = 'constant', cval = np.nan), index = tvec)

                if var == 'QVEGE':
                    retype = retype.iloc[1:]

                ax = axes[ii, jj]
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
    ax.legend(h, [f'plot{plot:02d}' for plot in plot_list], loc = 'upper right', bbox_to_anchor = (0.5, -.2), ncol = 6)
    fig.savefig(os.path.join(path_out, date, f'ts_bypft_{varset}.png'), dpi = 600., bbox_inches = 'tight')
    plt.close()