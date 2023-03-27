""" Following previous studies on evergreen needleleaf phenology, test
    10%, 15%, 25%, 50%, 90% of the LAI amplitude during the green-up
    and green-down phase. 

    The seasonal amplitude is the difference between maximum and minimum baseline.

    Moon et al., 2021. Multiscale assessment of land surface phenology from harmonized Landsat 8 and Sentinel-2, PlanetScope, and PhenoCam imagery. Remote Sensing of Environment.

    Liu et al., 2020. Using the red chromatic coordinate to characterize the phenology of forest canopy photosynthesis. Agricultural and Forest Meteorology.

    Richardson et al., 2018. Intercomparison of phenological transition dates derived from the PhenoCam Dataset V1.0 and MODIS satellite remote sensing. Scientific Reports. 
"""
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import itertools as it
from utils.constants import *
from utils.paths import *
from utils.plotting import *
import xarray as xr
from datetime import datetime
import matplotlib as mpl
from gcc_spruce_visualize import *
from glob import glob
from scipy.signal import savgol_filter
from scipy.stats import linregress


def get_model_dates(variable, pft, read = False):
    if not read:
        pct   = 25 # since GCC 25 was employed for GPP; ['15','25','50','90']
        years = range(2015,2021)

        data = pd.DataFrame(np.nan, 
                            index = pd.MultiIndex.from_product([chamber_list, years]),
                            columns = pd.MultiIndex.from_product([sims_names,
                                                                  ['temperature','CO2','SOS','EOS']]))
        for yy, model in zip(sims_prefix, sims_names):
            for chamber in chamber_list:
                hr = xr.open_mfdataset(os.path.join(path_run, f'{yy}_plot{chamber:02d}_US-SPR_ICB20TRCNPRDCTCBC', 'run', '*.h2.*.nc'), decode_times = True)
                tvec = pd.DatetimeIndex([datetime(year = t.year, month = t.month, day = t.day) for t in hr['time'].values])
                # use savgol filter to smooth the GPP
                filtered = pd.Series(savgol_filter(hr[variable][:, pft].values.copy(), 81, 3),
                                     index = tvec.year)
                hr.close()

                threshold = float(pct) * (filtered.groupby(filtered.index).max() - \
                                          filtered.groupby(filtered.index).min()) / 100 + \
                            filtered.groupby(filtered.index).min()
                filtered = ((filtered - threshold) >= 0)

                start = np.where(np.diff(np.insert(filtered.values, 0, 1).astype(int)) ==  1)[0]
                start = pd.Series(start, index = filtered.index[start])
                start = start.groupby(start.index).agg('first')
                end   = np.where(np.diff(np.insert(filtered.values, 0, 1).astype(int)) == -1)[0]
                end   = pd.Series(end, index = filtered.index[end])
                end   = end.groupby(end.index).agg('last')

                for ss,ee in zip(start, end):
                    year = tvec.year[ss]
                    data.loc[(chamber, year), 
                             (model,'temperature')] = chamber_levels[f'{chamber:02d}'][0]
                    data.loc[(chamber, year),
                             (model,'CO2'        )] = chamber_levels[f'{chamber:02d}'][1]
                    if tvec.dayofyear[ss] > 240:
                        data.loc[(chamber, year+1),
                                 (model,'SOS'      )] = tvec.dayofyear[ss] - 365
                    else:
                        data.loc[(chamber, year),
                                 (model,'SOS'      )] = tvec.dayofyear[ss]

                    year = tvec.year[ee]
                    data.loc[(chamber, year),
                             (model,'temperature')] = chamber_levels[f'{chamber:02d}'][0]
                    data.loc[(chamber, year),
                             (model,'CO2'        )] = chamber_levels[f'{chamber:02d}'][1]
                    if tvec.dayofyear[ee] < 120:
                        data.loc[(chamber, year-1),
                                 (model,'EOS'        )] = tvec.dayofyear[ee] + 365
                    else:
                        data.loc[(chamber, year),
                                 (model,'EOS'        )] = tvec.dayofyear[ee]

        data.to_csv(os.path.join(path_out, 'temperature_sensitivity', 'fit_model_phenology_{}_{}.csv'.format(variable, pft)))
    else:
        data = pd.read_csv(os.path.join(path_out, 'temperature_sensitivity', 'fit_model_phenology_{}_{}.csv'.format(variable, pft)),
                           index_col = [1,2], header = [0,1])
    return data


def get_obs_dates(variable, pft):
    hr = xr.open_dataset(os.path.join(path_intrim, 'spruce_validation_data.nc'))
    if pft == 2:
        pft_name = 'EN'
    elif pft == 3:
        pft_name = 'DN'
    elif pft == 11:
        pft_name = 'SH'
    else:
        pft_name = 'GR'
    if variable == 'TLAI':
        obs_dates = hr['pheno_dates_lai'].loc[:, :, :,
                                              pft_name].transpose('chamber','year','side')
    elif variable == 'GPP':
        sapflow_dates = hr['pheno_dates_sapflow'].loc[:, :, :,
                                                      pft_name].transpose('chamber','year','side')
        gcc_dates = hr['pheno_dates_gcc'].loc[:, :, :,
                                              pft_name].transpose('chamber','year','side')
        obs_dates = (sapflow_dates, gcc_dates)
    else:
        if pft in [2,3]:
            obs_dates = {'SOS': hr['sos_tree'].transpose('chamber','year'), 'EOS': hr['eos_tree'].transpose('chamber','year')}
        else:
            obs_dates = {'SOS': hr['sos_shrub'].transpose('chamber','year'), 'EOS': hr['eos_shrub'].transpose('chamber','year')}
    hr.close()
    return obs_dates


# for slides
read = False
mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.titlesize'] = 16


sims_prefix = ['20221212', '20221231', '20230214']
sims_names = ['Default', 'Evergreen', 'EvgrRoot']
clist = ['#de2d26', '#feb24c', '#3182bd']
co2_levels    = {'ambient' : [6, 20, 13, 8, 17],
                 'elevated': [19, 11, 4, 16, 10]}
pft_list = [2, 3, 11]
pft_list_names = ['Pima', 'Lala', 'Shrub']
for variable in ['TLAI', 'GPP', 'FROOTC']:
    for co2 in ['ambient','elevated']:
        fig, axes = plt.subplots(nrows = 2, ncols = len(pft_list), figsize = (len(pft_list)*4, 6), sharex = True, sharey = False)
        fig.subplots_adjust(hspace = 0.1, wspace = 0.1)        
        for i, side in enumerate(['SOS', 'EOS']):
            for j, (pft, pft_name) in enumerate(zip(pft_list, pft_list_names)):
                ax = axes[i, j]

                model_dates = get_model_dates(variable, pft, read)
                if variable == 'GPP':
                    sapflow_dates, gcc_dates = get_obs_dates(variable, pft)
                elif variable == 'TLAI':
                    obs_dates                = get_obs_dates(variable, pft)
                else:
                    root_dates               = get_obs_dates(variable, pft)

                if variable == 'TLAI':
                    x = np.array([chamber_levels[f'{x:02d}'][0] for x in co2_levels[co2]])
                    x = np.broadcast_to(x.reshape(-1,1), [len(x), obs_dates.shape[1]])
                    y = obs_dates.loc[co2_levels[co2], :, side].values

                    temp = np.isnan(x) | np.isnan(y)
                    x = x[~temp]
                    y = y[~temp]
                    if len(x) > 0:
                        h1 = ax_regress(ax, x, y, display = None,
                                        args_pt = {'color': 'k', 'marker': 'o', 'markersize': 3, 'lw': 0},
                                        args_ln = {'color': 'k', 'lw': 1.2},
                                        args_ci = {'color': 'k', 'alpha': 0.1})

                        res = linregress(x, y)
                        if res.pvalue <= 0.05:
                            fontweight = 'bold'
                        else:
                            fontweight = 'normal'

                        if (i == 1) & (j == 0):
                            ax.text(0.08, 0.02, f'{res.slope:.2f}', fontsize = 8, fontweight = fontweight, transform = ax.transAxes, color = 'k')
                        elif (i == 0) & (j == 0):
                            ax.text(0.08, 0.1, f'{res.slope:.2f}', fontsize = 8, fontweight = fontweight, transform = ax.transAxes, color = 'k')
                        else:
                            ax.text(0.08, 0.92, f'{res.slope:.2f}', fontsize = 8, fontweight = fontweight, transform = ax.transAxes, color = 'k')

                elif variable == 'GPP':
                    x = np.array([chamber_levels[f'{x:02d}'][0] for x in co2_levels[co2]])
                    x = np.broadcast_to(x.reshape(-1,1), [len(x), sapflow_dates.shape[1]])
                    y1 = sapflow_dates.loc[co2_levels[co2], :, side].values
                    y2 = gcc_dates.loc[co2_levels[co2], :, side].values

                    temp = np.isnan(x) | np.isnan(y1)
                    if (~temp).sum() > 0:
                        h1 = ax_regress(ax, x[~temp] - 0.25, y1[~temp], display = None,
                                        args_pt = {'color': 'k', 'marker': 'o', 'markersize': 3, 'lw': 0},
                                        args_ln = {'color': 'k', 'lw': 1.2},
                                        args_ci = {'color': 'k', 'alpha': 0.1})

                        res = linregress(x[~temp], y1[~temp])
                        if res.pvalue <= 0.05:
                            fontweight = 'bold'
                        else:
                            fontweight = 'normal'
                        ax.text(0.08, 0.92, f'{res.slope:.2f}', fontsize = 8, fontweight = fontweight, transform = ax.transAxes, color = 'k')

                    temp = np.isnan(x) | np.isnan(y2)
                    if (~temp).sum() > 0:
                        h2 = ax_regress(ax, x[~temp], y2[~temp], display = None,
                                        args_pt = {'color': '#756bb1', 'marker': 'o', 'markersize': 3, 'lw': 0},
                                        args_ln = {'color': '#756bb1', 'lw': 1.2},
                                        args_ci = {'color': '#756bb1', 'alpha': 0.1})

                        res = linregress(x[~temp], y2[~temp])
                        if res.pvalue <= 0.05:
                            fontweight = 'bold'
                        else:
                            fontweight = 'normal'
                        ax.text(0.08, 0.82, f'{res.slope:.2f}', fontsize = 8, fontweight = fontweight, transform = ax.transAxes, color = '#756bb1')

                else:
                    if side == 'SOS' and co2 == 'elevated':
                        h1, = ax.plot(chamber_levels['10'][0], root_dates[side].loc[10,:].loc[2019], 'o', 
                                      markersize = 3, color = 'k')
                        ax.plot(chamber_levels['19'][0], root_dates[side].loc[19,:].loc[2019], 'o',
                                markersize = 3, color = 'k')
                    else:
                        h1, = ax.plot(chamber_levels['06'][0], root_dates[side].loc[6,:].loc[2015], 'o', 
                                      markersize = 3, color = 'k')

                h3 = [None] * len(sims_names)
                for k, model in enumerate(sims_names):
                    if co2 == 'ambient':
                        x = np.array([chamber_levels[f'{x:02d}'][0] for x in model_dates.index.get_level_values(0) if chamber_levels[f'{x:02d}'][1] == 0]) + 0.25 * (k+1)
                        y = model_dates.loc[[x for x in chamber_list if chamber_levels[f'{x:02d}'][1] == 0], (model, side)].values
                    else:
                        x = np.array([chamber_levels[f'{x:02d}'][0] for x in model_dates.index.get_level_values(0) if chamber_levels[f'{x:02d}'][1] == 500]) + 0.25 * (k+1)
                        y = model_dates.loc[[x for x in chamber_list if chamber_levels[f'{x:02d}'][1] == 500], (model, side)].values
                    temp = np.isnan(x) | np.isnan(y)
                    if (~temp).sum() > 0:
                        h3[k] = ax_regress(ax, x[~temp], y[~temp], display = None,
                                        args_pt = {'color': clist[k], 'marker': 'o', 'markersize': 3, 'lw': 0},
                                        args_ln = {'color': clist[k], 'lw': 1.2},
                                        args_ci = {'color': clist[k], 'alpha': 0.1})

                        res = linregress(x[~temp], y[~temp])
                        if res.pvalue <= 0.05:
                            fontweight = 'bold'
                        else:
                            fontweight = 'normal'

                        if (variable != 'GPP') & (i == 1) & (j == 0):
                            ax.text(0.28 + 0.2 * k, 0.02, f'{res.slope:.2f}', fontsize = 8, fontweight = fontweight, transform = ax.transAxes, color = clist[k])
                        elif (variable == 'TLAI') & (i == 0) & (j == 0):
                            ax.text(0.28 + 0.2 * k, 0.1, f'{res.slope:.2f}', fontsize = 8, fontweight = fontweight, transform = ax.transAxes, color = clist[k])
                        else:
                            ax.text(0.28 + 0.2 * k, 0.92, f'{res.slope:.2f}', fontsize = 8, fontweight = fontweight, transform = ax.transAxes, color = clist[k])

                if j == 0:
                    if side == 'SOS':
                        ax.set_ylabel('SOS')
                    else:
                        ax.set_ylabel('EOS')
                else:
                    ax.set_yticklabels([])

                if i == 0:
                    ax.set_title(pft_name)

                if side == 'SOS':
                    ax.set_ylim([0, 240])
                else:
                    ax.set_ylim([200, 460])
                    ax.set_xlabel(r'$\Delta$ Warming ($^o$C)')

        if variable == 'GPP':
            ax.legend([h1, h2] + h3, ['Sapflow','GCC'] + sims_names, ncol = 5, loc = (-2.2,-0.5))
        else:
            ax.legend([h1] + h3, ['Ground obs'] + sims_names, ncol = 5, loc = (-2.2,-0.5))
        fig.savefig(os.path.join(path_out, 'temperature_sensitivity', 'temperature_sensitivity_{}_{}.png'.format(variable,co2)),
                    dpi = 600., bbox_inches = 'tight')
        plt.close(fig)
