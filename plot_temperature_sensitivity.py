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


def get_model_dates(variable, pft, read = False):
    if not read:
        pct   = 25 # since GCC 25 was employed for GPP; ['15','25','50','90']
        years = range(2015,2021)

        data = pd.DataFrame(np.nan, 
                            index = pd.MultiIndex.from_product([chamber_list, years]),
                            columns = pd.MultiIndex.from_product([sims_names,
                                                                  ['temperature','CO2','SOS',
                                                                   'EOS']]))
        for yy, model in zip(sims_prefix, sims_names):
            for chamber in chamber_list:
                hr = xr.open_mfdataset(os.path.join(path_run, yy + '_plot' + chamber + \
                    '_US-SPR_ICB20TRCNPRDCTCBC', 'run', '*.h2.*.nc'), decode_times = True)
                tvec = pd.DatetimeIndex([datetime(year = t.year, month = t.month,
                                                  day = t.day) for t in hr['time'].values])
                # use savgol filter to smooth the GPP
                filtered = pd.Series(savgol_filter(hr[variable][:, pft].values.copy(), 81, 3),
                                     index = tvec.year)
                hr.close()

                threshold = float(pct) * (filtered.groupby(filtered.index).max() - \
                                          filtered.groupby(filtered.index).min()) / 100 + \
                            filtered.groupby(filtered.index).min()
                filtered = ((filtered - threshold) > 0)

                start = np.where(np.diff(np.insert(filtered.values, 0, 1).astype(int)) ==  1)[0]
                start = pd.Series(start, index = filtered.index[start])
                start = start.groupby(start.index).agg('first')
                end   = np.where(np.diff(np.insert(filtered.values, 0, 1).astype(int)) == -1)[0]
                end   = pd.Series(end, index = filtered.index[end])
                end   = end.groupby(end.index).agg('last')

                for ss,ee in zip(start, end):
                    year = tvec.year[ss]
                    data.loc[(chamber, year), 
                             (model,'temperature')] = chamber_levels[chamber][0]
                    data.loc[(chamber, year),
                             (model,'CO2'        )] = chamber_levels[chamber][1]
                    if tvec.dayofyear[ss] > 240:
                        data.loc[(chamber, year+1),
                                 (model,'SOS'      )] = tvec.dayofyear[ss] - 365
                    else:
                        data.loc[(chamber, year),
                                 (model,'SOS'      )] = tvec.dayofyear[ss]

                    year = tvec.year[ee]
                    data.loc[(chamber, year),
                             (model,'temperature')] = chamber_levels[chamber][0]
                    data.loc[(chamber, year),
                             (model,'CO2'        )] = chamber_levels[chamber][1]
                    if tvec.dayofyear[ee] < 120:
                        data.loc[(chamber, year-1),
                                 (model,'EOS'        )] = tvec.dayofyear[ee] + 365
                    else:
                        data.loc[(chamber, year),
                                 (model,'EOS'        )] = tvec.dayofyear[ee]

        data.to_csv(os.path.join(path_out,
                                 'fit_model_phenology_{}_{}.csv'.format(variable, pft)))
    else:
        data = pd.read_csv(os.path.join(path_out,
                                        'fit_model_phenology_{}_{}.csv'.format(variable, pft)),
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
            obs_dates = hr['sos_tree'].transpose('chamber','year')
        else:
            obs_dates = hr['sos_shrub'].transpose('chamber','year')
    hr.close()
    return obs_dates


# for slides
read = False
mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.titlesize'] = 16


pft = 11
if pft == 2:
    sims_prefix = ['20211008', '20220103', '20220407_stopdate171_switch1']
    sims_names = ['Default', 'AltEvgrPheno', 'AltRoot_171_1']
elif pft == 3:
    sims_prefix = ['20211008', '20220103', '20220407_stopdate171_switch1']
    sims_names = ['Default', 'AltEvgrPheno', 'AltRoot_171_1']
elif pft == 11:
    sims_prefix = ['20211008', '20220103', '20220407_stopdate171_switch1']
    sims_names = ['Default', 'AltEvgrPheno', 'AltRoot_171_1']
co2_levels    = {'ambient' : [6, 20, 13, 8, 17],
                 'elevated': [19, 11, 4, 16, 10]}
for variable in ['TLAI', 'GPP', 'FROOTC']:
    model_dates = get_model_dates(variable, pft, read)
    if variable == 'GPP':
        sapflow_dates, gcc_dates = get_obs_dates(variable, pft)
    elif variable == 'TLAI':
        obs_dates                = get_obs_dates(variable, pft)
    else:
        root_dates               = get_obs_dates(variable, pft)

    for co2 in ['ambient','elevated']:
        fig, axes = plt.subplots(ncols = len(sims_names), nrows = 2, 
                                 figsize = (len(sims_names)*4, 6), sharex = True,
                                 sharey = False)
        fig.subplots_adjust(hspace = 0.1, wspace = 0.1)
        count = 0
        for side, model in it.product(['SOS','EOS'], sims_names):
            ax = axes.flat[count]

            if variable == 'TLAI':
                x = np.array([chamber_levels[f'{x:02d}'][0] for x in co2_levels[co2]])
                x = np.broadcast_to(x.reshape(-1,1), [len(x), obs_dates.shape[1]])
                y = obs_dates.loc[co2_levels[co2], :, side].values

                temp = np.isnan(x) | np.isnan(y)
                x = x[~temp]
                y = y[~temp]
                if len(x) > 0:
                    h1 = ax_regress(ax, x, y, display = None,
                                    args_pt = {'color': 'k', 'marker': 'o', 'markersize': 6,
                                               'lw': 0},
                                    args_ln = {'color': 'k', 'lw': 1.2},
                                    args_ci = {'color': 'k', 'alpha': 0.1})
            elif variable == 'GPP':
                x = np.array([chamber_levels[f'{x:02d}'][0] for x in co2_levels[co2]])
                x = np.broadcast_to(x.reshape(-1,1), [len(x), sapflow_dates.shape[1]])
                y1 = sapflow_dates.loc[co2_levels[co2], :, side].values
                y2 = gcc_dates.loc[co2_levels[co2], :, side].values

                temp = np.isnan(x) | np.isnan(y1)
                if (~temp).sum() > 0:
                    h1 = ax_regress(ax, x[~temp], y1[~temp], display = None,
                                    args_pt = {'color': 'k', 'marker': 'o', 'markersize': 6,
                                               'lw': 0},
                                    args_ln = {'color': 'k', 'lw': 1.2},
                                    args_ci = {'color': 'k', 'alpha': 0.1})

                temp = np.isnan(x) | np.isnan(y2)
                if (~temp).sum() > 0:
                    h2 = ax_regress(ax, x[~temp], y2[~temp], display = None,
                                    args_pt = {'color': '#756bb1', 'marker': 'o',
                                               'markersize': 6,'lw': 0},
                                    args_ln = {'color': '#756bb1', 'lw': 1.2},
                                    args_ci = {'color': '#756bb1', 'alpha': 0.1})
            else:
                if side == 'SOS' and co2 == 'elevated':
                    h1, = ax.plot(chamber_levels['10'][0], root_dates.loc[10,:].loc[2019], 'o', 
                                  markersize = 6, color = 'k')
                    ax.plot(chamber_levels['19'][0], root_dates.loc[19,:].loc[2019], 'o',
                            markersize = 6, color = 'k')

            x = np.array([chamber_levels[x][0] for x in model_dates.index.get_level_values(0)])
            y = model_dates.loc[:, (model, side)].values
            temp = np.isnan(x) | np.isnan(y)
            h3 = ax_regress(ax, x[~temp], y[~temp], display = None,
                            args_pt = {'color': 'g', 'marker': 'o', 'markersize': 6, 'lw': 0},
                            args_ln = {'color': 'g', 'lw': 1.2},
                            args_ci = {'color': 'g', 'alpha': 0.1})

            if np.mod(count, len(sims_names)) == 0:
                if side == 'SOS':
                    ax.set_ylabel('SOS')
                else:
                    ax.set_ylabel('EOS')
            else:
                ax.set_yticklabels([])
            if count < len(sims_names):
                ax.set_title(model)
            if side == 'SOS':
                ax.set_ylim([0, 240])
            else:
                ax.set_ylim([200, 460])
                ax.set_xlabel(r'$\Delta$ Warming ($^o$C)')

            count += 1
        if variable == 'GPP':
            if pft == 3:
                ax.legend([h3], ['Modeled'], ncol = 3, loc = (-1,-0.5))
            else:
                ax.legend([h1, h2, h3], ['Sapflow','GCC','Modeled'], ncol = 3, loc = (-1,-0.5))
        else:
            if pft == 3 and variable == 'TLAI':
                ax.legend([h3], ['Modeled'], ncol =2, loc = (-1,-0.5))
            else:
                ax.legend([h1, h3], ['Observed','Modeled'], ncol =2, loc = (-1,-0.5))
        fig.savefig(os.path.join(path_out, 
                                 'temperature_sensitivity_{}_{}_{}.png'.format(variable,co2,pft)),
                    dpi = 600., bbox_inches = 'tight')
        plt.close(fig)
