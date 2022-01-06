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


def get_model_dates(variable, read = False):
    if not read:
        pft = 2
        pct_list = ['15','25','50','90']
        years = range(2015,2021)

        data = pd.DataFrame(np.nan, 
                            index = pd.MultiIndex.from_product([pct_list, 
                                                                sorted(chamber_levels.keys()),
                                                                years]),
                            columns = pd.MultiIndex.from_product([['Default', 'Modified'],
                                                        ['temperature','CO2','rising','falling']]))
        for yy, model in zip(['20211008', '20211201'], ['Default', 'Modified']):
            for pct, chamber in it.product(pct_list, chamber_levels.keys()):

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
                    data.loc[(pct, chamber, year), 
                             (model,'temperature')] = chamber_levels[chamber][0]
                    data.loc[(pct, chamber, year),
                             (model,'CO2'        )] = chamber_levels[chamber][1]
                    data.loc[(pct, chamber, year),
                             (model,'rising'      )] = tvec.dayofyear[ss]

                    year = tvec.year[ee]
                    data.loc[(pct, chamber, year),
                             (model,'temperature')] = chamber_levels[chamber][0]
                    data.loc[(pct, chamber, year),
                             (model,'CO2'        )] = chamber_levels[chamber][1]
                    data.loc[(pct, chamber, year),
                             (model,'falling'        )] = tvec.dayofyear[ee]

        data.to_csv(os.path.join(path_out, 'fit_model_phenology_{}.csv'.format(variable)))
    else:
        data = pd.read_csv(os.path.join(path_out, 'fit_model_phenology_{}.csv'.format(variable)),
                           index_col = [0,1,2], header = [0,1])

    return data


def get_obs_dates(variable):
    if variable == 'TLAI':
        data = pd.read_excel(os.path.join(path_input, 'SPRUCE_budburst_summary.xlsx'),
                            engine = 'openpyxl', sheet_name = 'Sheet1', 
                            index_col = [0,1,2,3])
        return data.loc[:, ['buds breaking', 'leaves senescing']]
    else:
        sapflow_dates = parse_sapflow('EN')
        for chamber in chamber_levels.keys():
            for side in ['rising', 'falling']:
                sapflow_dates[int(chamber)][side] = pd.Series( \
                    pd.DatetimeIndex(sapflow_dates[int(chamber)][side]).dayofyear.values, 
                    index = pd.DatetimeIndex(sapflow_dates[int(chamber)][side]).year.values \
                    ).groupby(pd.DatetimeIndex(sapflow_dates[int(chamber)][side]).year.values).mean()
            sapflow_dates[int(chamber)] = pd.DataFrame(sapflow_dates[int(chamber)]).stack()
        sapflow_dates = pd.DataFrame(sapflow_dates)

        path_in = os.path.join(path_data, 'Vegetation', 'PhenoCam_V2_1674', 'data')
        gcc_dates = {}
        for chamber in chamber_levels.keys():
            flist = glob(os.path.join(path_in, 'data_record_4', 
                                      'spruceT%dP%s*_%s_*_3day.csv' % \
                                          (chamber_levels[chamber][0],chamber,'EN')))
            collect = {'rising': None, 'falling': None}
            for count, f in enumerate(flist):
                data = pd.read_csv(f, skiprows = 23, index_col = 0, parse_dates = True)
                temp = find_gcc_dates (data['smooth_gcc_90'])['transition_25']
                for side in ['rising', 'falling']:
                    if count == 0:
                        collect[side] = pd.DataFrame( \
                            {count: pd.DatetimeIndex(temp[side]).dayofyear.values}, 
                             index = pd.DatetimeIndex(temp[side]).year.values)
                    else:
                        collect[side] = pd.concat([collect[side], pd.DataFrame( \
                            {count: pd.DatetimeIndex(temp[side]).dayofyear.values}, 
                             index = pd.DatetimeIndex(temp[side]).year.values)], axis = 1)
            for side in ['rising', 'falling']:
                collect[side] = collect[side].stack()
            gcc_dates[chamber] = pd.DataFrame(collect).stack()
            gcc_dates[chamber].index.names = ['year','count','side']
            gcc_dates[chamber] = gcc_dates[chamber].groupby(['year','side']).mean()
        gcc_dates = pd.DataFrame(gcc_dates)

        return sapflow_dates, gcc_dates


# Plot
mpl.rcParams['font.size'] = 10
mpl.rcParams['axes.titlesize'] = 10
for variable in ['TLAI', 'GPP']:
    model_dates = get_model_dates(variable)
    if variable == 'GPP':
        sapflow_dates, gcc_dates = get_obs_dates(variable)
    else:
        obs_dates                = get_obs_dates(variable)

    co2_levels    = {'ambient': [ 6, 20, 13,  8, 17], 'elevated': [19, 11,  4, 16, 10]}
    thres   = 25 # since GCC 25 was employed for GPP
    fig, axes = plt.subplots(ncols = 2, nrows = 4, figsize = (8, 12), sharex = True,
                             sharey = False)
    count = 0
    for side, co2, model, in it.product(['rising', 'falling'], ['ambient', 'elevated'],
                                        ['Default','Modified']):
        ax = axes.flat[count]

        if variable == 'TLAI':
            temp = obs_dates.loc[np.array([c in co2_levels[co2] \
                                           for c in obs_dates.index.get_level_values('plot')]) & \
                (obs_dates.index.get_level_values('CO2 treatment') == co2), :]
            if side == 'rising':
                h1 = ax_regress(ax, temp.index.get_level_values('Temp treatment'),
                                temp['buds breaking'].values,
                                display = None,
                                args_pt = {'color': 'k', 'marker': 'o', 'markersize': 3, 'lw': 0},
                                args_ln = {'color': 'k', 'lw': 1.2},
                                args_ci = {'color': 'k', 'alpha': 0.1})
            else:
                h1 = ax_regress(ax, temp.index.get_level_values('Temp treatment'),
                                temp['leaves senescing'].values,
                                display = None,
                                args_pt = {'color': 'k', 'marker': 'o', 'markersize': 3, 'lw': 0},
                                args_ln = {'color': 'k', 'lw': 1.2},
                                args_ci = {'color': 'k', 'alpha': 0.1})
        else:
            temp = sapflow_dates.loc[(slice(None), side), co2_levels[co2]]
            x = np.broadcast_to([chamber_levels['%02d' % c][0] for c in temp.columns], temp.shape)
            h1 = ax_regress(ax, x.reshape(-1), temp.values.reshape(-1),
                            display = None,
                            args_pt = {'color': 'k', 'marker': 'o', 'markersize': 3, 'lw': 0},
                            args_ln = {'color': 'k', 'lw': 1.2},
                            args_ci = {'color': 'k', 'alpha': 0.1})

            temp = gcc_dates.loc[(slice(None), side), ['%02d' % c for c in co2_levels[co2]]]
            x = np.broadcast_to([chamber_levels[c][0] for c in temp.columns], temp.shape)
            h2 = ax_regress(ax, x.reshape(-1), temp.values.reshape(-1),
                            display = None,
                            args_pt = {'color': '#756bb1', 'marker': 'o', 'markersize': 3,'lw': 0},
                            args_ln = {'color': '#756bb1', 'lw': 1.2},
                            args_ci = {'color': '#756bb1', 'alpha': 0.1})

        temp2 = model_dates.loc[ \
            (model_dates.index.get_level_values(0) == str(thres)) & \
             np.array([int(c) in co2_levels[co2] for c in model_dates.index.get_level_values(1)]),
                                 model]
        h3 = ax_regress(ax, temp2['temperature'].values, temp2[side].values,
                        display = None,
                        args_pt = {'color': 'g', 'marker': 'o', 'markersize': 3, 'lw': 0},
                        args_ln = {'color': 'g', 'lw': 1.2},
                        args_ci = {'color': 'g', 'alpha': 0.1})

        if np.mod(count, 2) == 0:
            ax.set_ylabel(side + ' ' + co2)
        if count < 2:
            ax.set_title(model)
        if side == 'rising':
            ax.set_ylim([0, 190])
        else:
            ax.set_ylim([-220, 460])

        count += 1
    if variable == 'GPP':
        ax.legend([h1, h2, h3], ['Sapflow','GCC','Mod'], loc = [-0.5, -0.3], ncol = 3)
    else:
        ax.legend([h1, h3], ['Obs','Mod'], loc = [-0.5, -0.3], ncol = 2)
    fig.savefig(os.path.join(path_out, 'temperature_sensitivity_' + variable + '.png'),
                dpi = 600., bbox_inches = 'tight')
    plt.close(fig)
