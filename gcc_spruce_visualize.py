import pandas as pd
import os
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
import itertools as it
from utils.paths import *
from utils.plotting import *
from datetime import datetime, timedelta
import xarray as xr


####################################################################################################
# Functions
####################################################################################################
def parse_sapflow(pft):
    if pft == 'DN':
        species = 'Lala'
    elif pft == 'EN':
        species = 'Pima'
    raw = pd.read_excel(os.path.join(path_input, 'sapflow timing.xlsx'),
                        sheet_name = 'sapflow timing updated', engine = 'openpyxl')
    summary = {}
    for plots in np.unique(raw['plot']):
        summary[plots] = {'rising': [], 'falling': []}

        temp = raw.loc[(raw['plot'] == plots) & (raw['Species'] == species), :]
        for i in range(temp.shape[0]):
            for y in range(2016, 2020):
                if ~np.isnan(temp['%s Begins Day #' % y].iloc[i]):
                    summary[plots]['rising'].append( \
                        datetime.strptime('%s-12-31' % (y-1), '%Y-%m-%d') + \
                        timedelta(days = temp['%s Begins Day #' % y].iloc[i]))

                if ~np.isnan(temp['%s Ends Day #' % y].iloc[i]):
                    summary[plots]['falling'].append( \
                        datetime.strptime('%s-12-31' % (y-1), '%Y-%m-%d') + \
                        timedelta(days = temp['%s Ends Day #' % y].iloc[i]))

    return summary


def parse_gcc_dates(sitename, roi):
    gcc_dates = pd.read_csv(os.path.join(path_in, 'data_record_5', '3day', 
                                            'spruce%s_EN_%s_3day_transition_dates.csv' % (sitename,roi)),
                            skiprows = 16)
    summary = {}
    for th in ['transition_10', 'transition_25', 'transition_50']:
        summary[th] = {'rising': [], 'falling': []}
        for dir in ['rising', 'falling']:
            temp = gcc_dates.loc[(gcc_dates['direction'] == dir) & \
                                 (gcc_dates['gcc_value'] == 'gcc_90'), th]
            summary[th][dir] = [datetime.strptime(i, '%Y-%m-%d') for i in temp]
    return summary


def find_gcc_dates(gcc_90):
    summary = {}
    for p in [10, 25, 50]:
        summary['transition_%d' % p] = {'rising': [], 'falling': []}
        thres = gcc_90.groupby(gcc_90.index.year).max() * p/100 + \
                gcc_90.groupby(gcc_90.index.year).min() * (1-p/100)
        for y in np.unique(gcc_90.index.year):
            temp = gcc_90.loc[gcc_90.index.year == y] >= thres.loc[y]
            if temp.index[0].month < 3:
                summary['transition_%d' % p]['rising'].append(temp.index[temp].min())
            if temp.index[-1].month > 9:
                summary['transition_%d' % p]['falling'].append(temp.index[temp].max())
    return summary


def get_daymetT(file, sitename):
    """ The average of T_{max} and T_{min} must be above 5 degC for RCC-based start & end dates

    Liu et al. 2020. Using the red chromatic coordinate to characterize the phenology of forest canopy photosynthesis. Agricultural and Forest Meteorology. 285-286. 107910.
    """
    f = open(file)
    line = f.readline()
    while not 'lat' in line.lower():
        line = f.readline()
    lat = float(line.split('Lat: ')[-1])
    while not 'lon' in line.lower():
        line = f.readline()
    lon = float(line.split('Lon: ')[-1])
    f.close()

    hr = xr.open_mfdataset(os.path.join(path_intrim, 'Daymet', 'spruce%s' % sitename,
                                        'tmax_20[12]*.nc'))
    mat = np.power(hr['lat'].values - lat, 2) + np.power(hr['lon'].values - lon, 2)
    index = np.unravel_index(np.argmin(mat), mat.shape)
    tmax = hr['tmax'][:, index[0], index[1]].values
    hr.close()

    hr = xr.open_mfdataset(os.path.join(path_intrim, 'Daymet', 'spruce%s' % sitename,
                                        'tmin_20[12]*.nc'))
    tmin = hr['tmin'][:, index[0], index[1]].values
    hr.close()

    tindex = pd.Index(datetime.strptime('2010-01-01','%Y-%m-%d') + timedelta(days = i) \
        for i in range((datetime.strptime('2020-12-31','%Y-%m-%d') - \
                        datetime.strptime('2010-01-01','%Y-%m-%d')).days + 1))
    ta     = pd.Series(np.nan, tindex)

    feb29  = (tindex.month != 2) | (tindex.day != 29)
    ta.loc[tindex[ feb29]] = (tmax + tmin) / 2
    ta.loc[tindex[~feb29]] = (ta.loc[tindex[~feb29] - timedelta(days = 1)] + \
                              ta.loc[tindex[~feb29] + timedelta(days = 1)]) / 2
    return ta


def find_rcc_dates(rcc_series, file, sitename, roi, pft, overwrite = False):
    ta_file = os.path.join(path_intrim, 'Daymet', 'spruce%s' % sitename, 
                           'ta_extracted_%s.csv' % roi)
    if os.path.exists(ta_file) and not overwrite:
        ta = pd.read_csv(ta_file, index_col = 0, parse_dates = True).iloc[:,0]
    else:
        ta = get_daymetT(file, sitename)
        ta.to_csv(ta_file)
    ta = ta.loc[rcc_series.index]
    summary  = {'rising': [], 'falling': []}

    if pft == 'EN':
        for y in np.unique(rcc_series.index.year):
            rising  = rcc_series.loc[(rcc_series.index.year == y) & \
                                    (rcc_series.index.dayofyear < 150) & \
                                    (ta >= 0)]
            if len(rising) > 0:
                summary['rising' ].append(rising.index[np.argmax(rising)])
            falling  = rcc_series.loc[(rcc_series.index.year == y) & \
                                    (rcc_series.index.dayofyear > 215) & (ta <= 5)]
            if len(falling) > 0:
                localmax = np.where((falling >= np.insert(falling.values[:-1], 0, 0)) & \
                                    (falling >= np.append(falling.values[1:] , 0)))[0]
                localmin = np.where((falling <= np.insert(falling.values[:-1], 0, 1)) & \
                                    (falling <= np.append(falling.values[1:] , 1)))[0]
                firstmin = np.where(localmin > localmax[0])[0]
                if len(firstmin) > 0:
                    summary['falling'].append(falling.index[localmin[firstmin[0]]])
                else:
                    summary['falling'].append(falling.index[localmax[0]])
    else:
        for y in np.unique(rcc_series.index.year):
            rising  = rcc_series.loc[(rcc_series.index.year == y) & \
                                    (rcc_series.index.dayofyear < 150) & \
                                    (ta >= 5)]
            if len(rising) > 0:
                summary['rising' ].append(rising.index[np.argmax(rising)])
            falling  = rcc_series.loc[(rcc_series.index.year == y) & \
                                      (rcc_series.index.dayofyear > 215) & (ta <= 5)]
            if len(falling) > 0:
                summary['falling' ].append(falling.index[np.argmax(falling)])
    return summary


if __name__ == '__main__':
    ################################################################################################
    # Constants
    ################################################################################################
    pft = 'DN' # 'EN', 'DN', 'SH'
    tlevels = [0, 2.25, 4.5, 6.75, 9]
    aco2    = [ 6, 20, 13,  8, 17]
    eco2    = [19, 11,  4, 16, 10]

    path_in = os.path.join(path_data, 'Vegetation', 'PhenoCam_V2_1674', 'data')

    snow_ground = pd.read_csv(os.path.join(path_input,
                              '2015_20200331_Spruce_Snow_On_Ground_Flags.csv'),
                              index_col = 0, parse_dates = True)
    snow_trees  = pd.read_csv(os.path.join(path_input, 
                              '2015_20200331_Spruce_Snow_On_Trees_Flags.csv'),
                              index_col = 0, parse_dates = True)
    sapflow_dates = parse_sapflow(pft)


    ################################################################################################
    # Time series
    ################################################################################################
    for t,c in it.product(tlevels, aco2 + eco2):
        flist = glob(os.path.join(path_in, 'data_record_4',
                                  'spruceT%dP%02d*_%s_*_3day.csv' % (t,c,pft)))

        for f in flist:
            data = pd.read_csv(f, skiprows = 23, index_col = 0, parse_dates = True)

            sitename  = f.split('_')[-4].split('spruce')[-1]
            roi       = f.split('_')[-2]
            #gcc_dates = parse_gcc_dates(sitename, roi)
            gcc_dates = find_gcc_dates (data['smooth_gcc_90'])
            rcc_dates = find_rcc_dates (data['smooth_rcc_90'], f, sitename, roi, pft, False)

            #gcc1  = data['g_mean'] / (data['r_mean'] + data['g_mean'] + data['b_mean'])
            #rcc1  = data['r_mean'] / (data['r_mean'] + data['g_mean'] + data['b_mean'])
            #grvi1 = (data['g_mean'] - data['r_mean']) / (data['g_mean'] + data['r_mean'])

            fig, ax = plt.subplots(figsize = (8,5))
            #h1, = ax.plot(data.index,  gcc1, 'og', lw = 0, markersize = 0.5)
            #h2, = ax.plot(data.index,  rcc1, 'or', lw = 0, markersize = 0.5)

            #ax2 = ax.twinx()
            #ax2.yaxis.label.set_color('#8856a7')
            #ax2.tick_params(axis = 'y', colors = '#8856a7')
            #h3, = ax2.plot(data.index, grvi1,  'o', color = '#8856a7', lw = 0, markersize = 0.5)

            h4, = ax.plot(data.index[~np.isnan(data['gcc_90'])], 
                        data['gcc_90'][~np.isnan(data['gcc_90'])], ':g', lw = 0.5)
            h5, = ax.plot(data.index, data['smooth_gcc_90'], '-g', lw = 0.5)
            h6, = ax.plot(data.index[~np.isnan(data['rcc_90'])],
                        data['rcc_90'][~np.isnan(data['rcc_90'])], ':r', lw = 0.5)
            h7, = ax.plot(data.index, data['smooth_rcc_90'], '-r', lw = 0.5)

            for dates in list(snow_ground.index[snow_ground[sitename] > 0]) + \
                        list(snow_trees .index[snow_trees [sitename] > 0]):
                h8 = ax.axvspan(dates - timedelta(days = 0.5), dates + timedelta(days = 0.5),
                                facecolor = '#00FFFF', edgecolor = None)
        
            for dates in sapflow_dates[c]['rising'] + sapflow_dates[c]['falling']:
                h9 = ax.axvline(dates, color = 'k', lw = 0.5)

            h10 = [None, None, None]
            for i, th, mk in zip(range(3),
                                ['transition_10', 'transition_25', 'transition_50'],
                                ['o', '^', 'v']):
                h10[i], = ax.plot(gcc_dates[th]['rising'], [0.36] * len(gcc_dates[th]['rising']),
                                mk, markersize = 2, color = 'g')
                ax.plot(gcc_dates[th]['falling'], [0.36] * len(gcc_dates[th]['falling']), 
                        mk, markersize = 2, color = 'g')

            h11, = ax.plot(rcc_dates['rising'], [0.42] * len(rcc_dates['rising']),
                        'o', markersize = 2, color = 'r')
            ax.plot(rcc_dates['falling'], [0.42] * len(rcc_dates['falling']),
                    'o', markersize = 2, color = 'r')

            ax.legend([h4,h5,h6,h7,h8,h9] + h10 + [h11],
                    ['gcc_90','smooth_gcc_90','rcc_90','smooth_rcc_90','snow on ground or trees',
                    'sapflow start/end dates'] + \
                    ['gcc_90_transition_%d' % i for i in [10,25,50]] + \
                    ['rcc_90_transition'],
                    ncol = 3, loc = (0., 1.1))

            ax.set_xlim([data.index[0], data.index[-1]])

            fig.savefig(os.path.join(path_out, pft, '%s_%s.png' % (sitename,roi)), 
                        dpi = 600, bbox_inches = 'tight')
            plt.close(fig)

    ################################################################################################
    # Relationships between temperature/CO2 and transition dates 
    ################################################################################################
    dayofyear = pd.DataFrame('                       ', index = [],
                            columns = ['tlevel', 'co2', 'direction', 'year', 'roi',   
                                        'sapflow',                 
                                        'gcc_90_transition_10', 'gcc_90_transition_25', 
                                        'gcc_90_transition_50', 'rcc_90_transition'])
    count = 0
    for t, (c,co2) in it.product(tlevels, zip(aco2 + eco2, ['a']*len(aco2) + ['e']*len(eco2))):
        flist = glob(os.path.join(path_in, 'data_record_4', 'spruceT%dP%02d*_%s_*_3day.csv' % (t,c,pft)))
        for f in flist:
            data = pd.read_csv(f, skiprows = 23, index_col = 0, parse_dates = True)

            sitename  = f.split('_')[-4].split('spruce')[-1]
            roi       = f.split('_')[-2]
            #gcc_dates = parse_gcc_dates(sitename, roi)
            gcc_dates = find_gcc_dates (data['smooth_gcc_90'])
            rcc_dates = find_rcc_dates (data['smooth_rcc_90'], f, sitename, roi, pft, False)

            for y,dir in it.product(range(2015, 2020), ['rising', 'falling']):
                dayofyear.loc[count, 'tlevel']    = '%s' % t
                dayofyear.loc[count, 'co2']       = co2
                dayofyear.loc[count, 'direction'] = dir
                dayofyear.loc[count, 'roi']       = roi
                dayofyear.loc[count, 'year']      = str(y)

                if len(sapflow_dates[c][dir]) > 0:
                    n = np.where(pd.Index(sapflow_dates[c][dir]).year == y)[0]
                    if len(n) > 0:
                        dayofyear.loc[count, 'sapflow'] = \
                            sapflow_dates[c][dir][n[0]].strftime('%Y-%m-%d')
                for s in [10,25,50]:
                    n = np.where(pd.Index(gcc_dates['transition_%s' % s][dir]).year == y)[0]
                    if len(n) > 0:
                        dayofyear.loc[count, 'gcc_90_transition_%s' % s] = \
                            gcc_dates['transition_%s' % s][dir][n[0]].strftime('%Y-%m-%d')
                n = np.where(pd.Index(rcc_dates[dir]).year == y)[0]
                if len(n) > 0:
                    dayofyear.loc[count, 'rcc_90_transition'] = \
                        rcc_dates[dir][n[0]].strftime('%Y-%m-%d')

                count += 1
    dayofyear.to_csv(os.path.join(path_out, pft, 'summary_dayofyear.csv'), index = False)


    dayofyear = pd.read_csv(os.path.join(path_out, pft, 'summary_dayofyear.csv'))
    fig, axes = plt.subplots(ncols = 2, nrows = 2, figsize = (8, 8), sharex = True, sharey = False)
    for i,dir in enumerate(['rising', 'falling']):
        for j,co2 in enumerate(['a','e']):
            temp = dayofyear.loc[(dayofyear['direction'] == dir) & (dayofyear['co2'] == co2), :]

            np.random.seed(42)
            ax  = axes[i,j]
            h1, = ax.plot(temp['tlevel'] + np.random.rand(temp.shape[0]) * 0.5, 
                        pd.DatetimeIndex(temp['sapflow']).dayofyear, 
                        'ok', markersize = 3)
            ax_regress(ax, temp['tlevel'].values, pd.DatetimeIndex(temp['sapflow']).dayofyear.values, 
                    display = None, 
                    args_pt = {'markersize': 0, 'lw': 0},
                    args_ln = {'color': 'k', 'lw': 1.2},
                    args_ci = {'color': 'k', 'alpha': 0.1})

            h2  = [None] * 3
            for k, th, col in zip(range(3),
                                ['transition_10', 'transition_25', 'transition_50'],
                                ['#4daf4a', '#377eb8', '#984ea3']):
                h2[k], = ax.plot(temp['tlevel'] + np.random.rand(temp.shape[0]) * 0.5, 
                                pd.DatetimeIndex(temp['gcc_90_%s' % th]).dayofyear, 
                                'o', color = col, markersize = 3)
                ax_regress(ax, temp['tlevel'].values,
                        pd.DatetimeIndex(temp['gcc_90_%s' % th]).dayofyear.values,
                        display = None,
                        args_pt = {'markersize': 0, 'lw': 0},
                        args_ln = {'color': col, 'lw': 1.2},
                        args_ci = {'color': col, 'alpha': 0.1})

            h3, = ax.plot(temp['tlevel'] + np.random.rand(temp.shape[0]) * 0.5,
                        pd.DatetimeIndex(temp['rcc_90_transition']).dayofyear,
                        'or', markersize = 3)
            ax_regress(ax, temp['tlevel'].values,
                        pd.DatetimeIndex(temp['rcc_90_transition']).dayofyear.values,
                        display = None,
                        args_pt = {'markersize': 0, 'lw': 0},
                        args_ln = {'color': 'r', 'lw': 1.2},
                        args_ci = {'color': 'r', 'alpha': 0.1})

            if j == 0:
                ax.set_ylabel(dir)
            if i == 0:
                ax.set_title(co2 + 'co2')
                ax.set_ylim([0, 180])
            else:
                ax.set_xlabel('$\Delta$ temperature')
                ax.set_ylim([250,366])
    ax.legend([h1] + h2 + [h3], ['sapflow','gcc_90_transition_10', 'gcc_90_transition_25', 
                                'gcc_90_transition_50', 'rcc_90_transition'],
            ncol = 3, loc = (-1.2, -0.4))
    fig.savefig(os.path.join(path_out, pft, 'summary_dayofyear.png'), dpi = 600., bbox_inches = 'tight')
    plt.close(fig)

    ################################################################################################
    # Annual cycle
    ################################################################################################
    for t,c in it.product(tlevels, aco2 + eco2):
        flist = glob(os.path.join(path_in, 'data_record_4', 'spruceT%dP%02d*_%s_*_3day.csv' % (t,c,pft)))
        for f in flist:
            data = pd.read_csv(f, skiprows = 23, index_col = 0, parse_dates = True)
            # drop Feb 29
            data = data.loc[(data.index.month != 2) | (data.index.day != 29), :]

            sitename  = f.split('_')[-4].split('spruce')[-1]
            roi       = f.split('_')[-2]
            #gcc_dates = parse_gcc_dates(sitename, roi)
            gcc_dates = find_gcc_dates (data['smooth_gcc_90'])
            rcc_dates = find_rcc_dates (data['smooth_rcc_90'], f, sitename, roi, pft, False)

            fig, ax = plt.subplots(figsize = (8,5))

            adjusted = pd.DatetimeIndex([datetime(year = 1999, month = i.month, day = i.day) \
                                        for i in data.index]).dayofyear

            gcc_90 = data['smooth_gcc_90'].groupby(adjusted).mean()
            rcc_90 = data['smooth_rcc_90'].groupby(adjusted).mean()

            h5, = ax.plot(gcc_90.index, gcc_90, '-g', lw = 0.5)
            h6, = ax.plot(rcc_90.index, rcc_90, '-r', lw = 0.5)

            snow_dates = pd.Index(list(snow_ground.index[snow_ground[sitename] > 0]) + \
                                list(snow_trees .index[snow_trees [sitename] > 0]))
            snow_bin   = pd.Series(0, index = pd.date_range(snow_dates.min(), snow_dates.max()))
            snow_bin.loc[snow_dates] = 1
            # drop Feb 29
            snow_bin   = snow_bin[(snow_bin.index.month != 2) | (snow_bin.index.day != 29)]
            adjusted   = pd.DatetimeIndex([datetime(year = 1999, month = i.month, day = i.day) \
                                        for i in snow_bin.index]).dayofyear
            snow_pct   = snow_bin.groupby(adjusted).mean()
            ax2        = ax.twinx()
            ax2.yaxis.label.set_color('b')
            ax2.tick_params(axis = 'y', colors = 'b')
            for index, count in snow_pct.iteritems():
                h8 = ax2.bar(index, count, width = 1, facecolor = '#00FFFF', edgecolor = None,
                            alpha = 0.2)

            adjusted = []
            for i in sapflow_dates[c]['rising'] + sapflow_dates[c]['falling']:
                if (i.month == 2) & (i.day == 29):
                    adjusted.append(datetime(year = 1999, month = 2, day = 28))
                else:
                    adjusted.append(datetime(year = 1999, month = i.month, day = i.day))
            for a in pd.Index(adjusted).dayofyear:
                h9 = ax.axvline(a, color = 'k', lw = 0.5)

            h10 = [None, None, None]
            for j, th, mk in zip(range(3),
                                ['transition_10', 'transition_25', 'transition_50'],
                                ['o', '^', 'v']):
                adjusted = []
                for i in gcc_dates[th]['rising'] + gcc_dates[th]['falling']:
                    if (i.month == 2) & (i.day == 29):
                        adjusted.append(datetime(year = 1999, month = 2, day = 28))
                    else:
                        adjusted.append(datetime(year = 1999, month = i.month, day = i.day))
                h10[j], = ax.plot(pd.Index(adjusted).dayofyear, [0.36] * len(adjusted),
                                mk, markersize = 2, color = 'g')

            adjusted = []
            for i in rcc_dates['rising'] + rcc_dates['falling']:
                if (i.month == 2) & (i.day == 29):
                    adjusted.append(datetime(year = 1999, month = 2, day = 28))
                else:
                    adjusted.append(datetime(year = 1999, month = i.month, day = i.day))
            h11, = ax.plot(pd.Index(adjusted).dayofyear,
                        [0.42] * len(adjusted), 'o', markersize = 2, color = 'r')

            ax.legend([h5,h6,h8,h9] + h10 + [h11],
                    ['mean smooth_gcc_90','mean smooth_rcc_90','snow on ground or trees',
                    'sapflow start/end dates'] + \
                    ['gcc_90_transition_%d' % i for i in [10,25,50]] + \
                    ['rcc_90_transition'],
                    ncol = 3, loc = (0., 1.1))

            ax.set_xlim([0.5,365.5])
            ax2.set_ylabel('Fraction of snow days')

            fig.savefig(os.path.join(path_out, pft, 'cycle_%s_%s.png' % (sitename,roi)), 
                        dpi = 600, bbox_inches = 'tight')
            plt.close(fig)
