import pandas as pd
import os
from pandas.core.indexes.datetimes import DatetimeIndex
from utils.constants import *
from utils.paths import *
from utils.analysis import *
import numpy as np
from datetime import datetime, timedelta
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib as mpl
import itertools as it


lat, lon, tzn, elv = '47.563 N', '93.438 W', 'US/Central', 430


##################################################################################################
# Growing degree days model for budburst
##################################################################################################
pheno_obs = pd.read_excel(os.path.join(path_input, 'SPRUCE_budburst_summary.xlsx'),
                          engine = 'openpyxl', sheet_name = 'Sheet1')
pheno_obs = pheno_obs.loc[:, 
    ['year', 'plot', 'buds breaking']].set_index(['year', 'plot'])['buds breaking'].unstack()
pheno_obs.columns = ['%02d' % c for c in pheno_obs.columns]


annt2m = pd.DataFrame()
tsoi_collect   = {}
for fid in chamber_levels.keys():
    env       = pd.read_csv(os.path.join(path_input, 'WEW_Complete_Environ_20200420',
                                        'WEW PLOT_{}_Complete_Environ_20200417.csv'.format(fid)))
    env       = env.loc[(env['Year'] >= 2015) & (env['Year'] <= 2019), :]
    env.index = pd.DatetimeIndex(env['TIMESTAMP'])
    env       = env.replace(to_replace = '^\s+', value = np.nan, regex = True).sort_index()
    env       = env.loc[~env.index.duplicated(keep = 'first'), :]

    tseries   = (env['TA_2_0__1'].astype(float) + env['TA_2_0__2'].astype(float))/2

    #tmax      = tseries.groupby(env.index.year * 1000 + env.index.month * 100 + env.index.day).max()
    #tmin      = tseries.groupby(env.index.year * 1000 + env.index.month * 100 + env.index.day).min()
    #tavg      = tseries.groupby(env.index.year * 1000 + env.index.month * 100 + env.index.day).mean()
    annt2m[fid]     = tseries.groupby(env.index.year).mean()

    soil_layer_list = [200, 0, -5, -10, -20, -30, -40, -50, -100, -200]
    tsoi            = env.loc[:, 
        ['TS_%d__A%d' % (i,ind+1) for ind, i in enumerate(soil_layer_list[1:])]]
    tsoi.columns    = ['TS_%d' % i for i in soil_layer_list[1:]]

    # check: tseries.loc[(tsoi.index[0] - timedelta(days = 119, hours = 23, minutes = 30)):(tsoi.index[0])].mean()
    tsoi.loc[:, 'TS_200'] = tseries.rolling('21d').mean()
    tsoi.loc[tsoi.index < tseries.index[21*48], 'TS_200'] = np.nan
    tsoi = dict(tsoi.loc[:, ['TS_%d' % i for i in soil_layer_list]].dropna(axis = 0, how = 'any'))

    for k in tsoi.keys():
        if k not in tsoi_collect:
            tsoi_collect[k] = pd.DataFrame({fid: tsoi[k]})
        else:
            tsoi_collect[k][fid] = tsoi[k]


def crit_gdd(x, a, b):
    # crit_gdd  = pd.DataFrame(np.exp(4.8 + 0.13 * annt2m.values), index = annt2m.index + 1)
    return np.exp(a + b * x)
gdds = dict([('GDD_%d' % i,
              pd.DataFrame(index = [2016, 2017, 2018, 2019], columns = chamber_levels.keys()))
             for i in soil_layer_list])

for key, year in it.product(tsoi_collect.keys(), [2016, 2017, 2018, 2019]):
    # Find last solstice
    last_year = pd.DatetimeIndex([datetime(year-1, 12, 1) + timedelta(days = i) \
                                 for i in range((datetime(year,1,1)-datetime(year-1, 1, 1)).days)])

    dl_prev = 1e20
    for ind, yy, mm, dd in zip(range(len(last_year)), 
                                     last_year.year, last_year.month, last_year.day):
        dl = daylength(yy, mm, dd, lat, lon, elv, tzn)
        if dl >= dl_prev:
            break
        else:
            dl_prev = dl
    solstice = last_year[ind]

    # Accumulate GDD
    for fid in chamber_levels.keys():
        tsoi_accu = tsoi_collect[key].loc[(tsoi_collect[key].index >= solstice) & \
                                        (tsoi_collect[key].index <= (datetime(year-1, 12, 31) \
            + timedelta(days = float(pheno_obs.loc[year, fid])))), fid]
        tsoi_accu = (tsoi_accu * (tsoi_accu > 0.)).sum(axis = 0) * \
                    (tsoi_collect[key].index[1] - tsoi_collect[key].index[0]).total_seconds() / 86400
        gdds['GDD_' + key.split('_')[1]].loc[year, fid] = tsoi_accu


Optim = pd.DataFrame(index = ['a', 'b', 'RMSE'], columns = gdds.keys())

mpl.rcParams['font.size'] = 6
mpl.rcParams['axes.titlesize'] = 6
fig, axes = plt.subplots(3, 4, figsize = (9,9), sharex = True, sharey = True)
for i, lyr in enumerate(soil_layer_list):
    x = annt2m.loc[gdds['GDD_%d' % lyr].index - 1, :].values.reshape(-1)
    y = gdds['GDD_%d' % lyr].values.reshape(-1)
    popt, pcov = curve_fit(crit_gdd, x, y)
    ypred = crit_gdd(x, *popt)

    Optim.loc['a', 'GDD_%d' % lyr] = popt[0]
    Optim.loc['b', 'GDD_%d' % lyr] = popt[1]
    Optim.loc['RMSE', 'GDD_%d' % lyr] = np.sqrt(np.mean(np.power(y - ypred, 2)))

    ax = axes.flat[i]
    h1, = ax.plot(x, y, 'or', markersize = 3)
    h2, = ax.plot(np.linspace(6, 15, 40), crit_gdd(np.linspace(5, 15, 40), *popt), '-b')

    ax.text(0.6, 0.85, 'RMSE=%.2f' % Optim.loc['RMSE', 'GDD_%d' % lyr], transform = ax.transAxes)

    ax.set_title('Soil layer %s cm' % lyr)
    ax.set_xlabel('T2M previous year average')
    ax.set_ylabel('GDD at budburst')
ax.legend([h1,h2], ['Ground observation', 'Critical GDD model'], loc = (-1, -0.4))
fig.savefig(os.path.join(path_out, 'fit_ground_phenology.png'), dpi = 600.)
plt.close(fig)

Optim.to_csv(os.path.join(path_out, 'fit_ground_phenology.csv'))


##################################################################################################
# Critical daylength model for leaves senescing
##################################################################################################
pheno_obs = pd.read_excel(os.path.join(path_input, 'SPRUCE_budburst_summary.xlsx'),
                          engine = 'openpyxl', sheet_name = 'Sheet1')
pheno_obs = pheno_obs.loc[:, ['year', 'leaves senescing']].set_index('year', drop = True).dropna()
crit_dayl_all = []
for i in range(len(pheno_obs)):
    fall_date = datetime(pheno_obs.index[i]-1, 12, 31) + timedelta(float(pheno_obs.iloc[i]))
    crit_dayl = daylength(fall_date.year, fall_date.month, fall_date.day, lat, lon, elv, tzn)
    print('Critical day length = %f seconds' % crit_dayl)
    crit_dayl_all.append(crit_dayl)
print('Average critical day length = %f seconds' % np.mean(crit_dayl_all))