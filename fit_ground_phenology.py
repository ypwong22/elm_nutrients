import pandas as pd
import os
from pandas.core.indexes.datetimes import DatetimeIndex
from utils.constants import *
from utils.paths import *
from utils.analysis import *
from utils.phenofuncs import *
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib as mpl
import itertools as it
import xarray as xr


lat, lon, tzn, elv = '47.563 N', '93.438 W', 'US/Central', 430


##################################################################################################
# Load the start of season data and soil temperature
##################################################################################################
hr = xr.open_dataset(os.path.join(path_intrim, 'spruce_validation_data.nc'))
pheno_obs = hr['pheno_dates_lai'].loc[:, 'SOS', :, 'EN']
pheno_obs = pd.DataFrame(pheno_obs.values,
                         index = pheno_obs['year'],
                         columns = ['%02d' % i for i in pheno_obs['chamber']])
pheno_obs = pheno_obs.dropna(axis = 1, how = 'all' \
    ).dropna(axis = 0, how = 'all').sort_index(axis = 1)
hr.close()


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


##################################################################################################
# Growing degree days model for budburst
##################################################################################################
# match the annt2m's year to the next year
annt2m = annt2m.loc[(annt2m.index >= (pheno_obs.index[0]-1)) & (annt2m.index <= (pheno_obs.index[-1]-1)), :]
annt2m = annt2m.loc[:, pheno_obs.columns]

# Test the fitting outcome of different soil layers
Optim = pd.DataFrame(np.nan, columns = ['a', 'b', 'RMSE'], index = soil_layer_list)

mpl.rcParams['font.size'] = 6
mpl.rcParams['axes.titlesize'] = 6
fig, axes = plt.subplots(3, 4, figsize = (9,9), sharex = True, sharey = True)
for i,lyr in enumerate(soil_layer_list):
    model = GDD()

    popt, ypred, x, y = model.fit(pheno_obs, annt2m, tsoi_collect[f'TS_{lyr:d}'])

    Optim.loc['GDD_%d' % lyr, 'a'   ] = popt[0]
    Optim.loc['GDD_%d' % lyr, 'b'   ] = popt[1]
    Optim.loc['GDD_%d' % lyr, 'RMSE'] = np.sqrt(np.mean(np.power(y - ypred, 2)))

    ax = axes.flat[i]
    h1, = ax.plot(x, y, 'or', markersize = 3)
    h2, = ax.plot(np.linspace(6, 15, 40), model.predict(np.linspace(5, 15, 40)), '-b')

    ax.text(0.6, 0.85, 'RMSE=%.2f' % Optim.loc['GDD_%d' % lyr, 'RMSE'], transform = ax.transAxes)

    ax.set_title('Soil layer %s cm' % lyr)
    ax.set_xlabel('T2M previous year average')
    ax.set_ylabel('GDD at budburst')
ax.legend([h1,h2], ['Ground observation', 'Critical GDD model'], loc = (-1, -0.4))
fig.savefig(os.path.join(path_out, 'fit_ground_phenology.png'), dpi = 600.)
plt.close(fig)
Optim.to_csv(os.path.join(path_out, 'fit_ground_phenology.csv'))


##################################################################################################
# Alternating model for budburst
##################################################################################################




##################################################################################################
# Critical daylength model for leaves senescing
##################################################################################################
pheno_obs = pd.read_excel(os.path.join(path_input, 'SPRUCE_budburst_summary.xlsx'),
                          engine = 'openpyxl', sheet_name = 'EN')
pheno_obs = pheno_obs.loc[:, ['year', 'leaves senescing']].set_index('year', drop = True).dropna()
crit_dayl_all = []
for i in range(len(pheno_obs)):
    fall_date = datetime(pheno_obs.index[i]-1, 12, 31) + timedelta(float(pheno_obs.iloc[i]))
    crit_dayl = daylength(fall_date.year, fall_date.month, fall_date.day, lat, lon, elv, tzn)
    print('Critical day length = %f seconds' % crit_dayl)
    crit_dayl_all.append(crit_dayl)
print('Average critical day length = %f seconds' % np.mean(crit_dayl_all))
print('Minimum critical day length = %f seconds' % np.min(crit_dayl_all))
print('Minimum critical day length = %f seconds' % np.max(crit_dayl_all))
print('Std of day length = %f seconds' % np.std(crit_dayl_all))
