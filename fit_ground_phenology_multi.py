from pyPhenology import utils
import numpy as np
import pandas as pd
import os
from utils.constants import *
from utils.paths import *
from utils.analysis import *
from utils.phenofuncs import *
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib as mpl


pheno_obs = read_leaf_sos()['EN']
tsoi_collect, _ = read_obs_tsoi_daily()
soil_layer_list = ['2m', '10cm']

observations = pheno_obs.stack().to_frame('doy')
observations.index.names = ['year', 'site_id']
observations = observations.reset_index()

predictors = {}
for lyr in soil_layer_list:
    temp = tsoi_collect[lyr].loc[:, pheno_obs.columns]
    temp.loc['2016-02-29', :] = 0.5 * (temp.loc['2016-02-28', :].values + temp.loc['2016-03-01', :].values)
    temp.loc['2020-02-29', :] = 0.5 * (temp.loc['2020-02-28', :].values + temp.loc['2020-03-01', :].values)
    temp = temp.sort_index(axis = 0)

    temp = temp.stack().to_frame('temperature')
    temp.index.names = ['time', 'site_id']
    temp = temp.reset_index()

    temp['year'] = pd.DatetimeIndex(temp['time']).year
    temp['doy'] = pd.DatetimeIndex(temp['time']).dayofyear
    temp['daylength'] = [daylength_simple(doy, lat = 47.563) for  doy in temp['doy'].values]

    predictors[lyr] = temp.drop('time', axis = 1)

models_to_test = ['ThermalTime', 'M1', 'Uniforc', 'Sequential', 'Unichill', 'Alternating']

np.random.seed(100)
test = np.random.rand(observations.shape[0]) >= 0.85
observations_test = observations.loc[test, :]
observations_train = observations.loc[~test, :]

# AIC based off mean sum of squares
def aic(obs, pred, n_param):
    return len(obs) * np.log(np.mean((obs - pred)**2)) + 2*(n_param + 1)

def rmse(obs, pred):
    return np.sqrt(np.mean((obs - pred)**2))

for lyr in ['2m', '10cm']:
    print(lyr)

    best_aic=np.inf
    best_base_model = None
    best_base_model_name = None

    for model_name in models_to_test:
        Model = utils.load_model(model_name)
        model = Model()
        model.fit(observations_train, predictors[lyr], optimizer_params='practical')

        obs = observations_test.doy.values
        pred = model.predict(observations_test, predictors[lyr])

        model_aic = aic(obs = obs, pred = pred, n_param = len(model.get_params()))
        model_rmse = rmse(obs = obs, pred = pred)
        
        if model_aic < best_aic:
            best_model = model
            best_model_name = model_name
            best_aic = model_aic

        print('model {m} got an aic of {a} rmse of {b}'.format(m = model_name, a = model_aic, b = model_rmse))

    print('Best model: {m}'.format(m=best_model_name))
    print('Best model paramters:')
    print(best_model.get_params())

# 2m
# model ThermalTime got an aic of 40.1959427138737 rmse of 5.981452814975453
# model M1 got an aic of 40.97017082797774 rmse of 5.587684871413403
# model Uniforc got an aic of 41.09738142639188 rmse of 5.627314338711377
# model Sequential got an aic of 51.732627221969096 rmse of 7.280109889280518
# model Unichill got an aic of 49.034000882763166 rmse of 5.6075346137535735
# model Alternating got an aic of 43.5289709644902 rmse of 5.7638721552635275

# 10cm
# model ThermalTime got an aic of 52.531744018948565 rmse of 11.86966254317657
# model Sequential got an aic of 63.55898012115967 rmse of 14.043582955293932
# model Unichill got an aic of 60.18535867235907 rmse of 10.418999738725189
# model Alternating got an aic of 52.14707821681076 rmse of 9.303523824635242


"""
################################################################################
# Check the accumulation of GDD v.s. date of the year, until the day of leaf-out
################################################################################
lyr = '2m'

Model = utils.load_model('Alternating')
model = Model()
model.fit(observations_train, predictors[lyr], optimizer_params='practical')

a = model.get_params()['a']
b = model.get_params()['b']
c = model.get_params()['c']
threshold = model.get_params()['threshold']
t1 = model.get_params()['t1']

chill_days = ((tsoi_collect[lyr] < threshold) * 1).copy()
chill_days.loc[tsoi_collect[lyr].index.dayofyear < t1, :] = 0
chill_days = chill_days.loc[chill_days.index.year > 2015, :]
for year in range(2016, 2022):
    chill_days.loc[chill_days.index.year == year, :] = chill_days.loc[chill_days.index.year == year, :].cumsum(axis = 0)

# Accumulated growing degree days from Jan 1
gdd = tsoi_collect[lyr].copy()
gdd[gdd < threshold] = 0
gdd.loc[tsoi_collect[lyr].index.dayofyear < t1, :]= 0
gdd = gdd.loc[gdd.index.year > 2015, :]
for year in range(2016, 2022):
    gdd.loc[gdd.index.year == year, :] = gdd.loc[gdd.index.year == year, :].cumsum(axis = 0)

# Phenological event happens the first day gdd is > chill_day curve
chill_day_curve = a + b * np.exp(c * chill_days)
difference = gdd - chill_day_curve

# The estimate is equal to the first day that
# gdd - chill_day_curve > 0
before_leafout = difference <= 0

# plot until leaf out
# GDD hasn't even started accumulating when root starts growing
fig, axes = plt.subplots(4, 3, figsize = (12, 12), sharex = True, sharey = True)
for i, plot in enumerate(pheno_obs.columns):
    ax = axes.flat[i]
    temp = chill_days[plot].loc[before_leafout[plot]]
    temp = temp.loc[temp.index.year == 2018]
    ax.plot(temp.index, temp.values, 'ob', ms = 2)
    temp2 = gdd[plot].loc[before_leafout[plot]]
    temp2 = temp2.loc[temp2.index.year == 2018]
    ax.plot(temp2.index, temp2.values, 'or', ms = 2)
    plt.setp(ax.get_xticklabels(), rotation = 90)
for i in range(len(pheno_obs.columns), 12):
    ax = axes.flat[i]
    plt.setp(ax.get_xticklabels(), rotation = 90)
fig.savefig(os.path.join(path_out, 'fit_ground_phenology_multi.png'), dpi = 600., bbox_inches = 'tight')
plt.close(fig)

# using 5 degree threshold, plot the accumulated chiling days
# use 68 chilling threshold for shrub, 
# use 80 chilling threshold, or consecutive 3 days > 5 degree (i.e. unable to accumulate chilling anymore), for trees
fig, ax = plt.subplots(figsize = (10, 10))
for i, plot in enumerate(pheno_obs.columns):
    temp = chill_days[plot].loc[before_leafout[plot]]
    temp = temp.loc[temp.index.year == 2018]
    ax.plot(temp.index, temp.values, 'o', ms = 2)
ax.axvline(datetime(2018, 3, 10))
ax.axvline(datetime(2018, 4, 4))
ax.set_yticks(range(0, 121, 5))
fig.savefig(os.path.join(path_out, 'fit_ground_phenology_multi_chil.png'), dpi = 600., bbox_inches = 'tight')
plt.close(fig)

# plot simply the soil temperature series
# shrub: 5 day, 1 degree, using t21; tsoil doesn't look as good
# tree: 20 days, 4 degree, using 21; tsoil doesn't look as good
fig, ax = plt.subplots(figsize = (10, 10))
for i, plot in enumerate(pheno_obs.columns):
    temp = tsoi_collect[lyr][plot].loc[tsoi_collect[lyr].index.year > 2015].loc[before_leafout[plot]]
    temp = temp.loc[temp.index.year == 2018]
    ax.plot(temp.index, temp.values, 'o', ms = 2)
ax.grid(axis = 'y')
ax.axvline(datetime(2018, 3, 10))
ax.axvline(datetime(2018, 4, 4))
ax.set_yticks(range(-5, 15, 1))
fig.savefig(os.path.join(path_out, 'fit_ground_phenology_multi_t.png'), dpi = 600., bbox_inches = 'tight')
plt.close(fig)
"""