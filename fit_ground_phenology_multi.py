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


find_best = False
if find_best:

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
    # model Sequential got an aic of 51.3274140896035 rmse of 7.118052168020874
    # model Unichill got an aic of 49.19162312519754 rmse of 5.656854249492381
    # model Alternating got an aic of 43.25390711079871 rmse of 5.676462121975467

    # 10cm
    # model ThermalTime got an aic of 52.531744018948565 rmse of 11.86966254317657
    # model M1 got an aic of 54.531744018948565 rmse of 11.86966254317657
    # model Uniforc got an aic of 51.74164007929974 rmse of 10.16530045465127
    # model Sequential got an aic of 63.55898012115967 rmse of 14.043582955293932
    # model Unichill got an aic of 60.52877382725702 rmse of 10.619688214715994
    # model Alternating got an aic of 52.112351485622995 rmse of 9.285592184789413



lyr_select = '2m'


Model = utils.load_model('ThermalTime')
model = Model()
model.fit(observations_train, predictors[lyr_select], optimizer_params='practical')
print(model.get_params())
# {'t1': 82.89828273347462, 'T': 0.001321356324646139, 'F': 590.7274692372864}

Model = utils.load_model('Alternating')
model = Model()
model.fit(observations_train, predictors[lyr_select], optimizer_params='practical')
print(model.get_params())
# {'a': 281.3048492792021, 'b': 2249.8346532974097, 'c': -0.02412355913509323, 'threshold': 5, 't1': 1}