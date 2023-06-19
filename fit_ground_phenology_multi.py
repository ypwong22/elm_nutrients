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


# Because ELM starts accumulation on winter solstice, 
# roll the dates 12.21 -> 1.1: add 11 days
pheno_obs = pheno_obs + 11
for lyr in soil_layer_list:
    temp = tsoi_collect[lyr]
    temp.loc['2016-02-29', :] = 0.5 * (temp.loc['2016-02-28', :].values + temp.loc['2016-03-01', :].values)
    temp.loc['2020-02-29', :] = 0.5 * (temp.loc['2020-02-28', :].values + temp.loc['2020-03-01', :].values)
    tsoi_collect[lyr] = pd.DataFrame(temp.iloc[:-11, :].values, 
                                     index = temp.index[11:],
                                     columns = temp.columns)


observations = pheno_obs.stack().to_frame('doy')
observations.index.names = ['year', 'site_id']
observations = observations.reset_index()

predictors = {}
for lyr in soil_layer_list:
    temp = tsoi_collect[lyr].loc[:, pheno_obs.columns]
    temp = temp.sort_index(axis = 0)

    temp = temp.stack().to_frame('temperature')
    temp.index.names = ['time', 'site_id']
    temp = temp.reset_index()

    temp['year'] = pd.DatetimeIndex(temp['time']).year
    temp['doy'] = pd.DatetimeIndex(temp['time']).dayofyear

    # shift back 11 days to get the correct daylength
    doy_list_correct = (pd.DatetimeIndex(temp['time']) - timedelta(days = 11)).dayofyear
    temp['daylength'] = [daylength_simple(doy, lat = 47.563) for doy in doy_list_correct]

    predictors[lyr] = temp.drop('time', axis = 1)

#models_to_test = ['ThermalTime', 'M1', 'Alternating']
models_to_test = ['Uniforc', 'Sequential', 'Unichill']

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
            if model_name in ['ThermalTime', 'M1', 'Uniforc', 'Alternating']:
                model = Model(parameters = {'t1': 1})
            elif model_name in ['Sequential', 'Unichill']:
                model = Model(parameters = {'t0': 1})
            else:
                raise Exception('Not implemented')
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
    # model ThermalTime got an aic of 44.22816521661635 rmse of 7.483314773547883
    # model M1 got an aic of 46.22816521661635 rmse of 7.483314773547883
    # model Uniforc got an aic of 42.47119352422242 rmse of 6.073622386096199
    # model Sequential got an aic of 50.69378132429378 rmse of 6.871842709362768
    # model Unichill got an aic of 48.54386162857173 rmse of 5.456901847914967
    # model Alternating got an aic of 46.03495068808412 rmse of 6.624868971450597

    # 10cm
    # model ThermalTime got an aic of 60.632388253444255 rmse of 18.61600267392427
    # model M1 got an aic of 62.632388253444255 rmse of 18.61600267392427
    # model Uniforc got an aic of 52.59062974084688 rmse of 10.656244908763854
    # model Sequential got an aic of 62.764346656848566 rmse of 13.437096247164249
    # model Unichill got an aic of 59.780266748033604 rmse of 10.187137859957417
    # model Alternating got an aic of 64.72425507978949 rmse of 18.711256267581582


lyr_select = '2m'

# Note t1 - the start accumulation date, is counted from 12.21

Model = utils.load_model('ThermalTime')
model = Model(parameters = {'t1': 1})
model.fit(observations_train, predictors[lyr_select], optimizer_params='practical')
print(model.get_params())
# {'T': 4.08239397438806, 'F': 388.89698143143914, 't1': 1}
# {'T': 4.113410776245213, 'F': 386.9374360874435, 't1': 1}

Model = utils.load_model('Alternating')
model = Model(parameters = {'t1': 1})
model.fit(observations_train, predictors[lyr_select], optimizer_params='practical')
print(model.get_params())
# {'a': 349.6277253204969, 'b': 3169.8267073240854, 'c': -0.05347986748356126, 'threshold': 5, 't1': 1}
# {'a': 195.28932107981944, 'b': 3138.652468114228, 'c': -0.03495823449748503, 'threshold': 5, 't1': 1}

Model = utils.load_model('Uniforc')
model = Model(parameters = {'t1': 1})
model.fit(observations_train, predictors[lyr_select], optimizer_params='practical')
print(model.get_params())
# {'F': 66.94329645213861, 'b': -0.04412234589076469, 'c': 6.848386620467728, 't1': 1}
# {'F': 53.305745769355205, 'b': -0.03781848702688784, 'c': 17.92139941787829, 't1': 1}