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

models_to_test = ['ThermalTime', 'M1', 'Alternating'] # 'Uniforc', 'Sequential', 'Unichill', 'Alternating']

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
    # model ThermalTime got an aic of 38.677197400638 rmse of 5.497474167490214
    # model M1 got an aic of 42.1959427138737 rmse of 5.981452814975453
    # model Uniforc got an aic of 41.09738142639188 rmse of 5.627314338711377
    # model Sequential got an aic of 51.732627221969096 rmse of 7.280109889280518
    # model Unichill got an aic of 51.44650585471845 rmse of 6.411794687223781
    # model Alternating got an aic of 45.59128109448307 rmse of 6.4635731432217725

    # 10cm
    # model ThermalTime got an aic of 49.545980199572085 rmse of 10.055402085998905
    # model M1 got an aic of 54.531744018948565 rmse of 11.86966254317657
    # model Uniforc got an aic of 51.70284691578964 rmse of 10.143416036468626
    # model Sequential got an aic of 62.764346656848566 rmse of 13.437096247164249
    # model Unichill got an aic of 60.52877382725702 rmse of 10.619688214715994
    # model Alternating got an aic of 52.14707821681076 rmse of 9.303523824635242


lyr_select = '2m'

# Note t1 - the start accumulation date, is counted from 12.21

Model = utils.load_model('ThermalTime')
model = Model()
model.fit(observations_train, predictors[lyr_select], optimizer_params='practical')
print(model.get_params())
# {'t1': 88.96680755982837, 'T': -6.210506141223593, 'F': 995.5839297751559}

Model = utils.load_model('Alternating')
model = Model()
model.fit(observations_train, predictors[lyr_select], optimizer_params='practical')
print(model.get_params())
# {'a': 199.30927021132482, 'b': 3751.8473351447456, 'c': -0.03728095348592797, 'threshold': 5, 't1': 1}

Model = utils.load_model('Uniforc')
model = Model()
model.fit(observations_train, predictors[lyr_select], optimizer_params='practical')
print(model.get_params())
# {'t1': 83.27265975501095, 'F': 37.17482725282362, 'b': -0.13936174741659002, 'c': 7.587839240505523}