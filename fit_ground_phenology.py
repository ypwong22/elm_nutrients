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
from scipy.stats import linregress


plot_list =  ['04', '06','07', '08', '10', '11', '13', '16', '17', '19', '20', '21']

##################################################################################################
# Load the start of season data and soil temperature
##################################################################################################
pheno_obs = read_leaf_sos()['EN']
annt2m = read_sims_tair_annual() # observed values have missing values in 2015
tsoi_collect, _ = read_obs_tsoi_daily()
soil_layer_list = ['2m', '10cm']

# match columns
annt2m.columns = [f'{t:02g}' for t in annt2m.columns]
annt2m = annt2m.loc[:, pheno_obs.columns]
for lyr in soil_layer_list:
    tsoi_collect[lyr] = tsoi_collect[lyr].loc[:, pheno_obs.columns]
    tsoi_collect[lyr] = tsoi_collect[lyr].loc[tsoi_collect[lyr].index.year >= (annt2m.index[0] + 1), :]
    tsoi_collect[lyr] = tsoi_collect[lyr].loc[(tsoi_collect[lyr].index.year >= pheno_obs.index[ 0]) & \
                                              (tsoi_collect[lyr].index.year <= pheno_obs.index[-1]), :]

##################################################################################################
# Growing degree days model for budburst
##################################################################################################
# Test the fitting outcome of different soil layers
Optim = pd.DataFrame(np.nan, columns = ['a', 'b', 'RMSE'], index = soil_layer_list)

mpl.rcParams['font.size'] = 6
mpl.rcParams['axes.titlesize'] = 6
fig, axes = plt.subplots(3, 4, figsize = (9,9), sharex = True, sharey = True)
for i,lyr in enumerate(soil_layer_list):
    model = GDD(pheno_obs, tsoi_collect[lyr], annt2m)
    popt, rmse, ypred = model.fit()

    Optim.loc[lyr, 'a'   ] = popt[0]
    Optim.loc[lyr, 'b'   ] = popt[1]
    Optim.loc[lyr, 'RMSE'] = rmse

    ax = axes.flat[i]
    ax.plot(pheno_obs.values.reshape(-1), ypred.values.reshape(-1), 'or', markersize = 3)

    res = linregress(ypred.values.reshape(-1), pheno_obs.values.reshape(-1))
    xmin = pheno_obs.values.reshape(-1).min()
    xmax = pheno_obs.values.reshape(-1).max()

    ax.plot([xmin, xmax], [xmin * res.slope + res.intercept, xmax * res.slope + res.intercept], '-b')

    ax.text(0.6, 0.85, 'RMSE=%.2f' % Optim.loc[lyr, 'RMSE'], transform = ax.transAxes)

    ax.set_title('Soil layer %s cm' % lyr)
    ax.set_xlabel('Predicted budbreak')
    ax.set_ylabel('Observed budbreak')
fig.savefig(os.path.join(path_out, 'fit_ground_phenology.png'), dpi = 600.)
plt.close(fig)
Optim.to_csv(os.path.join(path_out, 'fit_ground_phenology.csv'))
annt2m.to_csv(os.path.join(path_out, 'fit_ground_phenology_annt2m.csv'))


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
