""" Calculate the Bias, uRMSE, and Corr against Paul's data 
    Plot the results
    Also plot the FROOTC_TO_LITTER time series against Colleen's data
"""
from utils.paths import *
from utils.constants import *
from utils.analysis import *
import numpy as np
import xarray as xr
import pandas as pd
import os
from scipy.stats import pearsonr
import seaborn as sns
import matplotlib.pyplot as plt
from glob import glob


expr_list = ['20221212', '20221231', '20230101']
var_list = ['ANPPtree', 'ANPPshrub', 'NPPmoss', 'BNPP', 'HR', 'NEE']


collection = pd.DataFrame(np.nan, index = pd.MultiIndex.from_product([['Bias (%)', 'uRMSE (%)', 'Corr'], chamber_list_names]),
                          columns = pd.MultiIndex.from_product([expr_list, var_list]))


hr = xr.open_dataset(os.path.join(path_intrim, 'spruce_validation_data.nc'))

for expr in expr_list:

    for plot, plot_name in zip(chamber_list, chamber_list_names):

        data = pd.read_csv(os.path.join(path_out, expr, f'plot{plot}_ts_extract.csv'), parse_dates = True, index_col = 0, header = [0, 1])

        for var in var_list:
            if var == 'ANPPtree': 
                obs = pd.Series(hr['annual_anpp_tree'].loc[:, plot].values, index = hr['year']).dropna()

                # hummock: 0.64, hollow: 0.36
                # pima: 0.36, lala: 0.14
                sim = ((data.loc[:, ('hummock', 'AGNPP_lala')] * 0.14 + \
                        data.loc[:, ('hummock', 'AGNPP_pima')] * 0.36) * 0.64 + \
                    (data.loc[:, ('hollow', 'AGNPP_lala')] * 0.14 + \
                        data.loc[:, ('hollow', 'AGNPP_pima')] * 0.36) * 0.36).groupby(data.index.year).sum()
            elif var == 'ANPPshrub':
                obs = pd.Series(hr['annual_anpp_shrub'].loc[:, plot].values, index = hr['year']).dropna()

                # hummock: 0.64, hollow: 0.36
                sim = (data.loc[:, ('hummock', 'AGNPP_shrub')] * 0.64 + \
                       data.loc[:, ('hollow', 'AGNPP_shrub')] * 0.36).groupby(data.index.year).sum() * 0.25
            elif var == 'NPPmoss':
                obs = pd.Series(hr['annual_npp_moss'].loc[:, plot].values, index = hr['year']).dropna()

                # hummock: 0.64, hollow: 0.36
                sim = (data.loc[:, ('hummock', 'NPP_moss')] * 0.64 + \
                       data.loc[:, ('hollow', 'NPP_moss')] * 0.36).groupby(data.index.year).sum()

                mossfrac = pd.read_excel('Sphagnum_fraction.xlsx', index_col = 0, skiprows = 1,
                                        engine = 'openpyxl').drop(['plot','Temp','CO2'], axis = 1)
                mossfrac[2015] = mossfrac[2016]
                mossfrac = mossfrac.drop(2021, axis = 1)

                sim = sim * mossfrac.loc[plot_name] / 100.

            elif var == 'BGNPP':
                obs = pd.Series(hr['annual_bnpp'].loc[:, plot].values, index = hr['year']).dropna()
                
                # hummock: 0.64, hollow: 0.36
                sim = ((data.loc[:, ('hummock', 'BGNPP_pima')] * 0.36 + \
                        data.loc[:, ('hummock', 'BGNPP_lala')] * 0.14 + \
                        data.loc[:, ('hummock', 'BGNPP_shrub')] * 0.25) * 0.64 + \
                       (data.loc[:, ('hollow', 'BGNPP_pima')] * 0.36 + \
                        data.loc[:, ('hollow', 'BGNPP_lala')] * 0.14 + \
                        data.loc[:, ('hollow', 'BGNPP_shrub')] * 0.25) * 0.36).groupby(data.index.year).sum()
            
            elif  var == 'HR':
                obs = pd.Series(hr['annual_rh'].loc[:, plot].values, index = hr['year']).dropna()

                sim = - (data.loc[:, ('hummock', 'HR')] * 0.64 + data.loc[:, ('hollow', 'HR')] * 0.36).groupby(data.index.year).sum()

            elif  var == 'NEE':
                obs = pd.Series(hr['annual_nee'].loc[:, plot].values, index = hr['year']).dropna()

                sim = - (data.loc[:, ('hummock', 'NEE')] * 0.64 + data.loc[:, ('hollow', 'NEE')] * 0.36).groupby(data.index.year).sum()

            ind = sim.index.intersection(obs.index)
            sim = sim.loc[ind]
            obs = obs.loc[ind]

            collection.loc[('Bias (%)', plot_name), (expr, var)] = 100 * (sim.mean() - obs.mean()) / np.abs(obs.mean())
            collection.loc[('uRMSE (%)', plot_name), (expr, var)] = 100 * np.sqrt(np.sum(np.power(sim - sim.mean() - obs + obs.mean(), 2))) / np.abs(obs.mean())
            collection.loc[('Corr', plot_name), (expr, var)], _ = pearsonr(sim, obs)

collection.to_csv(os.path.join(path_out, 'validation_metrics.csv'))


# Plot against Paul Hanson's data
fig, axes = plt.subplots(3, len(var_list), sharey = False, sharex = True, figsize = (17, 9))
for i, metric in enumerate(['Bias (%)', 'uRMSE (%)', 'Corr']):
    for j, var in enumerate(var_list):
        ax = axes[i, j]

        temp = collection.loc[metric, (slice(None), var)]
        temp.columns = temp.columns.droplevel(1)

        sns.barplot(data = temp, ax = ax, palette = 'husl')

        if i == 0:
            ax.set_title(var)
        if j == 0:
            ax.set_ylabel(metric)
        else:
            if (i == 2) or (j != (len(var_list)-1)):
                ax.set_yticklabels([])

        if i == 0:
            if j != (len(var_list) - 1):
                ax.set_ylim([-50, 150])
            else:
                ax.set_ylim([0, 350])
        elif i == 1:
            if j != (len(var_list) - 1):
                ax.set_ylim([0, 250])
            else:
                ax.set_ylim([0, 800])
        else:
            ax.set_ylim([-0.75, 0.8])

        if i != 2:
            ax.set_xticklabels([])
        else:
            ax.set_xticklabels(['Default', 'Evergreen', 'EvgrRoot'], rotation = 90)

fig.savefig(os.path.join(path_out, 'validation_metrics.png'), dpi = 600., bbox_inches = 'tight')
plt.close(fig)




# Plot against Colleen's data
# Seasonality
mort, _ = read_mortality()        
pft_list = [2, 3, 11]
pft_names = ['Pima', 'Lala', 'Shrub']
sims_prefix = ['20221212', '20221231', '20230101']
sims_names  = ['Default', 'Evergreen', 'EvgrRoot']
clist = ['#fcbba1', '#c6dbef', '#fc9272', '#9ecae1', '#fb6a4a', '#6baed6', '#de2d26', '#3182bd', '#a50f15', '#08519c']


varname = 'FROOTC_TO_LITTER'

fig, axes = plt.subplots(len(pft_list), len(sims_prefix), figsize = (12, 8), sharex = True, sharey = False)
fig.subplots_adjust(wspace = 0.06, hspace = 0.03)
for i, pft in enumerate(pft_list):
    for j, prefix in enumerate(sims_prefix):

        ax = axes[i, j]

        h = [None] * len(chamber_list)
        for k, plot in enumerate(chamber_list):
            path_data = os.path.join(os.environ['PROJDIR'], 'E3SM', 'output', f'{prefix}_plot{plot:02d}_US-SPR_ICB20TRCNPRDCTCBC', 'run')
            flist = sorted(glob(os.path.join(path_data, '*.h2.*.nc')))

            # osbolete # autocoversion misidentified PROOT_ONSET_OFFSET as a datetime variable because its unit "days"
            hr = xr.open_mfdataset(flist, decode_times = False)
            units, reference_date = hr.time.attrs['units'].split('since')
            tvec = pd.date_range(start = reference_date, end = '2021-01-01', freq='D')
            tvec = tvec[(tvec.month != 2) | (tvec.day != 29)]
            hr['time'] = tvec
            tvec = [datetime.strptime(x, '%Y-%m-%d %H') for x in hr['time'].indexes['time'].strftime('%Y-%m-%d %H')]
            retype = pd.DataFrame(hr[varname][:, pft].values * 0.64 + hr[varname][:,  pft + 17].values * 0.36, index = tvec)

            hr.close()

            if varname in ['GPP', 'QVEGE', 'QVEGT', 'LEAFC_TO_LITTER', 'FROOTC_TO_LITTER', 'LEAFC_STORAGE_TO_XFER', 'LEAFC_XFER_TO_LEAFC']:
                # s to year
                retype = retype * 24 * 3600 * 365
            if varname == 'GPP':
                retype2 = pd.DataFrame(savgol_filter(retype.values, 81, 3, mode = 'constant', cval = np.nan), index = tvec)
            if varname == 'QVEGE':
                retype = retype.iloc[1:]

            dayofyear = []
            for yy in np.unique(retype.index.year[:-1]):
                dayofyear = dayofyear + list(np.arange(1, 366))
            dayofyear = dayofyear + [1] # restart file
            dayofyear = np.array(dayofyear)
            retype = retype.groupby(dayofyear).mean()
            h[k], = ax.plot(retype.index, retype.values, color = clist[k])

        if j == 0:
            ax.set_ylabel(f'{varname}\ngC m-2 year-1\n{pft_names[i]}')
        else:
            ax.set_yticklabels([])
        if i == 0:
            ax.set_title(sims_names[j])
        elif i == 2:
            ax.set_xlabel('Day of year')

        if i == 0:
            ax.set_ylim([0, 600])
        elif i == 1:
            ax.set_ylim([0, 2000])
        else:
            ax.set_ylim([0, 3000])

        ax2 = ax.twinx()
        ax2.plot(mort.index.dayofyear, mort, '-k', lw = 2)
        if j == 2:
            ax2.set_ylabel('Mortality mm day-1')
        if j != 2:
            ax2.set_yticklabels([])

        ax.set_xlim([0, 365])

ax.legend(h, chamber_list_names, loc = [-2.2, -0.6], ncol = 5)
fig.savefig(os.path.join(path_out, f'validation_metrics_root_mortality.png'), dpi = 600.,  bbox_inches = 'tight')
plt.close(fig)