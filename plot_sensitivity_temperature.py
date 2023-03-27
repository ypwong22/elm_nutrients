""" Plot the sensitivity of PH's variables to temperature. """
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import itertools as it
import xarray as xr
import matplotlib as mpl
from utils.constants import *
from utils.paths import *
from utils.plotting import *
from datetime import datetime
from gcc_spruce_visualize import *
from glob import glob
from scipy.signal import savgol_filter
from scipy.stats import linregress


def plotter(ax, tbot, obs, sims):
    h = [None] * 4
    corr = np.full(4, np.nan)
    pval = np.full(4, np.nan)

    a = tbot.loc[co2_levels['ambient'], :].values.reshape(-1)
    b = obs .loc[co2_levels['ambient'], :].values.reshape(-1)
    c = sims.loc[co2_levels['ambient'], :].values.reshape(-1)

    h[0], = ax.plot(a, b, 'ob', markerfacecolor = 'b')
    h[1], = ax.plot(a, c, 'or', markerfacecolor = 'r')

    res = linregress(a[~np.isnan(b)], b[~np.isnan(b)])
    corr[0] = res.slope
    pval[0] = res.pvalue
    ax.plot([4, 17], res.intercept + res.slope * np.array([4, 17]), '-b')

    res = linregress(a, c)
    corr[1] = res.slope
    pval[1] = res.pvalue
    ax.plot([4, 17], res.intercept + res.slope * np.array([4, 17]), '-r')

    a = tbot.loc[co2_levels['elevated'], :].values.reshape(-1)
    b = obs .loc[co2_levels['elevated'], :].values.reshape(-1)
    c = sims.loc[co2_levels['elevated'], :].values.reshape(-1)

    h[2], = ax.plot(a, b, 'ob', markerfacecolor = 'none')
    h[3], = ax.plot(a, c, 'or', markerfacecolor = 'none')

    res = linregress(a[~np.isnan(b)], b[~np.isnan(b)])
    corr[2] = res.slope
    pval[2] = res.pvalue
    ax.plot([4, 17], res.intercept + res.slope * np.array([4, 17]), ':b')

    res = linregress(a, c)
    corr[3] = res.slope
    pval[3] = res.pvalue
    ax.plot([4, 17], res.intercept + res.slope * np.array([4, 17]), ':r')

    return h, corr, pval


# for slides
mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.titlesize'] = 16


sims_prefix = ['20221212', '20221231', '20230212']
sims_names = ['Default', 'Evergreen', 'EvgrRoot']
clist = ['#de2d26', '#fcfc00', '#3182bd']
var_list = ['ANPPtree', 'ANPPshrub', 'NPPmoss', 'BNPP', 'HR', 'NEE']
co2_levels    = {'ambient' : [6, 20, 13, 8, 17],
                 'elevated': [19, 11, 4, 16, 10]}


collect_tbot = pd.DataFrame(np.nan, index = chamber_list, columns = [2015, 2016, 2017, 2018, 2019, 2020])
for plot in chamber_list:
    path_data = os.path.join(os.environ['PROJDIR'], 'E3SM', 'output', f'20221212_plot{plot:02d}_US-SPR_ICB20TRCNPRDCTCBC', 'run')
    flist = sorted(glob(os.path.join(path_data, '*.h1.*.nc')))
    hr2 = xr.open_mfdataset(flist)
    tbot = (hr2['TBOT'][:-1, 0] * 0.64 + hr2['TBOT'][:-1, 1] * 0.36).resample({'time': '1Y'}).mean()
    tbot = tbot - 273.15
    hr2.close()
    collect_tbot.loc[plot, :] = tbot


hr = xr.open_dataset(os.path.join(path_intrim, 'spruce_validation_data.nc'))


for i, var in enumerate(var_list):
    
    if (i == 0) | (i == 3):
        fig, axes = plt.subplots(nrows = 3, ncols = len(sims_prefix), figsize = (15, 7.5), sharex = True, sharey = False)
        fig.subplots_adjust(hspace = 0.1, wspace = 0.05)
        count = 0

    collect_obs = pd.DataFrame(np.nan, index = chamber_list, columns = [2015, 2016, 2017, 2018, 2019, 2020])
    for plot in chamber_list:
        if var == 'ANPPtree': 
            obs = pd.Series(hr['annual_anpp_tree'].loc[:, plot].values, index = hr['year'])
        elif var == 'ANPPshrub':
            obs = pd.Series(hr['annual_anpp_shrub'].loc[:, plot].values, index = hr['year'])
        elif var == 'NPPmoss':
            obs = pd.Series(hr['annual_npp_moss'].loc[:, plot].values, index = hr['year'])
        elif var == 'BGNPP':
            obs = pd.Series(hr['annual_bnpp'].loc[:, plot].values, index = hr['year'])
        elif  var == 'HR':
            obs = pd.Series(hr['annual_rh'].loc[:, plot].values, index = hr['year'])
        elif  var == 'NEE':
            obs = pd.Series(hr['annual_nee'].loc[:, plot].values, index = hr['year'])
        collect_obs.loc[plot, :] = obs.loc[collect_obs.columns]


    for j, expr in enumerate(sims_prefix):
        collect_sim = pd.DataFrame(np.nan, index = chamber_list, columns = [2015, 2016, 2017, 2018, 2019, 2020])
        for plot, plot_name in zip(chamber_list, chamber_list_names):
            data = pd.read_csv(os.path.join(path_out_elm, expr, f'plot{plot}_ts_extract.csv'), parse_dates = True, index_col = 0, header = [0, 1])
            if var == 'ANPPtree': 
                # hummock: 0.64, hollow: 0.36
                # pima: 0.36, lala: 0.14
                sim = ((data.loc[:, ('hummock', 'AGNPP_lala')] * 0.14 + \
                        data.loc[:, ('hummock', 'AGNPP_pima')] * 0.36) * 0.64 + \
                    (data.loc[:, ('hollow', 'AGNPP_lala')] * 0.14 + \
                        data.loc[:, ('hollow', 'AGNPP_pima')] * 0.36) * 0.36).groupby(data.index.year).sum()
            elif var == 'ANPPshrub':
                # hummock: 0.64, hollow: 0.36
                sim = (data.loc[:, ('hummock', 'AGNPP_shrub')] * 0.64 + \
                    data.loc[:, ('hollow', 'AGNPP_shrub')] * 0.36).groupby(data.index.year).sum() * 0.25
            elif var == 'NPPmoss':
                # hummock: 0.64, hollow: 0.36
                sim = (data.loc[:, ('hummock', 'NPP_moss')] * 0.64 + \
                    data.loc[:, ('hollow', 'NPP_moss')] * 0.36).groupby(data.index.year).sum()

                mossfrac = pd.read_excel('Sphagnum_fraction.xlsx', index_col = 0, skiprows = 1,
                                        engine = 'openpyxl').drop(['plot','Temp','CO2'], axis = 1)
                mossfrac[2015] = mossfrac[2016]
                mossfrac = mossfrac.drop(2021, axis = 1)
                sim = sim * mossfrac.loc[plot_name] / 100.
            elif var == 'BGNPP':
                # hummock: 0.64, hollow: 0.36
                sim = ((data.loc[:, ('hummock', 'BGNPP_pima')] * 0.36 + \
                        data.loc[:, ('hummock', 'BGNPP_lala')] * 0.14 + \
                        data.loc[:, ('hummock', 'BGNPP_shrub')] * 0.25) * 0.64 + \
                    (data.loc[:, ('hollow', 'BGNPP_pima')] * 0.36 + \
                        data.loc[:, ('hollow', 'BGNPP_lala')] * 0.14 + \
                        data.loc[:, ('hollow', 'BGNPP_shrub')] * 0.25) * 0.36).groupby(data.index.year).sum()            
            elif  var == 'HR':
                sim = - (data.loc[:, ('hummock', 'HR')] * 0.64 + data.loc[:, ('hollow', 'HR')] * 0.36).groupby(data.index.year).sum()
            elif  var == 'NEE':
                sim = - (data.loc[:, ('hummock', 'NEE')] * 0.64 + data.loc[:, ('hollow', 'NEE')] * 0.36).groupby(data.index.year).sum()
            collect_sim.loc[plot, :] = sim

        ax = axes[count, j]
        h, corr, pval = plotter(ax, collect_tbot, collect_obs, collect_sim)

        clist = ['b', 'r', 'b', 'r']
        for k in range(4):
            if pval[k] <= 0.05:
                fontweight = 'bold'
            else:
                fontweight = 'normal'
            t = ax.text(0.05 + 0.19*k, 0.07, f'{corr[k]:.2f}', fontweight = fontweight, transform = ax.transAxes, fontsize = 9, color = clist[k])
            if k < 2:
                t.set_bbox(dict(facecolor = clist[k], alpha = 0.2, edgecolor = clist[k]))
            else:
                t.set_bbox(dict(facecolor = 'w', alpha = 0.2, edgecolor = clist[k]))

        if i == 0:
            ax.set_ylim([-20, 300])
        elif i == 1:
            ax.set_ylim([-20, 300])
        elif i == 2:
            ax.set_ylim([-60, 300])
        elif i == 3:
            ax.set_ylim([-80, 200])
        elif i == 4:
            ax.set_ylim([-900, -200])
        elif i == 5:
            ax.set_ylim([-1000, 700])

        if j == 0:
            ax.set_ylabel(var)
        else:
            ax.set_yticklabels([])

        if count == 2:
            ax.set_xlabel('TBOT ($^o$C)')

        if count == 0:
            ax.set_title(sims_names[j])

        if ((i == 2) | (i == (len(var_list) - 1))) and (j == 2):
            ax.legend(h, ['OBS aCO2', 'SIM aCO2', 'OBS eCO2', 'SIM eCO2'], loc = (-2.2, -0.65), ncol = 4)
            fig.savefig(os.path.join(path_out, f'validation_sensitivity_{i}.png'), dpi = 600., bbox_inches = 'tight')
            plt.close(fig)

    count = count + 1

hr.close()
