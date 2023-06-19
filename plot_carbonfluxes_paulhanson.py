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
from utils.analysis import *
from datetime import datetime
from gcc_spruce_visualize import *
from glob import glob
from scipy.signal import savgol_filter
from scipy.stats import linregress

q
def plotter(ax, tbot, obs, sims):
    h = [None] * 4
    corr = np.full(4, np.nan)
    pval = np.full(4, np.nan)

    # drop year 2020
    tbot = tbot.drop(2020, axis = 0)
    obs = obs.drop(2020, axis = 0)
    sims = sims.drop(2020, axis = 0)

    a = tbot.loc[:, co2_levels['ambient']].values.reshape(-1)
    b = obs .loc[:, co2_levels['ambient']].values.reshape(-1)
    c = sims.loc[:, co2_levels['ambient']].values.reshape(-1)

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

    a = tbot.loc[:, co2_levels['elevated']].values.reshape(-1)
    b = obs .loc[:, co2_levels['elevated']].values.reshape(-1)
    c = sims.loc[:, co2_levels['elevated']].values.reshape(-1)

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
mpl.rcParams['font.size'] = 14
mpl.rcParams['axes.titlesize'] = 14


#sims_prefix = ['20221212', '20230120', '20230505']  # '20230122', 20230121
#sims_names = ['Default', 'Optim', 'Optim EvgrRoot'] # 'Optim Evgr', 'Optim EvgrRoot'
sims_prefix = ['20221212', '20230120', '20230526', '20230601']
sims_names = ['Default', 'Optim XYS', 'Optim Scheme 2 Correct', 'Optim EvgrRoot']
var_list = ['ANPPtree', 'ANPPshrub', 'NPPmoss', 'BGNPP', 'HR', 'NEE']
co2_levels    = {'ambient' : [6, 20, 13, 8, 17],
                 'elevated': [19, 11, 4, 16, 10]}


collect_tbot = read_sims_tair_annual()


hr = xr.open_dataset(os.path.join(path_intrim, 'spruce_validation_data.nc'))

fig, axes = plt.subplots(nrows = len(var_list), ncols = len(sims_prefix), figsize = (15, 15), sharex = True, sharey = False)
fig.subplots_adjust(hspace = 0.1, wspace = 0.05)
for i, var in enumerate(var_list):
    collect_obs = pd.DataFrame(np.nan, index = [2015, 2016, 2017, 2018, 2019, 2020], columns = chamber_list)
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
        collect_obs.loc[:, plot] = obs.loc[collect_obs.index]

    for j, expr in enumerate(sims_prefix):
        collection_ts = read_extract_sims_ts(expr)

        collect_sim = pd.DataFrame(np.nan, index = [2015, 2016, 2017, 2018, 2019, 2020], columns = chamber_list)
        for plot, plot_name in zip(chamber_list, chamber_list_names):
            if var == 'ANPPtree': 
                # hummock: 0.64, hollow: 0.36
                # pima: 0.36, lala: 0.14
                sim = ((collection_ts.loc[:, (plot, 'AGNPP', 3, 'hummock')] * 0.14 + \
                        collection_ts.loc[:, (plot, 'AGNPP', 2, 'hummock')] * 0.36) * 0.64 + \
                       (collection_ts.loc[:, (plot, 'AGNPP', 3, 'hollow')] * 0.14 + \
                        collection_ts.loc[:, (plot, 'AGNPP', 2, 'hollow')] * 0.36) * 0.36).groupby(collection_ts.index.year).mean() * 86400 * 365
            elif var == 'ANPPshrub':
                # hummock: 0.64, hollow: 0.36
                sim = (collection_ts.loc[:, (plot, 'AGNPP', 11, 'hummock')] * 0.64 + \
                       collection_ts.loc[:, (plot, 'AGNPP', 11, 'hollow' )] * 0.36).groupby(collection_ts.index.year).mean() * 0.25 * 86400 * 365
            elif var == 'NPPmoss':
                # hummock: 0.64, hollow: 0.36
                sim = (collection_ts.loc[:, (plot, 'NPP', 12, 'hummock')] * 0.64 + \
                       collection_ts.loc[:, (plot, 'NPP', 12, 'hollow' )] * 0.36).groupby(collection_ts.index.year).mean() * 0.25 * 86400 * 365

                mossfrac = pd.read_excel(os.path.join(path_input, 'Sphagnum_fraction.xlsx'), index_col = 0, skiprows = 1,
                                         engine = 'openpyxl').drop(['plot','Temp','CO2'], axis = 1)
                mossfrac[2015] = mossfrac[2016]
                mossfrac = mossfrac.drop(2021, axis = 1)
                sim = sim * mossfrac.loc[plot_name] / 100.
            elif var == 'BGNPP':
                # hummock: 0.64, hollow: 0.36
                sim = ((collection_ts.loc[:, (plot, 'BGNPP',  2, 'hummock')] * 0.36 + \
                        collection_ts.loc[:, (plot, 'BGNPP',  3, 'hummock')] * 0.14 + \
                        collection_ts.loc[:, (plot, 'BGNPP', 11, 'hummock')] * 0.25) * 0.64 + \
                       (collection_ts.loc[:, (plot, 'BGNPP',  2, 'hollow')] * 0.36 + \
                        collection_ts.loc[:, (plot, 'BGNPP',  3, 'hollow')] * 0.14 + \
                        collection_ts.loc[:, (plot, 'BGNPP', 11, 'hollow')] * 0.25) * 0.36).groupby(collection_ts.index.year).mean() * 86400 * 365
            elif  var == 'HR':
                sim = - (collection_ts.loc[:, (plot, 'HR', 0, 'hummock')] * 0.64 + \
                         collection_ts.loc[:, (plot, 'HR', 0, 'hollow')] * 0.36).groupby(collection_ts.index.year).mean() * 86400 * 365
            elif  var == 'NEE':
                sim = - (collection_ts.loc[:, (plot, 'NEE', 0, 'hummock')] * 0.64 + \
                         collection_ts.loc[:, (plot, 'NEE', 0, 'hollow')] * 0.36).groupby(collection_ts.index.year).mean() * 86400 * 365

            collect_sim.loc[:, plot] = sim
    
        # import pdb; pdb.set_trace()

        ax = axes[i, j]
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
            ax.set_ylim([-50, 300])
        elif i == 1:
            ax.set_ylim([-50, 300])
        elif i == 2:
            ax.set_ylim([-120, 250])
        elif i == 3:
            ax.set_ylim([-50, 250])
        elif i == 4:
            ax.set_ylim([-900, -200])
        elif i == 5:
            ax.set_ylim([-1000, 700])

        if j == 0:
            ax.set_ylabel(f'{var}\ngC m-2 yr-1')
        else:
            ax.set_yticklabels([])

        if i == (len(var_list) - 1):
            ax.set_xlabel('Annual mean temperature ($^o$C)')

        if i == 0:
            ax.set_title(sims_names[j])

ax.legend(h, ['OBS_ACO2', 'MOD ACO2', 'OBS ECO2', 'MOD ECO2'], loc = (-2., -0.65), ncol = 4)
fig.savefig(os.path.join(path_out, f'plot_aboveground.png'), dpi = 600., bbox_inches = 'tight')
plt.close(fig)
hr.close()