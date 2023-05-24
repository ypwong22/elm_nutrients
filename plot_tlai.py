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
#sims_prefix = ['20230510', '20230122', '20230505', '20230509', '20230512']
#sims_names = ['Optim Scheme2', 'Optim Evgr', 'A', 'E', 'H']
sims_prefix = ['20221212', '20230120', '20230510', '20230518']
sims_names = ['Default', 'Optim', 'Optim Scheme 2', 'Optim Scheme 2 Correct input']
clist = ['#de2d26', '#fcfc00', '#3182bd']
co2_levels    = {'ambient' : [6, 20, 13, 8, 17],
                 'elevated': [19, 11, 4, 16, 10]}
pft_list = [2, 3, 11]
pft_names = ['spruce', 'larch', 'shrub']

collect_tbot = read_sims_tair_annual()

hr = xr.open_dataset(os.path.join(path_intrim, 'spruce_validation_data.nc'))

fig, axes = plt.subplots(nrows = 3, ncols = len(sims_prefix), figsize = (15, 8), sharex = True, sharey = True)
fig.subplots_adjust(hspace = 0.1, wspace = 0.05)

for i, (pft, pftname, pftname2) in enumerate(zip(pft_list, pft_names, ["EN", "DN", "SH"])):
    collect_obs = pd.DataFrame(np.nan, index = [2015, 2016, 2017, 2018, 2019, 2020], columns = chamber_list)
    for plot in chamber_list:
        obs = pd.Series(hr['annual_lai'].loc[:, plot, pftname2].values, index = hr['year'])
        collect_obs.loc[:, plot] = obs.loc[collect_obs.index]
    collect_obs = collect_obs

    for j, expr in enumerate(sims_prefix):
        collection_ts = read_extract_sims_ts(expr)

        collect_sim = pd.DataFrame(np.nan, index = [2015, 2016, 2017, 2018, 2019, 2020], columns = chamber_list)
        for plot, plot_name in zip(chamber_list, chamber_list_names):
            # hummock: 0.64, hollow: 0.36
            sim = collection_ts.loc[:, (plot, 'TLAI', pft, 'hummock')] * 0.64 + \
                  collection_ts.loc[:, (plot, 'TLAI', pft, 'hollow')] * 0.36
            sim = sim.groupby(pd.DatetimeIndex([datetime(y,m,1) for y, m in zip(sim.index.year, sim.index.month)])).mean()
            collect_sim.loc[:, plot] = sim.groupby(sim.index.year).max()

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
        #ax.set_ylim([-3, 10.5])
        ax.set_yticks([0, 2, 4, 6, 8, 10])

        if j == 0:
            ax.set_ylabel(pftname)

        if i == (len(pft_list) - 1):
            ax.set_xlabel('Annual mean temperature ($^o$C)')

        if i == 0:
            ax.set_title(sims_names[j])

ax.legend(h, ['OBS_ACO2', 'MOD ACO2', 'OBS ECO2', 'MOD ECO2'], ncol = 4, columnspacing = 1, loc = [-1, -0.5])
fig.savefig(os.path.join(path_out, f'plot_tlai.png'), dpi = 600., bbox_inches = 'tight')
plt.close(fig)
hr.close()