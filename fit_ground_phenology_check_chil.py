"""
Plot winter chilling and GDD v.s. start of observed fine root growth in minirhizotron.

For evergreen, ThermalTime & Alternating model.
For deciduous tree & shrub, Alternating model.
"""
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


lyr_select = '2m'

pheno_obs = read_leaf_sos()['EN']
tsoi_collect, _ = read_obs_tsoi_daily()

observations = pheno_obs.stack().to_frame('doy')
observations.index.names = ['year', 'site_id']
observations = observations.reset_index()


def get_cumu_therm(t1, T, F):
    gdd = tsoi_collect[lyr_select].copy()
    gdd[gdd < T] = 0
    gdd.loc[tsoi_collect[lyr_select].index.dayofyear < t1, :]= 0
    gdd = gdd.loc[gdd.index.year > 2015, :]
    for year in range(2016, 2022):
        gdd.loc[gdd.index.year == year, :] = gdd.loc[gdd.index.year == year, :].cumsum(axis = 0)

    before_leafout = gdd <= F
    return gdd, before_leafout


def get_cumu_alter(a, b, c, threshold, t1):
    chill_days = ((tsoi_collect[lyr_select] < threshold) * 1).copy()
    chill_days.loc[tsoi_collect[lyr_select].index.dayofyear < t1, :] = 0
    chill_days = chill_days.loc[chill_days.index.year > 2015, :]
    for year in range(2016, 2022):
        chill_days.loc[chill_days.index.year == year, :] = chill_days.loc[chill_days.index.year == year, :].cumsum(axis = 0)

    # Accumulated growing degree days from Jan 1
    gdd = tsoi_collect[lyr_select].copy()
    gdd[gdd < threshold] = 0
    gdd.loc[tsoi_collect[lyr_select].index.dayofyear < t1, :]= 0
    gdd = gdd.loc[gdd.index.year > 2015, :]
    for year in range(2016, 2022):
        gdd.loc[gdd.index.year == year, :] = gdd.loc[gdd.index.year == year, :].cumsum(axis = 0)

    # Phenological event happens the first day gdd is > chill_day curve
    chill_day_curve = a + b * np.exp(c * chill_days)
    difference = gdd - chill_day_curve

    # The estimate is equal to the first day that
    # gdd - chill_day_curve > 0
    before_leafout = difference <= 0

    return chill_days, gdd, before_leafout


########################################################################################
# Evergreen tree
########################################################################################
# ThermalTime
params = {'t1': 82.89828273347462, 'T': 0.001321356324646139, 'F': 590.7274692372864}
gdd, before_leafout = get_cumu_therm(**params)
# plot until leaf out
# GDD hasn't even started accumulating when root starts growing
fig, axes = plt.subplots(4, 3, figsize = (12, 12), sharex = True, sharey = True)
for i, plot in enumerate(pheno_obs.columns):
    ax = axes.flat[i]
    temp2 = gdd[plot].loc[before_leafout[plot]]
    temp2 = temp2.loc[temp2.index.year == 2018]
    ax.plot(temp2.index, temp2.values, 'or', ms = 2)
    plt.setp(ax.get_xticklabels(), rotation = 90)

    ax.axvline(datetime(2018, 3, 10)) # start of observed fine root growth, shrub
    ax.axvline(datetime(2018, 4, 4)) # start of observed fine root growth, tree
for i in range(len(pheno_obs.columns), 12):
    ax = axes.flat[i]
    plt.setp(ax.get_xticklabels(), rotation = 90)
fig.savefig(os.path.join(path_out, 'fit_ground_phenology_check_chil_evergreen_tt.png'), dpi = 600., bbox_inches = 'tight')
plt.close(fig)


# Alternating
params = {'a': 281.3048492792021, 'b': 2249.8346532974097, 'c': -0.02412355913509323, 'threshold': 5, 't1': 1}
chill_days, gdd, before_leafout = get_cumu_alter(**params)
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

    ax.axvline(datetime(2018, 3, 10)) # start of observed fine root growth, shrub
    ax.axvline(datetime(2018, 4, 4)) # start of observed fine root growth, tree
for i in range(len(pheno_obs.columns), 12):
    ax = axes.flat[i]
    plt.setp(ax.get_xticklabels(), rotation = 90)
fig.savefig(os.path.join(path_out, 'fit_ground_phenology_check_chil_evergreen_alt.png'), dpi = 600., bbox_inches = 'tight')
plt.close(fig)



"""
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


########################################################################################
# Deciduous tree
########################################################################################
params = {'a': 9, 'b': 2112, 'c': -0.04, 'threshold': 279.5 - 273.15, 't1': 1}
chill_days, gdd, before_leafout = get_cumu_alter(**params)
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

    ax.axvline(datetime(2018, 3, 10)) # start of observed fine root growth, shrub
    ax.axvline(datetime(2018, 4, 4)) # start of observed fine root growth, tree
for i in range(len(pheno_obs.columns), 12):
    ax = axes.flat[i]
    plt.setp(ax.get_xticklabels(), rotation = 90)
fig.savefig(os.path.join(path_out, 'fit_ground_phenology_check_chil_decid_tree_alt.png'), dpi = 600., bbox_inches = 'tight')
plt.close(fig)


########################################################################################
# Deciduous shrub
########################################################################################
params = {'a': 33, 'b': 1388, 'c': -0.02, 'threshold': 279.05 - 273.15, 't1': 1}
chill_days, gdd, before_leafout = get_cumu_alter(**params)
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

    ax.axvline(datetime(2018, 3, 10)) # start of observed fine root growth, shrub
    ax.axvline(datetime(2018, 4, 4)) # start of observed fine root growth, tree
for i in range(len(pheno_obs.columns), 12):
    ax = axes.flat[i]
    plt.setp(ax.get_xticklabels(), rotation = 90)
fig.savefig(os.path.join(path_out, 'fit_ground_phenology_check_chil_decid_shrub_alt.png'), dpi = 600., bbox_inches = 'tight')
plt.close(fig)
