"""
Plot GDD and (if applicable) winter chiling v.s. start of observed fine root growth in minirhizotron.

For evergreen, although the ThermalTime and Uniforc model both fit better than the Alternating model, 
    they do not start accumulating temperature until day (83 - 11) = 73, which is too
    late for the start of fine root growth. 

Therefore, 
    for evergreen, Alternating model.
    for deciduous tree & shrub, Alternating model.
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

    # shift the date forward by 11 days to make 12.21 (winter solstice) -> 1.1
    gdd.loc['2016-02-29', :] = 0.5 * (gdd.loc['2016-02-28', :].values + gdd.loc['2016-03-01', :].values)
    gdd.loc['2020-02-29', :] = 0.5 * (gdd.loc['2020-02-28', :].values + gdd.loc['2020-03-01', :].values)
    gdd = pd.DataFrame(gdd.iloc[:-11, :].values, index = gdd.index[11:], columns = gdd.columns)

    gdd = gdd - T
    gdd[gdd < 0] = 0
    gdd.loc[gdd.index.dayofyear < t1, :]= 0
    gdd = gdd.loc[gdd.index.year > 2015, :]
    for year in range(2016, 2022):
        gdd.loc[gdd.index.year == year, :] = gdd.loc[gdd.index.year == year, :].cumsum(axis = 0)

    before_leafout = gdd <= F

    # shift the dates back
    gdd.index = gdd.index - timedelta(days = 11)
    before_leafout.index = before_leafout.index - timedelta(days = 11)

    return gdd, before_leafout


def get_cumu_alter(a, b, c, threshold, t1):
    temp = tsoi_collect[lyr_select].copy()

    # shift the date forward by 11 days to make 12.21 (winter solstice) -> 1.1
    temp.loc['2016-02-29', :] = 0.5 * (temp.loc['2016-02-28', :].values + temp.loc['2016-03-01', :].values)
    temp.loc['2020-02-29', :] = 0.5 * (temp.loc['2020-02-28', :].values + temp.loc['2020-03-01', :].values)
    temp = pd.DataFrame(temp.iloc[:-11, :].values, 
                        index = temp.index[11:],
                        columns = temp.columns)

    chill_days = ((temp < threshold) * 1).copy()
    chill_days.loc[temp.index.dayofyear < t1, :] = 0
    chill_days = chill_days.loc[chill_days.index.year > 2015, :]
    for year in range(2016, 2022):
        chill_days.loc[chill_days.index.year == year, :] = chill_days.loc[chill_days.index.year == year, :].cumsum(axis = 0)

    # Accumulated growing degree days from winter solstice
    gdd = temp.copy()
    gdd = gdd - threshold
    gdd[gdd < 0] = 0
    gdd.loc[temp.index.dayofyear < t1, :]= 0
    gdd = gdd.loc[gdd.index.year > 2015, :]
    for year in range(2016, 2022):
        gdd.loc[gdd.index.year == year, :] = gdd.loc[gdd.index.year == year, :].cumsum(axis = 0)

    # Phenological event happens the first day gdd is > chill_day curve
    chill_day_curve = a + b * np.exp(c * chill_days)
    difference = gdd - chill_day_curve

    # The estimate is equal to the first day that
    # gdd - chill_day_curve > 0
    before_leafout = difference <= 0

    # shift the dates back
    gdd.index = gdd.index - timedelta(days = 11)
    chill_days.index = chill_days.index - timedelta(days = 11)
    before_leafout.index = before_leafout.index - timedelta(days = 11)

    return chill_days, gdd, before_leafout



########################################################################################
# Check the GDD & chilling accumulation in the Alternating model
########################################################################################








########################################################################################
# Evergreen tree
########################################################################################
# ThermalTime
params = {'t1': 88.96680755982837, 'T': -6.210506141223593, 'F': 995.5839297751559}
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
params = {'a': 199.30927021132482, 'b': 3751.8473351447456, 'c': -0.03728095348592797, 'threshold': 5, 't1': 1}
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


"""# using 5 degree threshold, plot the accumulated chiling days
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
    temp = tsoi_collect[lyr_select][plot].loc[tsoi_collect[lyr_select].index.year > 2015].loc[before_leafout[plot]]
    temp = temp.loc[temp.index.year == 2018]
    ax.plot(temp.index, temp.values, 'o', ms = 2)
ax.grid(axis = 'y')
ax.axvline(datetime(2018, 3, 10))
ax.axvline(datetime(2018, 4, 4))
ax.set_yticks(range(-5, 15, 1))
fig.savefig(os.path.join(path_out, 'fit_ground_phenology_multi_t.png'), dpi = 600., bbox_inches = 'tight')
plt.close(fig)
"""


# Using varying degree threshold, plot the accumulated growing degree days
# when the fine root starts to grow
# Because t1 = 88, the GDD can never start 
thres_list = np.arange(-20, 5.1, 0.5)
gdd_tree_list = []
for t in thres_list:
    params = {'t1': 88.96680755982837, 'T': t, 'F': 995.5839297751559}
    gdd, before_leafout = get_cumu_therm(**params)

    gdd_tree = gdd.loc[datetime(2018, 3, 10), :]

    # average over the chambers
    gdd_tree = pd.Series({'mean': gdd_tree.mean(), 'ambient': gdd_tree.loc['06']})

    gdd_tree_list.append(gdd_tree)
gdd_tree_list = pd.DataFrame(gdd_tree_list)
fig, ax = plt.subplots(figsize = (10, 10))
ax.plot(thres_list, gdd_tree_list['mean'], 'ob', label = 'chamber average')
ax.plot(thres_list, gdd_tree_list['ambient'], 'or', label = 'ambient chamber [06]')
ax.set_ylabel('Accumulated GDD')
ax.set_xlabel('Temperature threshold for GDD')
ax.legend()
fig.savefig(os.path.join(path_out, 'fit_ground_phenology_multi_varyGDD_thermaltime.png'), dpi = 600., bbox_inches = 'tight')
plt.close(fig)


thres_list = np.arange(-10, 5.1, 0.5)
gdd_tree_list = []
chil_tree_list = []
for t in thres_list:
    params = {'a': 199.30927021132482, 'b': 3751.8473351447456, 'c': -0.03728095348592797, 'threshold': t, 't1': 1}
    chill_days, gdd, before_leafout = get_cumu_alter(**params)

    gdd_tree = gdd.loc[datetime(2018, 3, 10), :]
    chil_tree = chill_days.loc[datetime(2018, 3, 10), :]

    # average over the chambers
    gdd_tree = pd.Series({'mean': gdd_tree.mean(), 'ambient': gdd_tree.loc['06']})
    chil_tree = pd.Series({'mean': chil_tree.mean(), 'ambient': chil_tree.loc['06']})

    gdd_tree_list.append(gdd_tree)
    chil_tree_list.append(chil_tree)
gdd_tree_list = pd.DataFrame(gdd_tree_list)
chil_tree_list = pd.DataFrame(chil_tree_list)
fig, axes = plt.subplots(1, 2, figsize = (16, 10))
axes[0].plot(thres_list, gdd_tree_list['mean'], 'ob', label = 'chamber average')
axes[0].plot(thres_list, gdd_tree_list['ambient'], 'or', label = 'ambient chamber [06]')
axes[0].set_ylabel('Accumulated GDD')
axes[0].set_xlabel('Temperature threshold for GDD')
axes[1].plot(thres_list, chil_tree_list['mean'], 'ob', label = 'chamber average')
axes[1].plot(thres_list, chil_tree_list['ambient'], 'or', label = 'ambient chamber [06]')
axes[1].set_ylabel('Accumulated chilling')
axes[1].set_xlabel('Temperature threshold for GDD')
axes[1].legend()
fig.savefig(os.path.join(path_out, 'fit_ground_phenology_multi_varyGDD_alternating.png'), dpi = 600., bbox_inches = 'tight')
plt.close(fig)



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


thres_list = np.arange(-10, 5.1, 0.5)
gdd_tree_list = []
chil_tree_list = []
for t in thres_list:
    params = {'a': 9, 'b': 2112, 'c': -0.04, 'threshold': t, 't1': 1}
    chill_days, gdd, before_leafout = get_cumu_alter(**params)

    gdd_tree = gdd.loc[datetime(2018, 3, 10), :]
    chil_tree = chill_days.loc[datetime(2018, 3, 10), :]

    # average over the chambers
    gdd_tree = pd.Series({'mean': gdd_tree.mean(), 'ambient': gdd_tree.loc['06']})
    chil_tree = pd.Series({'mean': chil_tree.mean(), 'ambient': chil_tree.loc['06']})

    gdd_tree_list.append(gdd_tree)
    chil_tree_list.append(chil_tree)
gdd_tree_list = pd.DataFrame(gdd_tree_list)
chil_tree_list = pd.DataFrame(chil_tree_list)
fig, axes = plt.subplots(1, 2, figsize = (16, 10))
axes[0].plot(thres_list, gdd_tree_list['mean'], 'ob', label = 'chamber average')
axes[0].plot(thres_list, gdd_tree_list['ambient'], 'or', label = 'ambient chamber [06]')
axes[0].set_ylabel('Accumulated GDD')
axes[0].set_xlabel('Temperature threshold for GDD')
axes[1].plot(thres_list, chil_tree_list['mean'], 'ob', label = 'chamber average')
axes[1].plot(thres_list, chil_tree_list['ambient'], 'or', label = 'ambient chamber [06]')
axes[1].set_ylabel('Accumulated chilling')
axes[1].set_xlabel('Temperature threshold for GDD')
axes[1].legend()
fig.savefig(os.path.join(path_out, 'fit_ground_phenology_multi_varyGDD_alternating_decid_tree.png'),
            dpi = 600., bbox_inches = 'tight')
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


thres_list = np.arange(-10, 5.1, 0.5)
gdd_shrub_list = []
chil_shrub_list = []
for t in thres_list:
    params = {'a': 33, 'b': 1388, 'c': -0.02, 'threshold': t, 't1': 1}
    chill_days, gdd, before_leafout = get_cumu_alter(**params)

    gdd_shrub = gdd.loc[datetime(2018, 3, 10), :]
    chil_shrub = chill_days.loc[datetime(2018, 3, 10), :]

    # average over the chambers
    gdd_shrub = pd.Series({'mean': gdd_shrub.mean(), 'ambient': gdd_shrub.loc['06']})
    chil_shrub = pd.Series({'mean': chil_shrub.mean(), 'ambient': chil_shrub.loc['06']})

    gdd_shrub_list.append(gdd_shrub)
    chil_shrub_list.append(chil_shrub)
gdd_shrub_list = pd.DataFrame(gdd_shrub_list)
chil_shrub_list = pd.DataFrame(chil_shrub_list)
fig, axes = plt.subplots(1, 2, figsize = (16, 10))
axes[0].plot(thres_list, gdd_shrub_list['mean'], 'ob', label = 'chamber average')
axes[0].plot(thres_list, gdd_shrub_list['ambient'], 'or', label = 'ambient chamber [06]')
axes[0].set_ylabel('Accumulated GDD')
axes[0].set_xlabel('Temperature threshold for GDD')
axes[1].plot(thres_list, chil_shrub_list['mean'], 'ob', label = 'chamber average')
axes[1].plot(thres_list, chil_shrub_list['ambient'], 'or', label = 'ambient chamber [06]')
axes[1].set_ylabel('Accumulated chilling')
axes[1].set_xlabel('Temperature threshold for GDD')
axes[1].legend()
fig.savefig(os.path.join(path_out, 'fit_ground_phenology_multi_varyGDD_alternating_decid_shrub.png'),
            dpi = 600., bbox_inches = 'tight')
plt.close(fig)