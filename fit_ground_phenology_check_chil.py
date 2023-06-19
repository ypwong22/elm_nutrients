"""
Plot GDD and (if applicable) winter chiling v.s. start of observed fine root growth in minirhizotron.


For evergreen, although the Uniforc model fits best, it requires a low 
    base temperature (< -15 degC) for accumulation to reach 66 (same as leaf). 

Also, for all the PFTs, It is not feasible to adjust the GDD threshold (requires negative values)
   or the temperature threshold (requires < -7 degC) within reasonable ranges, to model
   roots out timing. 

Therefore, use a chiling threshold in the alternating model for all 3 PFTs
    (cumulated number of days < 4 degC)
    evergreen: 100
    deciduous: 100
    shrub: 80

    In some years, this chiling threshold may not be reached because weather is warm. 
    Then let fine root grow at the latest at the same time as leaf out. 
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
    chill_day_curve.index = chill_day_curve.index - timedelta(days = 11)
    before_leafout.index = before_leafout.index - timedelta(days = 11)

    return chill_day_curve, chill_days, gdd, before_leafout


def get_cumu_uniforc(t1, F, b, c):
    gdd = tsoi_collect[lyr_select].copy()

    # shift the date forward by 11 days to make 12.21 (winter solstice) -> 1.1
    gdd.loc['2016-02-29', :] = 0.5 * (gdd.loc['2016-02-28', :].values + gdd.loc['2016-03-01', :].values)
    gdd.loc['2020-02-29', :] = 0.5 * (gdd.loc['2020-02-28', :].values + gdd.loc['2020-03-01', :].values)
    gdd = pd.DataFrame(gdd.iloc[:-11, :].values, index = gdd.index[11:], columns = gdd.columns)

    gdd = 1 / (1 + np.exp(b * (gdd - c)))
    gdd.loc[gdd.index.dayofyear < t1, :]= 0
    gdd = gdd.loc[gdd.index.year > 2015, :]
    for year in range(2016, 2022):
        gdd.loc[gdd.index.year == year, :] = gdd.loc[gdd.index.year == year, :].cumsum(axis = 0)

    before_leafout = gdd <= F

    # shift the dates back
    gdd.index = gdd.index - timedelta(days = 11)
    before_leafout.index = before_leafout.index - timedelta(days = 11)

    return gdd, before_leafout


def check_gdd_chil_root(year, month, date, params):
    # thres_list = np.arange(-15, 5.1, 0.5)
    if params['a'] < 100:
        thres_list = np.linspace(-20, 51, 21)
    else:
        thres_list = np.linspace(-800, 200, 41)
    gdd_tree_list = []
    chil_tree_list = []
    for t in thres_list:
        params = params.copy()
        params['a'] = t
        chill_day_curve, gdd, _ = get_cumu_alter(**params)

        gdd_tree = gdd.loc[datetime(year, month, date), :]
        chil_tree = chill_day_curve.loc[datetime(year, month, date), :]

        # average over the chambers
        gdd_tree = pd.Series({'mean': gdd_tree.mean(), 'ambient': gdd_tree.loc['06']})
        chil_tree = pd.Series({'mean': chil_tree.mean(), 'ambient': chil_tree.loc['06']})

        gdd_tree_list.append(gdd_tree)
        chil_tree_list.append(chil_tree)
    gdd_tree_list = pd.DataFrame(gdd_tree_list)
    chil_tree_list = pd.DataFrame(chil_tree_list)
    fig, ax = plt.subplots(figsize = (8, 10))
    ax.plot(thres_list, gdd_tree_list['mean'] - chil_tree_list['mean'], 'ob', label = 'chamber average')
    ax.plot(thres_list, gdd_tree_list['ambient'] - chil_tree_list['ambient'], 'or', label = 'ambient chamber [06]')
    ax.set_ylabel('Accumulated GDD \u2212 chil threshold at root-out')
    # ax.set_xlabel('Temperature threshold for GDD')
    ax.set_xlabel('parameter a')
    ax.legend()
    return fig


def check_chil_root(year, month, date, params):
    _, chill_days, _, _ = get_cumu_alter(**params)

    chil_tree = chill_days.loc[datetime(year, month, date), :]

    # plot for each chamber
    fig, ax = plt.subplots(figsize = (8, 10))
    ax.plot(range(len(chil_tree)), chil_tree, 'ob', label = 'chamber average')
    ax.set_xticks(range(len(chil_tree)))
    ax.set_xticklabels(chil_tree.index)
    ax.set_ylabel('Accumulated chiling at root-out')
    ax.set_xlabel('Chamber ID')
    ax.legend()
    return fig


def check_uniforc_root(year, month, date, params):
    # thres_list = np.arange(-15, 5.1, 0.5)
    thres_list = np.linspace(params['F'] * 0.2, params['F'], 11)
    gdd_tree_list = []
    for t in thres_list:
        params = params.copy()
        params['F'] = t
        gdd, _ = get_cumu_uniforc(**params)

        gdd_tree = gdd.loc[datetime(year, month, date), :]

        # minus threshold
        gdd_tree = gdd_tree - params['F']

        # average over the chambers
        gdd_tree = pd.Series({'mean': gdd_tree.mean(), 'ambient': gdd_tree.loc['06']})

        gdd_tree_list.append(gdd_tree)
    gdd_tree_list = pd.DataFrame(gdd_tree_list)
    fig, ax = plt.subplots(figsize = (8, 10))
    ax.plot(thres_list, gdd_tree_list['mean'], 'ob', label = 'chamber average')
    ax.plot(thres_list, gdd_tree_list['ambient'], 'or', label = 'ambient chamber [06]')
    ax.set_ylabel('Accumulated GDD \u2212 threshold at root-out')
    # ax.set_xlabel('Temperature threshold for GDD')
    ax.set_xlabel('GDD threshold')
    ax.legend()
    return fig


########################################################################################
# Check the GDD & chilling accumulation in the Alternating model
########################################################################################
"""
# Varying the GDD threshold (a/F) or temperature threshold for GDD (threshold) do not work well for shrubs

# Evergreen tree
# params = {'a': 199.30927021132482, 'b': 3751.8473351447456, 'c': -0.03728095348592797, 'threshold': 5, 't1': 1}
params = {'a': 199.30927021132482, 'b': 3751.8473351447456, 'c': -0.03728095348592797, 'threshold': -5, 't1': 1}
fig = check_gdd_chil_root(2018, 4, 4, params)
fig.savefig(os.path.join(path_out, 'fit_ground_phenology_multi_varyGDD_alternating_evgr.png'), dpi = 600.,
            bbox_inches = 'tight')
plt.close(fig)

# params = {'F': 66.94329645213861, 'b': -0.04412234589076469, 'c': 6.848386620467728, 't1': 1}
params = {'F': 66.94329645213861, 'b': -0.04412234589076469, 'c': -5, 't1': 1}
fig = check_uniforc_root(2018, 4, 4, params)
fig.savefig(os.path.join(path_out, 'fit_ground_phenology_multi_varyGDD_uniforc_evgr.png'), dpi = 600.,
            bbox_inches = 'tight')
plt.close(fig)

# Deciduous tree
# params = {'a': 9, 'b': 2112, 'c': -0.04, 'threshold': 279.5 - 273.15, 't1': 1}
params = {'a': 9, 'b': 2112, 'c': -0.04, 'threshold': -5, 't1': 1}
fig = check_gdd_chil_root(2018, 4, 4, params)
fig.savefig(os.path.join(path_out, 'fit_ground_phenology_multi_varyGDD_alternating_decid_tree.png'), dpi = 600.,
            bbox_inches = 'tight')
plt.close(fig)

# Deciduous shrub
params = {'a': 33, 'b': 1388, 'c': -0.02, 'threshold': -5, 't1': 1}
fig = check_gdd_chil_root(2018, 3, 10, params)
fig.savefig(os.path.join(path_out, 'fit_ground_phenology_multi_varyGDD_alternating_decid_shrub.png'), dpi = 600.,
            bbox_inches = 'tight')
plt.close(fig)
"""


# Evergreen tree
params = {'a': 199.30927021132482, 'b': 3751.8473351447456, 'c': -0.03728095348592797, 'threshold': 5, 't1': 1}
fig = check_chil_root(2018, 4, 4, params)
fig.savefig(os.path.join(path_out, 'fit_ground_phenology_multi_chil_evgr.png'), dpi = 600.,
            bbox_inches = 'tight')
plt.close(fig)


# Deciduous tree
params = {'a': 9, 'b': 2112, 'c': -0.04, 'threshold': 279.5 - 273.15, 't1': 1}
fig = check_chil_root(2018, 4, 4, params)
fig.savefig(os.path.join(path_out, 'fit_ground_phenology_multi_chil_decid_tree.png'), dpi = 600.,
            bbox_inches = 'tight')
plt.close(fig)


# Deciduous shrub
params = {'a': 33, 'b': 1388, 'c': -0.02, 'threshold': 279.05 - 273.15, 't1': 1}
fig = check_chil_root(2018, 3, 10, params)
fig.savefig(os.path.join(path_out, 'fit_ground_phenology_multi_chil_decid_shrub.png'), dpi = 600.,
            bbox_inches = 'tight')
plt.close(fig)


########################################################################################
# Use selected temperature threshold to see what happens
########################################################################################
# Evergreen tree
params = {'a': 199.30927021132482, 'b': 3751.8473351447456, 'c': -0.03728095348592797, 'threshold': 5, 't1': 1}
_, _, _, before_leafout = get_cumu_alter(**params)
params = {'a': 199.30927021132482, 'b': 3751.8473351447456, 'c': -0.03728095348592797, 'threshold': 5, 't1': 1}
_, chill_days, gdd, before_rootout = get_cumu_alter(**params)
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

    ax.axvline(datetime(2018, 4, 4)) # start of observed fine root growth, tree

    temp3 = before_rootout[plot].loc[(before_rootout[plot].index.year == 2018) & before_leafout[plot]]
    ax.axvline(temp2.index[temp3][-1], ls = '--', color = 'k') # estimated fine root growth
# all chamber average
ax = axes.flat[i+1]
temp = chill_days.mean(axis = 1)
temp = temp.loc[(temp.index.year == 2018) & (temp.index.month < 5) ]
ax.plot(temp.index, temp.values, 'ob', ms = 4)
temp2 = gdd.mean(axis = 1)
temp2 = temp2.loc[(temp2.index.year == 2018) & (temp2.index.month < 5) ]
ax.plot(temp2.index, temp2.values, 'or', ms = 4)
ax.axvline(datetime(2018, 4, 4)) # start of observed fine root growth, tree
for i in range(len(pheno_obs.columns), 12):
    ax = axes.flat[i]
    plt.setp(ax.get_xticklabels(), rotation = 90)
fig.savefig(os.path.join(path_out, 'fit_ground_phenology_check_chil_evergreen_alt.png'), dpi = 600., bbox_inches = 'tight')
plt.close(fig)


"""# Evergreen tree
params = {'F': 66.94329645213861, 'b': -0.04412234589076469, 'c': 6.848386620467728, 't1': 1}
_, before_leafout = get_cumu_uniforc(**params)
params = {'F': 66.94329645213861, 'b': -0.04412234589076469, 'c': -6, 't1': 1}
gdd, before_rootout = get_cumu_uniforc(**params)
# plot until leaf out
fig, axes = plt.subplots(4, 3, figsize = (12, 12), sharex = True, sharey = True)
for i, plot in enumerate(pheno_obs.columns):
    ax = axes.flat[i]
    temp2 = gdd[plot].loc[before_leafout[plot]]
    temp2 = temp2.loc[temp2.index.year == 2018]
    ax.plot(temp2.index, temp2.values, 'or', ms = 2)
    plt.setp(ax.get_xticklabels(), rotation = 90)

    ax.axvline(datetime(2018, 4, 4)) # start of observed fine root growth, tree

    temp3 = before_rootout[plot].loc[(before_rootout[plot].index.year == 2018) & before_leafout[plot]]
    ax.axvline(temp2.index[temp3][-1], ls = '--', color = 'k') # estimated fine root growth
# all chamber average
ax = axes.flat[i+1]
temp2 = gdd.mean(axis = 1)
temp2 = temp2.loc[(temp2.index.year == 2018) & (temp2.index.month < 5) ]
ax.plot(temp2.index, temp2.values, 'or', ms = 4)
ax.axvline(datetime(2018, 4, 4)) # start of observed fine root growth, tree
for i in range(len(pheno_obs.columns), 12):
    ax = axes.flat[i]
    plt.setp(ax.get_xticklabels(), rotation = 90)
fig.savefig(os.path.join(path_out, 'fit_ground_phenology_check_gdd_evergreen_uniforc.png'), dpi = 600., bbox_inches = 'tight')
plt.close(fig)
"""


# Deciduous tree
params = {'a': 9, 'b': 2112, 'c': -0.04, 'threshold': 279.5 - 273.15, 't1': 1}
_, _, _, before_leafout = get_cumu_alter(**params)
params = {'a': 9, 'b': 2112, 'c': -0.04, 'threshold': 279.5 - 273.15, 't1': 1}
_, chill_days, gdd, before_rootout = get_cumu_alter(**params)
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

    ax.axvline(datetime(2018, 4, 4)) # start of observed fine root growth, tree

    temp3 = before_rootout[plot].loc[(before_rootout[plot].index.year == 2018) & before_leafout[plot]]
    ax.axvline(temp2.index[temp3][-1], ls = '--', color = 'k') # estimated fine root growth
# all chamber average
ax = axes.flat[i+1]
temp = chill_days.mean(axis = 1)
temp = temp.loc[(temp.index.year == 2018) & (temp.index.month < 5) ]
ax.plot(temp.index, temp.values, 'ob', ms = 4)
temp2 = gdd.mean(axis = 1)
temp2 = temp2.loc[(temp2.index.year == 2018) & (temp2.index.month < 5) ]
ax.plot(temp2.index, temp2.values, 'or', ms = 4)
ax.axvline(datetime(2018, 4, 4)) # start of observed fine root growth, tree
for i in range(len(pheno_obs.columns), 12):
    ax = axes.flat[i]
    plt.setp(ax.get_xticklabels(), rotation = 90)
fig.savefig(os.path.join(path_out, 'fit_ground_phenology_check_chil_decid_tree_alt.png'), dpi = 600., bbox_inches = 'tight')
plt.close(fig)

# Deciduous shrub
params = {'a': 33, 'b': 1388, 'c': -0.02, 'threshold': 279.05 - 273.15, 't1': 1}
_, _, _, before_leafout = get_cumu_alter(**params)
params = {'a': 33, 'b': 1388, 'c': -0.02, 'threshold': 279.05 - 273.15, 't1': 1}
_, chill_days, gdd, before_rootout = get_cumu_alter(**params)
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

    temp3 = before_rootout[plot].loc[(before_rootout[plot].index.year == 2018) & before_leafout[plot]]
    ax.axvline(temp2.index[temp3][-1], ls = '--', color = 'k') # estimated fine root growth
# all chamber average
ax = axes.flat[i+1]
temp = chill_days.mean(axis = 1)
temp = temp.loc[(temp.index.year == 2018) & (temp.index.month < 5) ]
ax.plot(temp.index, temp.values, 'ob', ms = 4)
temp2 = gdd.mean(axis = 1)
temp2 = temp2.loc[(temp2.index.year == 2018) & (temp2.index.month < 5) ]
ax.plot(temp2.index, temp2.values, 'or', ms = 4)
ax.axvline(datetime(2018, 3, 10)) # start of observed fine root growth, shrub
for i in range(len(pheno_obs.columns), 12):
    ax = axes.flat[i]
    plt.setp(ax.get_xticklabels(), rotation = 90)
fig.savefig(os.path.join(path_out, 'fit_ground_phenology_check_chil_decid_shrub_alt.png'), dpi = 600., bbox_inches = 'tight')
plt.close(fig)