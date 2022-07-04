""" Compare the RMSE of the evaluation against spruce_validation_data.nc across multiple 
    simulations """
import itertools as it
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import numpy as np
import os
from glob import glob
from datetime import datetime, timedelta
from utils.constants import *


var_list = ['co2','ch4','anpp_tree','anpp_shrub','npp_moss', 'bnpp', 'runoff', 'zwt']
title_list = ['Daily CO$_2$ flux', 'Daily CH$_4$ flux', 'Annual ANPP$_{tree}$', 
              'Annual ANPP$_{shrub}$', 'Annual NPP$_{moss}$', 'Annual BNPP',
              'Daily Runoff', 'Daily Water Table']
unit_list = ['gC m-2 day-1', 'gC m-2 day-1', 'gC m-2 yr-1', 'gC m-2 yr-1', 'gC m-2 yr-1',
             'gC m-2 yr-1', 'mm day-1', 'm below hollow surface']

stat = 'RMSE' # 'RMSE','Corr','Bias

ds_collect = pd.DataFrame(np.nan, index = var_list,
                          columns = pd.MultiIndex.from_product([sims_prefix, chamber_list]))

for prefix, chamber in it.product(sims_prefix, chamber_list):
    ds = pd.read_csv(os.path.join('/lustre/haven/proj/UTK0134/Phenology_ELM',
                                  'output_elm', prefix, 'plot' + chamber + '_RMSE.csv'),
                     index_col = 0)
    ds_collect.loc[:, (prefix,chamber)] = ds.loc[var_list, stat]
ds_collect = pd.DataFrame(ds_collect)
ds_collect.stack().to_csv('/lustre/haven/proj/UTK0134/Phenology_ELM/output_elm/' +
                          f'{stat}_summary.csv')


# Plot individual chambers
fig, axes = plt.subplots(3, 3, figsize = (12,12), sharex = True, sharey = False)
fig.subplots_adjust(wspace = 0.2)
for i, var in enumerate(var_list):
    ax = axes.flat[i]
    temp = ds_collect.loc[var, :].unstack().T
    if var == 'runoff':
        temp = temp.dropna(axis = 0, how = 'all')
    # ax.boxplot(temp)
    #for j,col in enumerate(temp.columns):
    #    ax.bar(np.arange(0.5, len(chamber_list) + 0.5) + j/10,
    #           temp.loc[chamber_list, col], align = 'edge', width = 1/10)
    ax.set_prop_cycle('color', ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
                                '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a'])
    ax.plot(temp.loc[chamber_list,:])
    ax.set_title(title_list[i])
    if stat != 'Corr':
        ax.set_ylabel(unit_list[i])
    else:
        ax.set_ylim([-1,1])
    if i >= 5:
        ax.set_xlabel('Chamber')
axes[2,2].axis('off')
ax.legend(sims_names, ncol = 2, loc = [1.2,0.5])
fig.savefig(f'/lustre/haven/proj/UTK0134/Phenology_ELM/output_elm/{stat}_summary.png',
            dpi = 600.)
plt.close()


# Plot the difference from the defaullt model
fig, axes = plt.subplots(2, 4, figsize = (16,6), sharex = True, sharey = False)
fig.subplots_adjust(wspace = 0.2, hspace = 0.23)
for i, var in enumerate(var_list):
    ax = axes.flat[i]
    temp = ds_collect.loc[var, :].unstack().T
    if var == 'runoff':
        temp = temp.dropna(axis = 0, how = 'all')
    if stat == 'Bias':
        # diff in absollute bias
        temp = temp.iloc[:,1:].abs() - np.abs(temp.iloc[:,0].values.reshape(-1,1))
    else:
        temp = temp.iloc[:,1:] - temp.iloc[:,0].values.reshape(-1,1)
    ax.boxplot(temp)
    ax.set_prop_cycle('color', ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
                                '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a'])
    if stat != 'Corr':
        ax.set_title(title_list[i] + '\n' + unit_list[i])
    else:
        ax.set_title(title_list[i])
    ax.set_xticks(range(1,10))
    ax.set_xticklabels(sims_names[1:], rotation = 90)
    if i in [0,4]:
        ax.set_ylabel('$\Delta$ Sim, Default')
fig.savefig(f'/lustre/haven/proj/UTK0134/Phenology_ELM/output_elm/{stat}_summary2.png',
            dpi = 600., bbox_inches = 'tight')
plt.close()


# Plot the difference from the defaullt model, selected
fig, axes = plt.subplots(2, 2, figsize = (9,6), sharex = True, sharey = False)
fig.subplots_adjust(wspace = 0.2, hspace = 0.23)
for i, var in enumerate(var_list[2:-2], 2):
    ax = axes.flat[i-2]
    temp = ds_collect.loc[var, :].unstack().T
    if var == 'runoff':
        temp = temp.dropna(axis = 0, how = 'all')
    if stat == 'Bias':
        # diff in absollute bias
        temp = temp.iloc[:,1:].abs() - np.abs(temp.iloc[:,0].values.reshape(-1,1))
    else:
        temp = temp.iloc[:,1:] - temp.iloc[:,0].values.reshape(-1,1)
    ax.boxplot(temp)
    ax.set_prop_cycle('color', ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
                                '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a'])
    if stat != 'Corr':
        ax.set_title(title_list[i] + '\n' + unit_list[i])
    else:
        ax.set_title(title_list[i])
    ax.set_xticks(range(1,10))
    ax.set_xticklabels(sims_names[1:], rotation = 90)
    if (i-2) in [0,2]:
        ax.set_ylabel('$\Delta$ Sim, Default')
fig.savefig(f'/lustre/haven/proj/UTK0134/Phenology_ELM/output_elm/{stat}_summary3.png',
            dpi = 600., bbox_inches = 'tight')
plt.close()
