import itertools as it
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import numpy as np
import os
from glob import glob
from datetime import datetime, timedelta

level = 'ambient'
date  = '20211201' # '20211201' (modified), '20211008' (original)
path_data = os.path.join(os.environ['SCRATCHDIR'], 'E3SM', 'output', 
                         date + '_' + level + '_US-SPR_ICB20TRCNPRDCTCBC','run')
path_out = os.path.join(os.environ['PROJDIR'], 'Phenology_ELM', 'output_elm')


###################################################################################################
# Plot the annual cycle of vegetation variables
###################################################################################################
var_list = ['TLAI', 'GPP', 'QVEGE', 'QVEGT',
            'LEAFC_STORAGE', 'LEAFC_XFER', 'LEAFC_STORAGE_TO_XFER', 'LEAFC_XFER_TO_LEAFC']

# -------------------------------------------------------------------------------------------------
# PFT level
# -------------------------------------------------------------------------------------------------
pft_list = {'needleleaf evergreen boreal tree': 2, 
            'needleleaf deciduous boreal tree': 3, 
            'broadleaf deciduous boreal shrub': 11}
clist = ['#1b9e77', '#d95f02', '#7570b3']

flist = glob(os.path.join(path_data, '*h2*.nc'))

fig, axes = plt.subplots(4, 2, figsize = (8, 12), sharex = True)
for ii in range(len(var_list)):
    hr = xr.open_mfdataset(flist, decode_times = True)

    ax = axes.flat[ii]

    count = 0
    nlist = []
    for jj, kk in pft_list.items():
        tvec = [datetime.strptime(x, '%Y-%m-%d %H') for x in \
                hr['time'].indexes['time'].strftime('%Y-%m-%d %H')]
        retype = pd.DataFrame(hr[var_list[ii]][:, kk].values, index = tvec)

        if var_list[ii] == 'TLAI':
            retype = retype.groupby(retype.index.month * 100 + retype.index.day).mean()
        else:
            # # s to day, add up to month, divide by years
            retype = retype.groupby(retype.index.month * 100 + retype.index.day).sum() * 24 * 3600 / len(flist)
        retype.index = range(1,366)

        ax.plot(retype.index, retype.values, '-o', color = clist[count], markersize = 1)
        ax.set_title(var_list[ii])
        ax.set_xlim([0.5, 365.5])
        #ax.set_xticks(range(1,13))

        count += 1
        nlist.append(jj)
    hr.close()
ax.legend(nlist, loc = 'upper right', bbox_to_anchor = (-.6, -.2))
fig.savefig(os.path.join(path_out, date + '_' + level + '_1.png'), dpi = 600., bbox_inches = 'tight')
plt.close()


# -------------------------------------------------------------------------------------------------
# Grid level
# -------------------------------------------------------------------------------------------------
var_list = ['TLAI', 'GPP', 'QVEGE', 'QVEGT']

grid_list = {'hollow': 0, 'hummock': 1}
clist = ['#1b9e77', '#d95f02']

flist = glob(os.path.join(path_data, '*h0*.nc'))

fig, axes = plt.subplots(2, 2, figsize = (6, 6), sharex = True)
for ii in range(len(var_list)):
    hr = xr.open_mfdataset(flist, decode_times = True)

    ax = axes.flat[ii]

    count = 0
    nlist = []
    for jj, kk in grid_list.items():
        tvec = [datetime.strptime(x, '%Y-%m-%d %H') for x in \
                hr['time'].indexes['time'].strftime('%Y-%m-%d %H')]
        retype = pd.DataFrame(hr[var_list[ii]][:, kk].values, index = tvec)
        if var_list[ii] == 'TLAI':
            retype = retype.groupby(retype.index.month * 100 + retype.index.day).mean()
        else:
            # # s to day, add up to month, divide by years
            retype = retype.groupby(retype.index.month * 100 + retype.index.day).sum() * 24 * 3600 / len(flist)
        retype.index = range(1,366)

        ax.plot(retype.index, retype.values, '-o', color = clist[count], markersize = 1)
        ax.set_title(var_list[ii])
        ax.set_xlim([0.5, 365.5])
        #ax.set_xticks(range(1,13))

        count += 1
        nlist.append(jj)

    hr.close()
ax.legend(nlist, loc = 'upper right', bbox_to_anchor = (-.8, -.2))
fig.savefig(os.path.join(path_out, date + '_' + level + '_1-1.png'), dpi = 600., bbox_inches = 'tight')
plt.close()

###################################################################################################
# Plot the annual cycle of air temperature
###################################################################################################
clist = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02']
fig, ax = plt.subplots()
count = 0
for yy in range(2016,2019):
    flist = glob(os.path.join(path_data, '*h1*' + str(yy) + '*.nc'))
    hr = xr.open_mfdataset(flist, decode_times = True)
    ax.plot(range(hr['TSA'].values.shape[0]),
            hr['TSA'].values[:,0] - 273.15, clist[count])
    hr.close()
    ax.set_xlabel('Day of year')
    ax.set_ylabel('Temperature 2m (degC)')
    count += 1
ax.legend([str(yy) for yy in range(2015, 2021)], loc = 'upper right', 
          bbox_to_anchor = (1.2, 1.))
fig.savefig(os.path.join(path_out, date + '_' + level + '_2.png'), dpi = 600.,
            bbox_inches = 'tight')
plt.close()

###################################################################################################
# Plot the diurnal cycle of air temperature and atmospheric incident solar radiation
###################################################################################################
var_list = ['TSA', 'FSDS']
clist = ['#1b9e77', '#d95f02']

flist = os.path.join(path_data, '*h0*.nc')
hr = xr.open_mfdataset(flist, decode_times = True)

fig, axes = plt.subplots(4, 2, sharex = True, figsize = (6, 8))
for ii,jj in it.product(range(len(var_list)), range(4)):
    if var_list[ii] == 'TSA':
        var = hr[var_list[ii]] - 273.15
        unit = 'degC'
    else:
        var = hr[var_list[ii]]
        unit = var.units

    # round time up
    tvec = hr['time'].indexes['time'].to_datetimeindex()
    tvec = [datetime.strptime(x.strftime('%Y-%m-%d %H:%M:%S')[:13],
                              '%Y-%m-%d %H') for x in tvec]

    var = pd.DataFrame(var[:,0].values, index = tvec)
    if jj == 0:
        var = var.loc[(var.index.month == 12) | (var.index.month <= 2), :]
    elif jj == 1:
        var = var.loc[(var.index.month >= 3) & (var.index.month <= 5), :]
    elif jj == 2:
        var = var.loc[(var.index.month >= 6) & (var.index.month <= 8), :]
    else:
        var = var.loc[(var.index.month >= 9) & (var.index.month <= 11), :]
    var = var.groupby(var.index.hour).mean()

    ax = axes[jj,ii]
    ax.plot(var.index, var.values, '-o', color = clist[ii])
    ax.set_title(var_list[ii] + ', ' + unit + ', season = ' + str(jj))
fig.savefig(os.path.join(path_out, date + '_' + level + '_3.png'), dpi = 600., 
            bbox_inches = 'tight')
plt.close()

hr.close()
