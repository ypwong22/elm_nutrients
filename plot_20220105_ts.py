import itertools as it
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import numpy as np
import os
from glob import glob
from datetime import datetime, timedelta
from scipy.signal import savgol_filter
from scipy.stats import pearsonr
from optparse import OptionParser


parser = OptionParser()
parser.add_option("--level", dest="level", default='', \
                  help = 'Choose the chamber to analyze.')
parser.add_option("--read", dest="read", action='store_true', \
                  default = False, help = 'Read from existing csv.')
parser.add_option('--add', dest = 'add', default = '', \
                  help = 'Update the file to add a few variables.')
(options,args) = parser.parse_args()


level = options.level
read  = options.read
if options.add == '':
    add   = []
else:
    add   = options.add.split(',')


print(level, read, add)

###################################################################################################

path_out = os.path.join(os.environ['PROJDIR'], 'Phenology_ELM', 'output_elm')
# HR_vr and H2OSOI are vertically resolved; PARVEGLN is not filled
var_list = ['TLAI', 'LAISHA', 'LAISUN', 'FSA', 'PSNSUN', 'PSNSHA', 
            'QVEGT', 'FSH_G', 'FSH_V', 'FSH', 'FGR', 'FSM', 'TBOT',
            'GPP', 'FROOTC', 'AR', 'LITFALL', 'NEE', 'SOILC_HR', 'HR_vr',
            'FROST_TABLE', 'H2OSOI', 'QINTR', 'FSDS', 'DENIT', 'CH4PROD']


if not read:
    data_collect = {}
    for date,model in zip(['20211201', '20211008'], ['modified', 'original']):
        path_data = os.path.join(os.environ['SCRATCHDIR'], 'E3SM', 'output', 
                                 date + '_' + level + '_US-SPR_ICB20TRCNPRDCTCBC','run')
        flist = glob(os.path.join(path_data, '*h0*.nc'))
        hr = xr.open_mfdataset(flist, decode_times = True)
        for var in var_list:
            tvec = [datetime.strptime(x, '%Y-%m-%d %H') for x in \
                    hr['time'].indexes['time'].strftime('%Y-%m-%d %H')]
    
            for kk, loc in zip([0,1], ['hummock', 'hollow']):
                if var in ['H2OSOI', 'HR_vr']:
                    thickgrnd = np.insert(np.insert(hr['levgrnd'].values.reshape(1, -1), 0, 
                                                    - hr['levgrnd'].values[0]), -1, 40)
                    thickgrnd = np.diff((thickgrnd[1:] + thickgrnd[:-1])/ 2)
                    temp = np.sum(hr[var][:, :, kk].values * thickgrnd, axis = 1) / \
                           thickgrnd.sum()
                    data_collect[(loc, model, var + ';' + hr[var].attrs['units'])] = \
                        pd.Series(temp, index = tvec)
                else:
                    data_collect[(loc, model, var + ';' + hr[var].attrs['units'])] = \
                        pd.Series(hr[var][:, kk].values, index = tvec)
        hr.close()
    data_collect = pd.DataFrame(data_collect)
    data_collect.to_csv(os.path.join(path_out, '20220105', level + '.csv'))

if len(add) > 0:
    data_collect = pd.read_csv(os.path.join(path_out, '20220105', level + '.csv'),
                               header = [0,1,2],
                               index_col = 0, parse_dates = True)
    for date,model in zip(['20211201', '20211008'], ['modified', 'original']):
        for kk, loc in zip([0,1], ['hummock', 'hollow']):
            path_data = os.path.join(os.environ['SCRATCHDIR'], 'E3SM', 'output', 
                                     date + '_' + level + '_US-SPR_ICB20TRCNPRDCTCBC','run')
            flist = glob(os.path.join(path_data, '*h0*.nc'))
            hr = xr.open_mfdataset(flist, decode_times = True)
            for var in add:
                tvec = [datetime.strptime(x, '%Y-%m-%d %H') for x in \
                        hr['time'].indexes['time'].strftime('%Y-%m-%d %H')]
                for kk, loc in zip([0,1], ['hummock', 'hollow']):
                    if var in ['H2OSOI', 'HR_vr']:
                        thickgrnd = np.insert(np.insert(hr['levgrnd'].values.reshape(1, -1), 0, 
                                                        - hr['levgrnd'].values[0]), -1, 40)
                        thickgrnd = np.diff((thickgrnd[1:] + thickgrnd[:-1])/ 2)
                        temp = np.sum(hr[var][:, :, kk].values * thickgrnd, axis = 1) / \
                               thickgrnd.sum()
                        data_collect[(loc, model, var + ';' + hr[var].attrs['units'])] = \
                            pd.Series(temp, index = tvec)
                    else:
                        data_collect[(loc, model, var + ';' + hr[var].attrs['units'])] = \
                            pd.Series(hr[var][:, kk].values, index = tvec)
            hr.close()
    data_collect.to_csv(os.path.join(path_out, '20220105', level + '.csv'))


data_collect = pd.read_csv(os.path.join(path_out, '20220105', level + '.csv'), header = [0,1,2],
                           index_col = 0, parse_dates = True)
data_collect = data_collect.groupby(data_collect.index.date).mean()
data_collect.to_csv(os.path.join(path_out, '20220105', level + '_daily.csv'))
for var in var_list + add:
    fig, axes = plt.subplots(2, 1, figsize = (12, 6))
    for ind, grid in enumerate(['hummock', 'hollow']):
        ax = axes.flat[ind]
        which = np.array([(x[0] == grid) and (x[1] == 'modified') and (x[2].split(';')[0] == var) \
                          for x in data_collect.columns])
        which2 = np.array([(x[0] == grid) and (x[1] == 'original') \
                           and (x[2].split(';')[0] == var) for x in data_collect.columns])

        ax.plot(data_collect.index, savgol_filter(data_collect.loc[:, which].values[:,0],
                                                  61, 3),
                '-b', lw = 0.5)
        ax.plot(data_collect.index, savgol_filter(data_collect.loc[:, which2].values[:,0],
                                                  61, 3),
                '-r', lw = 0.5)
        ax.set_title(grid + ' ' + data_collect.columns.values[which][0][2])
        ax.legend(['modified', 'original'])
    fig.savefig(os.path.join(path_out, '20220105', level + '_' + var + '.png'), dpi = 600., 
                bbox_inches = 'tight')
    plt.close(fig)


for var in ['TLAI','FSDS']:
    fig, axes = plt.subplots(2, 2, figsize = (12, 12))
    for ind,grid in enumerate(['hummock', 'hollow']):
        for ind2,model in enumerate(['modified', 'original']):
            ax = axes[ind, ind2]
            if var == 'TLAI':
                x = data_collect[(grid, model, 'TLAI;none')].values
            else:
                x = data_collect[(grid, model, 'FSDS;W/m^2')].values
            y = data_collect[(grid, model, 'GPP;gC/m^2/s')].values
            ax.plot(x, y, 'o', markersize = 3)

            r, pval = pearsonr(x, y)
            ax.text(0.05, 0.85, '{:.4f}, p = {:.3f}'.format(r, pval),
                    transform = ax.transAxes)

    fig.savefig(os.path.join(path_out, '20220105', level + '_' + var + \
                             '_photosyntheis.png'), dpi = 600., 
                bbox_inches = 'tight')
    plt.close(fig)


fig, axes = plt.subplots(2, 2, figsize = (12, 12), sharex = True, sharey = True)
for ind,grid in enumerate(['hummock', 'hollow']):
    for ind2,model in enumerate(['modified', 'original']):
        ax = axes[ind, ind2]
        x = data_collect[(grid, model, 'TLAI;none')].values
        y = data_collect[(grid, model, 'FSDS;W/m^2')].values
        z = data_collect[(grid, model, 'GPP;gC/m^2/s')].values
        cf = ax.scatter(x, y, c = z, s = 3, cmap = 'Spectral')
        plt.colorbar(cf, ax = ax)
        ax.set_title(grid + ' ' + model)
        ax.set_xlabel('TLAI;none')
        ax.set_ylabel('FSDS;W/m^2')
fig.savefig(os.path.join(path_out, '20220105', level + '_photosyntheis.png'), dpi = 600., 
            bbox_inches = 'tight')
plt.close(fig)
