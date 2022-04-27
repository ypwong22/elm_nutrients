""" Compare the time series of multiple runs on the same graph. """
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
from utils.constants import *


path_out = os.path.join(os.environ['PROJDIR'], 'Phenology_ELM', 'output_elm')


# Time series
for varname, unit in zip(var_list, unit_list):
    print(varname, unit)

    data_collect = {}
    for prefix,model in zip(sims_prefix, sims_names):
        for plot in chamber_list:
            data_list = pd.read_csv(os.path.join(path_out, prefix,
                                                 f'plot{plot}_ts_extract.csv'),
                                    index_col = 0, parse_dates = True, header = [0,1])
            for loc in ['hummock', 'hollow']:
                data_collect[(f'P{plot}', loc, model)] = data_list[(loc, varname)]
    data_collect = pd.DataFrame(data_collect)

    # separately plot the elevated CO2 and ambient CO2 chambers
    for co2 in [0,1]:
        if co2 == 0:
            ch_list = chamber_list[:5]
        else:
            ch_list = chamber_list[5:]

        fig, axes = plt.subplots(5, 2, figsize = (16, 16), sharex = True, sharey = True)
        for i, loc in enumerate(['hummock', 'hollow']):
            for j, plot in enumerate(ch_list):
                ax = axes[j,i]

                data = data_collect[(f'P{plot}', loc)]
                h = [None] * len(sims_names)
                for k, m in enumerate(sims_names):
                    temp = savgol_filter(data_collect[(f'P{plot}', loc, m)].values, 61, 3)
                    h[k], = ax.plot(data_collect.index, temp, '-', lw = 0.5)

                    if i == 0:
                        ax.set_ylabel(f'P{plot} ({chamber_levels[plot][0]}$^o$C, {chamber_levels[plot][1]})')
                        if j == 2:
                            ax.text(-0.4, 0.5, f'{varname} ({unit})', fontsize = 14,
                                    transform = ax.transAxes, rotation = 90)
                    if j == 0:
                        ax.set_title(loc.capitalize())
        ax.legend(h, sims_names, loc = [-1, -1], ncol = 3)
        fig.savefig(os.path.join(path_out, 'ts_multiruns', f'{varname}_{co2}.png'), dpi = 600., 
                    bbox_inches = 'tight')
        plt.close(fig)


# Sensitivity of GPP to the variables
for prefix, model in zip(sims_prefix, sims_names):
    for loc, suffix in it.product(['hummock', 'hollow'], ['_pima', '_lala', '_shrub']):
        fig, axes = plt.subplots(len(chamber_list), 4, figsize = (3 * len(chamber_list), 16))
        for j,ch in enumerate(chamber_list):
            data_list = pd.read_csv(os.path.join(path_out, prefix, f'plot{ch}_ts_extract.csv'),
                                    index_col = 0, parse_dates = True, header = [0,1])
            for i, (varname,unit) in enumerate([('FSDS','W m-2'), 
                                                ('TLAI',''),
                                                ('FROOTC','gC m-2')]):
                ax = axes[j,i]
                if varname == 'FSDS':
                    x = data_list[(loc, varname)].values
                else:
                    x = data_list[(loc, varname + suffix)].values
                y = data_list[(loc, 'GPP' + suffix)].values
                ax.plot(x, y, 'o', markersize = 3)
                r, pval = pearsonr(x, y)
                ax.text(0.05, 0.85, '{:.4f}, p = {:.3f}'.format(r, pval), 
                        transform = ax.transAxes)

                if j == 0:
                    ax.set_title(varname)
                elif j == (len(chamber_list)-1):
                    ax.set_xlabel(unit)
                if i == 0:
                    ax.set_ylabel('GPP (gC m-2)')
                if i == 2:
                    ax.text(-0.5, 0.5, 
                            f'P{ch} ({chamber_levels[ch][0]}$^o$C, {chamber_levels[ch][1]})', rotation = 90, transform = ax.transAxes)

        for j,ch in enumerate(chamber_list):
            data_list = pd.read_csv(os.path.join(path_out, prefix,
                                                 f'plot{ch}_ts_extract.csv'),
                                    index_col = 0, parse_dates = True, header = [0,1])

            ax = axes[j,3]

            x = data_list[(loc, 'TLAI' + suffix)].values
            # y = data_list[(loc, 'FSDS')].values
            y = data_list[(loc, 'FROOTC' + suffix)].values
            z = data_list[(loc, 'GPP' + suffix)].values
            cf = ax.scatter(x, y, c = z, s = 3, cmap = 'Spectral')
            plt.colorbar(cf, ax = ax)
            ax.set_title('GPP (gC m-2)')
            ax.set_xlabel('TLAI')
            #ax.set_ylabel('FSDS (W m-2)')
            ax.set_ylabel('FROOTC (gC m-2)')

        fig.savefig(os.path.join(path_out, 'ts_multiruns', model + '_' + loc + \
                                 '_photosyntheis' + suffix + '.png'), dpi = 600., 
                    bbox_inches = 'tight')
        plt.close(fig)
