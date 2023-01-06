""" Plot against the data in spruce_validation_data.nc """
import itertools as it
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import numpy as np
import os
from glob import glob
from datetime import datetime, timedelta
from utils.constants import *
from optparse import OptionParser
from scipy.stats import pearsonr


"""
parser = OptionParser()
parser.add_option("--sim_id", dest="sim_id", default = '', help = 'Index to the simulation to extract.')
parser.add_option('--plot_id', dest = 'plot_id', default = '', help = 'Index to the chamber to extract.')
(options,args) = parser.parse_args()

prefix  = sims_prefix[int( options.sim_id)]
plot    = chamber_list[int(options.plot_id)]
print(prefix + '_plot' + plot)
"""


prefix = '20221231'
for plot in [4, 10, 11, 16, 19, 6, 8, 13, 17, 20]:
    #############################################################
    # Some constants
    #############################################################
    # weights for the grid
    gridweights = pd.Series([0.64, 0.36], index = ['hummock', 'hollow'])


    # moss fractions
    def get_mossfrac():
        mossfrac = pd.read_excel('Sphagnum_fraction.xlsx', index_col = 0, skiprows = 1,
                                engine = 'openpyxl').drop(['plot','Temp','CO2'], axis = 1)
        mossfrac = mossfrac / 100. # percentage -> fraction
        mossfrac[2015] = mossfrac[2016]
        grid_to_plot = {'T0.00': 'P06', 'T2.25': 'P20',
                        'T4.50': 'P13', 'T6.75': 'P08',
                        'T9.00': 'P17', 'T0.00CO2': 'P19',
                        'T2.25CO2': 'P11', 'T4.50CO2': 'P04',
                        'T6.75CO2': 'P16', 'T9.00CO2': 'P10',
                        'TAMB'    : 'P07'}
        mossfrac.index = [int(grid_to_plot[i].replace('P','')) for i in mossfrac.index]
        return mossfrac.T.sort_index()
    mossfrac       = get_mossfrac()

    mossfrac_daily = pd.DataFrame(np.nan, index = pd.date_range('2015-01-01', '2021-12-31'),
                                columns = mossfrac.columns)
    for yy in range(2015,2022):
        mossfrac_daily.loc[mossfrac_daily.index.year == yy, :] = mossfrac.loc[yy, :].values


    # open the validation data file
    folder = os.path.join(os.path.join(os.environ['PROJDIR'], 'E3SM', 'output', f'{prefix}_plot{plot:02d}_US-SPR_ICB20TRCNPRDCTCBC', 'run'))
    validation_file = os.path.join(os.environ['PROJDIR'], 'Phenology_ELM', 'intermediate', 'spruce_validation_data.nc')
    hrv = xr.open_dataset(validation_file)


    ###########################################################
    # CH4 and CO2 fluxes
    ###########################################################
    def get_run_co2():
        # gC m-2 day-1, whole grid

        data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'Phenology_ELM', 'output_elm', prefix, f'plot{plot}_ts_extract.csv'),
                        index_col = 0, header = [0,1])

        rh = data.loc[:, (['hummock','hollow'], 'HR')]
        rh.columns = rh.columns.droplevel(1)
        rh = (rh * gridweights).sum(axis = 1)

        # NPP; use only shrub and sedge PFTs
        # shrub: 0.25, sedge: depends on the year
        npp = (data.loc[:,('hummock','NPP_shrub')] * 0.64 + \
            data.loc[:,('hollow' ,'NPP_shrub')] * 0.36) * 0.25 + \
            (data.loc[:,('hummock','NPP_moss' )] * 0.64 + \
            data.loc[:,('hollow' ,'NPP_moss' )] * 0.36) * \
            mossfrac_daily.loc[sim_tvec, int(plot)].values

        return pd.Series(npp.values - rh.values, index = sim_tvec)


    def get_run_ch4():
        data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'Phenology_ELM', 'output_elm', prefix, f'plot{plot}_ts_extract.csv'), index_col = 0, header = [0,1])

        # gC m-2 day-1
        ch4 = data.loc[:, (['hummock','hollow'], 'CH4PROD')]
        ch4.columns = ch4.columns.droplevel(1)
        ch4 = (ch4 * gridweights).sum(axis = 1)

        return pd.Series(ch4.values, sim_tvec)


    sim_co2 = get_run_co2()
    # sim_ch4 = get_run_ch4()

    obs_co2 = pd.Series(hrv['daily_nee'].loc[:, int(plot), :].mean(dim = 'isday').values, 
                        index = hrv['time'].to_index()) * 86400 # s-1 -> day-1
    # obs_ch4 = pd.Series(hrv['daily_ch4'].loc[:, int(plot), :].mean(dim = 'isday').values, 
    #                     index = hrv['time'].to_index()) * 86400 # s-1 -> day-1

    co2_rmse  = np.sqrt(np.nanmean(np.power(sim_co2 - obs_co2 , 2)))
    # ch4_rmse  = np.sqrt(np.nanmean(np.power(sim_ch4 - obs_ch4, 2)))

    temp      = pd.concat([sim_co2, obs_co2], axis = 1, join = 'inner').dropna(axis = 0, how = 'any')
    co2_rho   = np.corrcoef(temp.values.T)[0,1]
    co2_bias  = temp.values[:,0].mean() - temp.values[:,1].mean() # sim - obs
    # temp      = pd.concat([sim_ch4, obs_ch4], axis = 1, join = 'inner').dropna(axis = 0, how = 'any')
    # ch4_rho   = np.corrcoef(temp.values.T)[0,1]
    # ch4_bias  = temp.values[:,0].mean() - temp.values[:,1].mean()


    # fig, axes = plt.subplots(2, 1, figsize = (12, 12))
    # ax = axes.flat[0]
    fig, ax = plt.subplots(figsize = (6, 6))
    h, = ax.plot(sim_co2.index, sim_co2, '-b', lw = 0.5)
    ax.plot(obs_co2.index, obs_co2.values, 'ok', markersize = 8)
    h2 = ax.errorbar(obs_co2.index, obs_co2.values,
                    yerr = hrv['daily_nee_se'].loc[:,int(plot),:].mean(dim = 'isday').values * 86400,
                    color = 'k', ecolor = 'k', elinewidth = 1)
    ax.legend([h, h2], ['Model', 'Obs±SD'])
    ax.set_title('CO$_2$')
    ax.set_ylabel('gC m-2 day-1')
    ax.text(0.05, 0.8, f'RMSE = {co2_rmse:.04e}', transform = ax.transAxes)
    """
    ax = axes.flat[1]
    h, = ax.plot(sim_ch4.index, sim_ch4, '-b', lw = 0.5)
    ax.plot(obs_ch4.index, obs_ch4.values, 'ok', markersize = 8)
    h2 = ax.errorbar(obs_ch4.index, obs_ch4.values,
                    yerr = hrv['daily_ch4_se'].loc[:, int(plot),
                                                    :].mean(dim = 'isday').values * 86400,
                    color = 'k', ecolor = 'k', elinewidth = 1)
    ax.set_title('CH$_4$')
    ax.set_ylabel('gC m-2 day-1')
    ax.text(0.05, 0.8, f'RMSE = {ch4_rmse:.04e}', transform = ax.transAxes)
    """
    fig.savefig(os.path.join(os.environ['PROJDIR'], 'Phenology_ELM', 'output_elm', prefix, f'plot{plot:02d}_daily_nee.png'), dpi = 600.,
                bbox_inches = 'tight')
    plt.close()


    ###########################################################
    # Annual NPP
    ###########################################################
    def get_anpp_tree():
        data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'Phenology_ELM', 'output_elm', prefix, f'plot{plot}_ts_extract.csv'),
                        index_col = 0, header = [0,1], parse_dates = True)
        # gC m-2 day-1 -> yr-1
        npp = (data.loc[:,('hummock','NPP_pima')] * 0.64 + \
            data.loc[:,('hollow' ,'NPP_pima')] * 0.36) * 0.36 + \
            (data.loc[:,('hummock','NPP_lala' )] * 0.64 + \
            data.loc[:,('hollow' ,'NPP_lala' )] * 0.36) * 0.14
        return npp.groupby(npp.index.year).sum()

    def get_anpp_shrub():
        data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'Phenology_ELM', 'output_elm', prefix, f'plot{plot}_ts_extract.csv'),
                        index_col = 0, header = [0,1], parse_dates = True)
        # gC m-2 day-1 -> yr-1
        npp = (data.loc[:,('hummock','NPP_shrub')] * 0.64 + \
            data.loc[:,('hollow' ,'NPP_shrub')] * 0.36) * 0.25
        return npp.groupby(npp.index.year).sum()

    def get_npp_moss():
        data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'Phenology_ELM', 'output_elm', prefix, f'plot{plot}_ts_extract.csv'),
                        index_col = 0, header = [0,1], parse_dates = True)
        # gC m-2 day-1 -> yr-1
        npp = (data.loc[:,('hummock','NPP_moss')] * 0.64 + \
            data.loc[:,('hollow' ,'NPP_moss')] * 0.36) * \
            mossfrac_daily.loc[sim_tvec, int(plot)]
        return npp.groupby(npp.index.year).sum()

    def get_bnpp():
        data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'Phenology_ELM', 'output_elm', prefix, f'plot{plot}_ts_extract.csv'),
                        index_col = 0, header = [0,1], parse_dates = True)
        # gC m-2 day-1 -> yr-1
        npp = data.loc[:,('hummock', 'BGNPP')] * 0.64 + data.loc[:, ('hollow' ,'BGNPP')] * 0.36
        return npp.groupby(npp.index.year).sum()


    sim_anpp_tree    = get_anpp_tree()
    sim_anpp_shrub   = get_anpp_shrub()
    sim_npp_moss     = get_npp_moss()
    sim_bnpp         = get_bnpp()

    obs_anpp_tree    = pd.Series(hrv['annual_anpp_tree'].loc[:, int(plot)].values, 
                                index = hrv['year'].to_index())
    obs_anpp_shrub   = pd.Series(hrv['annual_anpp_shrub'].loc[:, int(plot)].values, 
                                index = hrv['year'].to_index())
    obs_npp_moss     = pd.Series(hrv['annual_npp_moss'].loc[:, int(plot)].values, 
                                index = hrv['year'].to_index())
    obs_bnpp         = pd.Series(hrv['annual_bnpp'].loc[:, int(plot)].values, 
                                index = hrv['year'].to_index())

    anpp_tree_rmse     = np.sqrt(np.nanmean(np.power(sim_anpp_tree  - obs_anpp_tree , 2)))
    anpp_shrub_rmse    = np.sqrt(np.nanmean(np.power(sim_anpp_shrub - obs_anpp_shrub, 2)))
    npp_moss_rmse      = np.sqrt(np.nanmean(np.power(sim_npp_moss   - obs_npp_moss  , 2)))
    bnpp_rmse          = np.sqrt(np.nanmean(np.power(sim_bnpp       - obs_bnpp      , 2)))

    temp             = pd.concat([sim_anpp_tree, obs_anpp_tree], axis = 1,
                                join = 'inner').dropna(axis = 0, how = 'any')
    anpp_tree_rho    = np.corrcoef(temp.values.T)[0,1]
    anpp_tree_bias   = temp.values[:,0].mean() - temp.values[:,1].mean()
    temp             = pd.concat([sim_anpp_shrub, obs_anpp_shrub], axis = 1,
                                join = 'inner').dropna(axis = 0, how = 'any')
    anpp_shrub_rho   = np.corrcoef(temp.values.T)[0,1]
    anpp_shrub_bias  = temp.values[:,0].mean() - temp.values[:,1].mean()
    temp             = pd.concat([sim_npp_moss, obs_npp_moss], axis = 1,
                                join = 'inner').dropna(axis = 0, how = 'any')
    npp_moss_rho     = np.corrcoef(temp.values.T)[0,1]
    npp_moss_bias    = temp.values[:,0].mean() - temp.values[:,1].mean()
    temp             = pd.concat([sim_bnpp, obs_bnpp], axis = 1,
                                join = 'inner').dropna(axis = 0, how = 'any')
    bnpp_rho         = np.corrcoef(temp.values.T)[0,1]
    bnpp_bias        = temp.values[:,0].mean() - temp.values[:,1].mean()


    fig, axes = plt.subplots(2, 2, figsize = (12, 12))
    ax = axes.flat[0]
    h, = ax.plot(sim_anpp_tree.index, sim_anpp_tree, '-b', lw = 0.5)
    h2, = ax.plot(obs_anpp_tree.index, obs_anpp_tree.values, 'ok', markersize = 8)
    ax.legend([h, h2], ['Model', 'Obs'])
    ax.set_title('ANPP$_{tree}$')
    ax.set_ylabel('gC m-2 yr-1')
    ax.text(0.05, 0.8, f'RMSE = {anpp_tree_rmse:.04e}', transform = ax.transAxes)
    ax = axes.flat[1]
    h, = ax.plot(sim_anpp_shrub.index, sim_anpp_shrub, '-b')
    h2, = ax.plot(obs_anpp_shrub.index, obs_anpp_shrub.values, 'ok', markersize = 8)
    ax.set_title('ANPP$_{shrub}$')
    ax.set_ylabel('gC m-2 yr-1')
    ax.text(0.05, 0.8, f'RMSE = {anpp_shrub_rmse:.04e}', transform = ax.transAxes)
    ax = axes.flat[2]
    h, = ax.plot(sim_npp_moss.index, sim_npp_moss, '-b')
    h2, = ax.plot(obs_npp_moss.index, obs_npp_moss.values, 'ok', markersize = 8)
    ax.set_title('NPP$_{moss}$')
    ax.set_ylabel('gC m-2 yr-1')
    ax.text(0.05, 0.8, f'RMSE = {npp_moss_rmse:.04e}', transform = ax.transAxes)
    ax = axes.flat[3]
    h, = ax.plot(sim_bnpp.index, sim_bnpp, '-b')
    h2, = ax.plot(obs_bnpp.index, obs_bnpp.values, 
                'ok', markersize = 8)
    ax.set_title('BNPP')
    ax.set_ylabel('gC m-2 yr-1')
    ax.text(0.05, 0.8, f'RMSE = {bnpp_rmse:.04e}', transform = ax.transAxes)
    fig.savefig(os.path.join(os.environ['PROJDIR'], 'Phenology_ELM', 'output_elm', prefix, f'plot{plot:02d}_annual_npp.png'), dpi = 600.,
                bbox_inches = 'tight')
    plt.close()


    ###########################################################
    # Daily outflow (runoff)
    ###########################################################
    def get_runoff():
        data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'Phenology_ELM', 'output_elm', prefix, f'plot{plot}_ts_extract.csv'),
                        index_col = 0, header = [0,1])
        # mm/s -> mm/day
        runoff = data.loc[:, (['hummock','hollow'], 'QOVER')] * 86400
        runoff.columns = runoff.columns.droplevel(1)
        runoff = (runoff * gridweights).sum(axis = 1)
        return pd.Series(runoff.values, index = sim_tvec)


    sim_runoff   = get_runoff()
    obs_runoff   = pd.Series(hrv['daily_runoff'].loc[:, int(plot)].values,
                            index = hrv['time'].to_index()) * 86400 # mm/s -> mm/day

    runoff_rmse  = np.sqrt(np.nanmean(np.power(sim_runoff - obs_runoff, 2)))
    temp         = pd.concat([sim_runoff, obs_runoff], axis = 1,
                            join = 'inner').dropna(axis = 0, how = 'any')
    runoff_rho   = np.corrcoef(temp.values.T)[0,1]
    runoff_bias  = temp.values[:,0].mean() - temp.values[:,1].mean()

    fig, ax = plt.subplots(figsize = (12, 6))
    h, = ax.plot(sim_runoff.index, np.log(sim_runoff + 1), '-b', lw = 0.5)
    h2, = ax.plot(obs_runoff.index, np.log(obs_runoff.values + 1), 'ok', markersize = 4)
    ax.semilogy()
    ax.set_ylabel('log(1 + mm day-1)')
    ax.text(0.05, 0.9, f'RMSE = {runoff_rmse:.04e} (mm day-1)', transform = ax.transAxes)
    fig.savefig(os.path.join(os.environ['PROJDIR'], 'Phenology_ELM', 'output_elm', prefix, f'plot{plot:02d}_daily_runoff.png'), dpi = 600.,
                bbox_inches = 'tight')
    plt.close()

    ###########################################################
    # Water table (zwt)
    ###########################################################
    def get_zwt():
        data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'Phenology_ELM', 'output_elm', prefix, f'plot{plot}_ts_extract.csv'),
                        index_col = 0, header = [0,1])
        # hollow and hummock already adjusted
        zwt = data.loc[:, (['hummock','hollow'], 'ZWT')]
        zwt.columns = zwt.columns.droplevel(1)
        zwt = (zwt * gridweights).sum(axis = 1)
        return pd.Series(zwt.values, index = sim_tvec)


    sim_zwt = get_zwt()
    obs_zwt = pd.Series(hrv['daily_zwt'].loc[:, int(plot)].values,
                        index = hrv['time'].to_index())

    zwt_rmse  = np.sqrt(np.nanmean(np.power(sim_zwt - obs_zwt, 2)))
    temp      = pd.concat([sim_zwt, obs_zwt], axis = 1,
                        join = 'inner').dropna(axis = 0, how = 'any')
    zwt_rho   = np.corrcoef(temp.values.T)[0,1]
    zwt_bias  = temp.values[:,0].mean() - temp.values[:,1].mean()

    fig, ax = plt.subplots(figsize = (12, 6))
    h, = ax.plot(sim_zwt.index, sim_zwt, '-b', lw = 0.5)
    h2, = ax.plot(obs_zwt.index, obs_zwt, 'ok', markersize = 4)
    ax.set_ylabel('m (relative to hollow surface)')
    ax.text(0.8, 0.9, f'RMSE = {zwt_rmse:.04e} (m)', transform = ax.transAxes)
    fig.savefig(os.path.join(os.environ['PROJDIR'], 'Phenology_ELM', 'output_elm', prefix, f'plot{plot:02d}_daily_watertable.png'), dpi = 600.,
                bbox_inches = 'tight')
    plt.close()

    ###########################################################
    # Collect the rmse
    ###########################################################
    RMSE = pd.DataFrame({
        'RMSE': pd.Series({'co2': co2_rmse,  # 'ch4': ch4_rmse,
                        'anpp_tree': anpp_tree_rmse,
                        'anpp_shrub': anpp_shrub_rmse,
                        'npp_moss': npp_moss_rmse,
                        'bnpp': bnpp_rmse,
                        'runoff': runoff_rmse,
                        'zwt': zwt_rmse}),
        'Corr': pd.Series({'co2': co2_rho, # 'ch4': ch4_rho,
                        'anpp_tree': anpp_tree_rho,
                        'anpp_shrub': anpp_shrub_rho,
                        'npp_moss': npp_moss_rho,
                        'bnpp': bnpp_rho,
                        'runoff': runoff_rho,
                        'zwt': zwt_rho}),
        'Bias': pd.Series({'co2': co2_bias, # 'ch4': ch4_bias,
                        'anpp_tree': anpp_tree_bias,
                        'anpp_shrub': anpp_shrub_bias,
                        'npp_moss': npp_moss_bias,
                        'bnpp': bnpp_bias,
                        'runoff': runoff_bias,
                        'zwt': zwt_bias})
        })
    RMSE.to_csv(os.path.join(os.environ['PROJDIR'], 'Phenology_ELM', 'output_elm', prefix, f'plot{plot:02d}_RMSE.csv'))

    ###########################################################
    hrv.close()
