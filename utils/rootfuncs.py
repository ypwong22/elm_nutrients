import os
from .paths import *
import pandas as pd
import numpy as np
import xarray as xr
from .phenofuncs import _find_solstice
from glob import glob
import itertools as it
import warnings


def get_observation():

    # Convert the g biomass to gC
    cconc = pd.read_csv(os.path.join(path_input, 'FRED3_cleaned.csv'))

    picea_rootC = cconc.loc[(cconc['Plant taxonomy_Accepted genus_TPL'] == 'Picea') & (cconc['Plant Taxonomy_Accepted species_TPL'] != 'engelmannii'), 
                            ['Plant taxonomy_Accepted genus_TPL', 'Plant Taxonomy_Accepted species_TPL', 'Root C content']].dropna(axis = 0)
    picea_rootC = picea_rootC['Root C content'].mean() / 1000 # mg/g

    larch_rootC = cconc.loc[(cconc['Plant taxonomy_Accepted genus_TPL'] == 'Larix') & (cconc['Plant Taxonomy_Accepted species_TPL'] == 'decidua'), 
                            ['Plant taxonomy_Accepted genus_TPL', 'Plant Taxonomy_Accepted species_TPL', 'Root C content']].dropna(axis = 0)
    larch_rootC = larch_rootC['Root C content'].mean() / 1000.

    shrub_rootC = cconc.loc[(cconc['Plant taxonomy_Accepted family_TPL'] == 'Ericaceae'), 
                            ['Plant taxonomy_Accepted genus_TPL', 'Plant Taxonomy_Accepted species_TPL', 'Root C content']].dropna(axis = 0)
    shrub_rootC = shrub_rootC['Root C content'].mean() / 1000.


    # Soren's minirhizotron data
    minirhizotron = pd.read_csv(os.path.join(path_input, 'soren_root_prod_mort_growth_20230316.csv'))
    minirhizotron['end_date'] = pd.DatetimeIndex(minirhizotron['end_date']) 
    minirhizotron['start_date'] = pd.DatetimeIndex(minirhizotron['start_date'])
    minirhizotron = minirhizotron.loc[minirhizotron['pft'] != 'sedge']

    # Take the absolute value of mortality rate
    minirhizotron['m_g_d'] = minirhizotron['m_g_d'].abs()

    # Convert from g biomass to gC
    # pima: 0.36, lala: 0.14
    minirhizotron.loc[minirhizotron['pft'] == 'tree', 'm_g_d'] *= (picea_rootC * 0.36 + larch_rootC * 0.14) / 0.5
    minirhizotron.loc[minirhizotron['pft'] == 'tree', 'g_g_d'] *= (picea_rootC * 0.36 + larch_rootC * 0.14) / 0.5
    minirhizotron.loc[minirhizotron['pft'] == 'shrub', 'm_g_d'] *= shrub_rootC
    minirhizotron.loc[minirhizotron['pft'] == 'shrub', 'g_g_d'] *= shrub_rootC

    # 2014-2017 ingrowth core data
    ingrowth = pd.read_csv(os.path.join(path_input, 'spruce_root_ingrowth_2014_2017_20200408.csv'), header = [0, 1], encoding = 'Windows-1252')
    ingrowth.columns = [f'{i} {j}' if not 'Unnamed' in j else f'{i}' for i, j in zip(ingrowth.columns.get_level_values(0), ingrowth.columns.get_level_values(1))]
    ingrowth['start_date yyyy-mm-dd'] = pd.DatetimeIndex(ingrowth['start_date yyyy-mm-dd'])
    ingrowth['end_date yyyy-mm-dd'] = pd.DatetimeIndex(ingrowth['end_date yyyy-mm-dd'])
    ingrowth['root_biomass g m-2 day-1'] = ingrowth['root_biomass g/m2/season'] / np.array([(d.days + 1) for d in (ingrowth['end_date yyyy-mm-dd'] - ingrowth['start_date yyyy-mm-dd'])])
    ingrowth = ingrowth.loc[ingrowth['pft'] != 'graminoid', :]
    ingrowth.loc[ingrowth['pft'] == 'spruce', 'root_biomass g m-2 day-1'] = ingrowth.loc[ingrowth['pft'] == 'spruce', 'root_biomass g m-2 day-1'].values * picea_rootC
    ingrowth.loc[ingrowth['pft'] == 'shrub', 'root_biomass g m-2 day-1'] = ingrowth.loc[ingrowth['pft'] == 'shrub', 'root_biomass g m-2 day-1'].values * shrub_rootC
    ingrowth.loc[ingrowth['pft'] == 'larch', 'root_biomass g m-2 day-1'] = ingrowth.loc[ingrowth['pft'] == 'larch', 'root_biomass g m-2 day-1'].values * larch_rootC

    return minirhizotron, ingrowth



"""
Give set of parameters
* Calculate fine root growth and death rates using default ELM outputs 
* Use `postproc` to compare the calculated values to observations
* Save the metrics to .csv file
"""
def calculator(root_onset, root_dormant, tsoi, tsoi3, swc, psoi, downreg, zwt, sucsat, watsat, bsw, params): # tvec, dt
    """
     Give set of parameters, calculate fine root growth and death rates using default ELM outputs 
    tvec is DatetimeIndex
    dt is the size of the time step (in seconds)
    The drivers are 1-D numpy arrays indexed by tvec in the first dimension
    """
    # calculate the soil matric potential at wd_thres
    psi_crit = - sucsat * np.power(np.maximum(params['wd_thres'] / watsat, 0.01), -bsw)

    td = np.minimum(params['td_max'], params['td_offset'] + params['td_scale'] * np.power(tsoi - params['td_base'], 2))

    wd_le_base = 0.5 + 1 / np.pi * np.arctan(params['wd_scale'] * np.pi * (np.abs(psoi) - params['wd_base']))
    wd_at_base = 0.5 + 1 / np.pi * np.arctan(params['wd_scale'] * np.pi * (np.abs(psi_crit) - params['wd_base']))
    wd_ge_base = (swc - params['wd_thres']) / (1 - params['wd_thres']) * (1 - wd_at_base) + wd_at_base
    wd = np.where(swc >= params['wd_thres'], wd_ge_base, wd_le_base)

    # mortality rate
    mortality = 1 / (params['froot_long'] * 365 * 86400) * 2 * np.maximum(td, wd)

    # compensatory growth rate
    td_0 = np.minimum(params['td_max'], params['td_offset'] + params['td_scale'] * np.power(params['td_base'], 2))
    m0 = 1 / (params['froot_long'] * 365 * 86400) * 2 * td_0

    growth_compens = np.where(tsoi3 < 0, np.maximum(mortality - m0, 0), 0)

    fnmin = 1 + np.power(downreg, params['downreg_a'])
    fw = np.minimum(1, np.maximum(1 - np.power(zwt, params['zwt_a']) / params['zwt_max'], params['zwt_min']))

    # tvec = [pd.date_range(flist_pft[0].split('.')[-2].split('-')[0] + '-01-01 00:00:00', flist_pft[-1].split('.')[-2].split('-')[0] + '-12-31 23:00:00', freq = '1H')
    # tvec = tvec[~((tvec.month == 2) & (tvec.day == 29))]
    #filt = (tvec.year == 2015) & (tvec.month == 9) & (tvec.day == 6) & (tvec.hour == 23)
    #print(tsoi[filt] + 273.15, td[filt], params['td_base'], params['td_max'], params['td_scale'], params['td_offset'])
    #print(swc[filt], abs(psi_crit), psoi[filt], params['wd_thres'], wd_at_base, wd[filt], sucsat, watsat, bsw)
    #print(fnmin[filt], fw[filt])

    # transfer growth rate
    growth_onset = fnmin * fw * params['f_stor'] / 86400
    growth_onset = np.where(tsoi3 >= 0, growth_onset, 0)

    growth_onset = np.where(root_onset > 0, growth_onset, 0.)
    growth_compens = np.where(root_dormant > 0, growth_compens, 0)

    return growth_onset, growth_compens, mortality


def extractor(casename, params):
    path_files = os.path.join(os.environ['PROJDIR'], 'E3SM', 'output', casename, 'run')

    flist_const = sorted(glob(path_files + '/*.h0.*.nc'))[0]
    flist_col = sorted(glob(path_files + '/*.h1.*.nc'))[:-1]
    flist_pft = sorted(glob(path_files + '/*.h2.*.nc'))[:-1]

    pft_list = [2, 3, 11]
    dt = 3600 # hourly

    if casename.startswith('test_'):
        tvec = pd.date_range(flist_pft[0].split('.')[-2].split('-')[0] + '-01-01 00:00:00', flist_pft[-1].split('.')[-2].split('-')[0] + '-12-31 23:00:00', freq = '1H')
    else:
       tvec = pd.date_range(flist_pft[0].split('.')[-2].split('-')[0] + '-01-01 00:00:00', flist_pft[-1].split('.')[-2].split('-')[0] + '-12-31 00:00:00', freq = '1D')
    tvec = tvec[~((tvec.month == 2) & (tvec.day == 29))]

    template = pd.DataFrame(np.nan, index = tvec, columns = pd.MultiIndex.from_product([['2','3','11'], ['hummock', 'hollow']]))
    growth_onset = template.copy() # 1/s
    growth_compens = template.copy() # 1/s
    mortality = template.copy() # 1/s
    root_growth = template.copy() # gC m-2 s-1
    frootc_to_litter = template.copy() # gC m-2 s-1

    hr = xr.open_mfdataset(flist_const)
    sucsat = pd.Series(hr['SUCSAT'][2, :].values, index = ['hummock', 'hollow'])
    watsat = pd.Series(hr['WATSAT'][2, :].values, index = ['hummock', 'hollow'])
    bsw = pd.Series(hr['BSW'][2, :].values, index = ['hummock', 'hollow'])
    dzsoi = hr['DZSOI'].values
    hr.close()

    hr = xr.open_mfdataset(flist_col)
    hr2 = xr.open_mfdataset(flist_pft)
    for pft, loc in it.product(pft_list, ['hummock','hollow']):
        if loc == 'hummock':
            col_id = 0
            pft_id = pft
        else:
            col_id = 1
            pft_id = pft + 17
        
        if pft == 2:
            params['froot_long'] = 3.5 * params['froot_long']
        elif pft == 3:
            params['froot_long'] = 1 * params['froot_long']
        else:
            params['froot_long'] = 1 * params['froot_long']


        if casename.startswith('test_'):
            root_onset = (hr2['ONSET_FLAG_ROOT'][:, pft_id].values > 0)
            root_dormant = (hr2['DORMANT_FLAG_ROOT'][:, pft_id].values > 0)
        else:
            # flag: between winter solstice and date of onset, or not
            leaf_onset = hr2['ONSET_FLAG'][:, pft_id].values
            root_onset = np.zeros(len(leaf_onset))
            root_dormant = np.ones(len(leaf_onset))
            for yy in np.unique(tvec.year):
                if np.sum((tvec.year == yy) & (leaf_onset > 0)) > 0:
                    day0 = np.where((tvec.year == yy) & (leaf_onset > 0))[0][0]
                    day_sol = _find_solstice([yy]).values[0].dayofyear # np.where((tvec.year == yy) & ((tvec.month * 100 + tvec.day) >= 1221))[0][0]
                    root_onset[day0:day_sol] = 1
                    root_dormant[day0:day_sol] = 0

        rootfr = hr2['ROOTFR'][:, :, pft_id].values

        tsoi = np.sum((hr['TSOI'][:, :, col_id].values - 273.15) * rootfr, axis = 1)
        tsoi3 = hr['TSOI'][:, 2, col_id].values - 273.15
        swc = np.sum(hr['SOILLIQ'][:, :, col_id].values / dzsoi[:, col_id].reshape(1, -1) / 1000 * rootfr, axis = 1)
        psoi = np.nanmax(hr['SMP'][:, :, col_id].values, axis = 1)  # mm; need to convert to MPa
        zwt = hr['H2OSFC'][:, col_id].values

        downreg = hr2['DOWNREG'][:, pft_id].values

        growth_onset.loc[:, (str(pft), loc)], growth_compens.loc[:, (str(pft), loc)], mortality.loc[:, (str(pft), loc)] = \
            calculator(root_onset, root_dormant, tsoi, tsoi3, swc, psoi, downreg, zwt, sucsat.loc[loc], watsat.loc[loc], bsw.loc[loc], params)

        # Convert onset rate & compense rate & mortality rate to growth and litterfall rate gC m-2 s-1 -> gC m-2 day-1
        root_growth.loc[:, (str(pft), loc)] = 86400 * np.minimum(growth_onset.loc[:, (str(pft), loc)] + growth_compens.loc[:, (str(pft), loc)], 1 / dt) * hr2['FROOTC_STORAGE'][:, pft_id].values
        frootc_to_litter.loc[:, (str(pft), loc)] = 86400 *mortality.loc[:, (str(pft), loc)] * hr2['FROOTC'][:, pft_id].values

    hr.close()
    hr2.close()
    return growth_onset, growth_compens, mortality, root_growth, frootc_to_litter


""" Use the default ELM version, convert to the same dates as Soren's observations """
def convert(params, minirhizotron, ingrowth):
    tvec = pd.date_range('2015-01-01 00:00:00', '2020-12-31 00:00:00', freq = '1D')
    tvec = tvec[~((tvec.month == 2) & (tvec.day == 29))]
    root_growth_collect = pd.DataFrame(index = tvec,
                                    columns = pd.MultiIndex.from_product([['hummock', 'hollow'], chamber_list, ['2', '3', '11']]))
    frootc_to_litter_collect = pd.DataFrame(index = tvec, 
                                            columns = pd.MultiIndex.from_product([['hummock', 'hollow'], chamber_list, ['2', '3', '11']]))
    warnings.filterwarnings('ignore')
    for plot in chamber_list:
        _, _, _, root_growth, frootc_to_litter = extractor(f'20221212_plot{plot:02g}_US-SPR_ICB20TRCNPRDCTCBC', params)
        for col in root_growth.columns:
            root_growth_collect.loc[:, (col[1], plot, col[0])] = root_growth[col].values
            frootc_to_litter_collect.loc[:, (col[1], plot, col[0])] = frootc_to_litter[col].values
    warnings.filterwarnings('default')
    root_growth_collect = root_growth_collect['hummock'] * 0.64 + root_growth_collect['hollow'] * 0.36
    frootc_to_litter_collect = frootc_to_litter_collect['hummock'] * 0.64 + frootc_to_litter_collect['hollow'] * 0.36


    def _expand(x):
        data = 0.64 * x.loc[x.index.get_level_values('topo') == 'hummock', :].values + \
            0.36 * x.loc[x.index.get_level_values('topo') == 'hollow', :].values
        return pd.Series(data.reshape(-1), index = x.columns)


    minirhizotron_average = minirhizotron.drop(['temp', 'co2', 'year', 'time_step', 'npp_km_d', 'm_km_d', 'g_km_d', 'npp_g_d', 'swc', 'tube'], 
                                            axis = 1).set_index(['plot', 'topo', 'pft', 'start_date', 'end_date'])
    idx = pd.MultiIndex.from_product([minirhizotron_average.index.levels[0], minirhizotron_average.index.levels[1], minirhizotron_average.index.levels[2], 
                                    minirhizotron_average.index.levels[3], minirhizotron_average.index.levels[4]], names = ['plot', 'topo', 'pft', 'start_date', 'end_date'])
    minirhizotron_average = minirhizotron_average.reindex(idx).groupby(['plot', 'pft', 'start_date', 'end_date']).apply(_expand)
    minirhizotron_average = minirhizotron_average.dropna(axis = 0, how = 'all')


    ingrowth2 = ingrowth.loc[season == 'summer', :].copy()
    ingrowth2 = ingrowth2.groupby(['topog hummock or hollow', 'plot', 'start_date yyyy-mm-dd', 'end_date yyyy-mm-dd', 'pft']).sum()[['root_biomass g m-2 day-1']]
    # hummock-hollow averages
    ingrowth2 = (ingrowth2.loc['hollow', :] * 0.64 + ingrowth2.loc['hummock', :] * 0.36).dropna(axis = 0)
    ingrowth2 = ingrowth2.iloc[:, 0].unstack()

    minirhizotron_sim = pd.DataFrame(np.nan, index = minirhizotron_average.index, columns = minirhizotron_average.columns)
    for i, _ in minirhizotron_average.iterrows():
        if i[0] in [7, 21]:
            continue # temporarily skip TAMB
        if (i[2].year < 2015) | (i[3].year > 2020):
            continue
        if i[1] == 'tree':
            minirhizotron_sim.loc[i, 'g_g_d'] = root_growth_collect.loc[i[2]:i[3], (i[0], '2')].mean() * 0.36 + root_growth_collect.loc[i[2]:i[3], (i[0], '3')].mean() * 0.14
            minirhizotron_sim.loc[i, 'm_g_d'] = frootc_to_litter_collect.loc[i[2]:i[3], (i[0], '2')].mean() * 0.36 + frootc_to_litter_collect.loc[i[2]:i[3], (i[0], '3')].mean() * 0.14
        elif i[1] == 'shrub':
            minirhizotron_sim.loc[i, 'g_g_d'] = root_growth_collect.loc[i[2]:i[3], (i[0], '11')].mean() * 0.25
            minirhizotron_sim.loc[i, 'm_g_d'] = frootc_to_litter_collect.loc[i[2]:i[3], (i[0], '11')].mean() * 0.25

    ingrowth_sim = pd.DataFrame(np.nan, index = ingrowth2.index, columns = ingrowth2.columns)
    for i, _ in ingrowth2.iterrows():
        if i[0] in [7, 21]:
            continue # temporarily skip TAMB
        if (i[1].year < 2015) | (i[2].year > 2020):
            continue
        ingrowth_sim.loc[i, 'spruce'] = root_growth_collect.loc[i[1]:i[2], (i[0], '2')].mean()
        ingrowth_sim.loc[i, 'larch'] = root_growth_collect.loc[i[1]:i[2], (i[0], '3')].mean()
        ingrowth_sim.loc[i, 'shrub'] = root_growth_collect.loc[i[1]:i[2], (i[0], '11')].mean()

    return minirhizotron_sim, ingrowth_sim, minirhizotron_average, ingrowth2
