"""
Process the small collar measured shrub-moss community NEE
Use May 15 - Oct 15 growing season mean value, with fitting RMSE as the uncertainty level
Unit: gC m-2 day-1
"""
import os
import matplotlib.pyplot as plt
import pandas as pd
import itertools as it
import numpy as np
from utils.nee_flux_partitioning import fit_q10, q10_function


########################################################################
# (1) Collect the small collar measurements into well-formatted 
#     single dataframe
########################################################################

def read_2023():
    ###data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ELM_Nutrients', 'input', 'Data and Guide for SPRUCE.104', 
    ###                                'SPRUCE_GHG_Flux_15_min_2023.csv'), 
    ###                index_col = 0, na_values = -9999, parse_dates=True)
    ###data = data.groupby(['date_time_30_min_CST', 'plot', 'warming_treatment', 'CO2_treatment', 'collar_number']).mean()
    ###data = data.drop(['julian_day','meanCH4_conc','meanCO2_conc','CH4_flux'], axis = 1).reset_index()

    data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ELM_Nutrients', 'input', 'Data and Guide for SPRUCE.104', 
                                    'SPRUCE_GHG_Flux_15_min_2023.csv'), na_values = -9999)
    data = data.drop(['date_time_30_min_CST', 'julian_day', 'meanCO2_conc', 'meanCH4_conc', 'CH4_flux'], axis = 1)
    data = data.groupby(['date_time_CST', 'plot', 'warming_treatment', 'CO2_treatment', 'collar_number']).mean().reset_index()
    data.rename({'date_time_CST': 'Timestamp'}, axis = 1, inplace = True)

    # Reformat the data into consistency with 2022: the collar-specific measurements are 
    # indicated by Collar_1 & Collar_2 suffixes instead of a collar_number column

    # Columns that vary per collar
    collar_varying = ['CO2_flux', 'PAR', 'collar_soil_temp', 'collar_VWC']

    # Columns that are identical across collars (shared)
    shared_cols = [c for c in data.columns if c not in collar_varying + ['collar_number']]

    # Split by collar
    c1 = data[data['collar_number'] == 1].copy()
    c2 = data[data['collar_number'] == 2].copy()

    # Rename varying columns with suffix
    c1 = c1.rename(columns={col: f'{col}_Collar_1' for col in collar_varying})
    c2 = c2.rename(columns={col: f'{col}_Collar_2' for col in collar_varying})

    # Drop collar_number before merging
    c1 = c1.drop(columns=['collar_number'])
    c2 = c2.drop(columns=['collar_number'])[['Timestamp', 'plot', 'warming_treatment', 'CO2_treatment'] + [f'{col}_Collar_2' for col in collar_varying]]

    # Merge on shared identifying columns
    ###merge_keys = ['date_time_30_min_CST', 'plot', 'warming_treatment', 'CO2_treatment']
    merge_keys = ['Timestamp', 'plot', 'warming_treatment', 'CO2_treatment']
    result = c1.merge(c2, on=merge_keys, how='outer')

    result['Timestamp'] = pd.to_datetime(result['Timestamp'])

    return result


def read_2022():
    data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ELM_Nutrients', 'input', 'Data and Guide for SPRUCE.104', 
                                    'SPRUCE_GHG_Flux_30_min_2022.csv'), na_values = -9999)
    data = data.groupby(['date_time_30_min_CST', 'plot', 'warming_treatment', 'CO2_treatment']).mean().reset_index()

    # 2022 data only has single PAR
    data['PAR_Collar_1'] = data['PAR']
    data['PAR_Collar_2'] = data['PAR']

    data['Timestamp'] = pd.to_datetime(data['date_time_30_min_CST'])

    data = data.drop(['PAR','DOY','CH4_flux_Collar_1','CH4_flux_Collar_2','date_time_30_min_CST'], axis = 1)

    return data


# It turns out 2023 data is much more complete than 2022; 2022 only has a few months
# Focus on 2023 then
#data = pd.concat([read_2023(), read_2022()], axis = 0, ignore_index=True).sort_values(by=['warming_treatment','CO2_treatment','Timestamp'])
data = read_2023().sort_values(by=['warming_treatment','CO2_treatment','Timestamp'])


########################################################################
# (2) Examine the data availability across the two years
# day and nighttime average separately
########################################################################
warming_list = data['warming_treatment'].drop_duplicates().sort_values().values
co2_list = data['CO2_treatment'].drop_duplicates().values

fig, axes = plt.subplots(4, 3, figsize = (16, 12), sharex = True, sharey = False)
count = 0
for warming, co2 in it.product(warming_list, co2_list):
    subset = data.query('warming_treatment == @warming and CO2_treatment == @co2').set_index('Timestamp')

    subset_1 = subset['CO2_flux_Collar_1'][(subset['PAR_Collar_1'] < 50) & (subset['CO2_flux_Collar_1'] > 0)]
    subset_2 = subset['CO2_flux_Collar_1'][subset['PAR_Collar_1'] > 50]

    subset_3 = subset['CO2_flux_Collar_2'][(subset['PAR_Collar_2'] < 50) & (subset['CO2_flux_Collar_2'] > 0)]
    subset_4 = subset['CO2_flux_Collar_2'][subset['PAR_Collar_2'] > 50]

    # convert unit from umol m-2 s-1 => gC m-2 day-1
    subset_1 *= (86400 * 12 * 1e-6)
    subset_2 *= (86400 * 12 * 1e-6)
    subset_3 *= (86400 * 12 * 1e-6)
    subset_4 *= (86400 * 12 * 1e-6)

    #
    ax = axes.flat[count]
    ax.plot(subset_1.index, subset_1.values, 'or', label = 'Night 1', markersize = 1)
    ax.plot(subset_2.index, subset_2.values, 'ob', label = 'Day 1', markersize = 1)
    ax.plot(subset_3.index, subset_3.values, 'o', color = '#FFC0CB', label = 'Night 2', markersize = 1)
    ax.plot(subset_4.index, subset_4.values, 'o', color = '#00FFFF', label = 'Day 2', markersize = 1)
    ax.legend()

    ax.set_title(f'{warming} {co2}')

    count = count + 1

    plt.setp(ax.get_xticklabels(), rotation = 90)

for i in range(count, len(axes.flat)):
    axes.flat[i].axis('off')

plt.savefig(os.path.join(os.environ['PROJDIR'], 'ELM_Nutrients', 'output', 'extract',
                            'extract_obs_small_collar', f'nee_rawdata.png'), dpi=300)
plt.close(fig)


########################################################################
# (3) Examine the Q10 of nighttime Reco
########################################################################
warming_list = data['warming_treatment'].drop_duplicates().sort_values().values
co2_list = data['CO2_treatment'].drop_duplicates().values

print('CO2,warming,Q10,RMSE')
fig, axes = plt.subplots(2, 5, figsize = (20, 12), sharex = True, sharey = False)
count = 0
for co2, warming in it.product(co2_list, warming_list):
    subset = data.query('warming_treatment == @warming and CO2_treatment == @co2').set_index('Timestamp')

    # restrict to PAR < 50, following Jon Stelling et al. 2026
    t1 = pd.concat([subset['collar_soil_temp_Collar_1'][subset['PAR_Collar_1'] <= 50], 
                    subset['collar_soil_temp_Collar_2'][subset['PAR_Collar_2'] <= 50]])
    #t1 = pd.concat([subset['TS_ 10__B3'][subset['PAR_Collar_1'] <= 50], 
    #                subset['TS_ 10__B3'][subset['PAR_Collar_2'] <= 50]])
    subset_1 = pd.concat([subset['CO2_flux_Collar_1'][subset['PAR_Collar_1'] <= 50],
                          subset['CO2_flux_Collar_2'][subset['PAR_Collar_2'] <= 50]])

    # further restrict to Reco > 0, following Jon Stelling et al. 2026
    filt = subset_1.values > 0
    t1 = t1.values[filt]
    subset_1 = subset_1.values[filt]

    #
    ax = axes.flat[count]
    ax.plot(t1, subset_1, 'o', label = 'Night Obs.', markersize = 3)

    popt,pcov,names = fit_q10(t1, subset_1, T_ref=15)
    tmin, tmax = np.nanmin(t1), np.nanmax(t1)
    t_grid = np.linspace(tmin, tmax, 50)
    fitted = q10_function(t_grid, *popt, T_ref=15)
    ax.plot(t_grid, fitted, '-b', label = 'fitted')

    subset_predicted = q10_function(t1, *popt, T_ref = 15)
    rmse = np.sqrt(np.nanmean(np.power(subset_1 - subset_predicted, 2)))
    print(f'{co2},{warming},{popt[1]},{rmse}')

    ax.text(0.05, 0.85, f'Q10={popt[1]:.3f}\n base={popt[0]:.3f}\n RMSE={rmse:.3f}', transform = ax.transAxes)

    ax.set_title(f'{warming} {co2}')

    count = count + 1

    plt.setp(ax.get_xticklabels(), rotation = 90)

for i in range(count, len(axes.flat)):
    axes.flat[i].axis('off')

plt.savefig(os.path.join(os.environ['PROJDIR'], 'ELM_Nutrients', 'output', 'extract',
                            'extract_obs_small_collar', f'reco_Q10.png'), dpi=300)
plt.close(fig)

