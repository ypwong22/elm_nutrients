import os
import xarray as xr
import pandas as pd
import numpy as np
from glob import glob
from utils.constants import chamber_list_complete_dict
from utils.analysis import get_sim_carbonfluxes


grid_to_plot = {
    "T0.00": "P06",
    "T2.25": "P20",
    "T4.50": "P13",
    "T6.75": "P08",
    "T9.00": "P17",
    "T0.00CO2": "P19",
    "T2.25CO2": "P11",
    "T4.50CO2": "P04",
    "T6.75CO2": "P16",
    "T9.00CO2": "P10",
    "TAMB": "P07",
}
plot_to_grid = dict([(b,a) for a,b in grid_to_plot.items()])


pft_stride = 17
plot_list = ['P04', 'P06', 'P07', 'P08', 'P10', 'P11', 'P13', 'P16', 'P17', 'P19', 'P20']

var_list_extra = ['ZWT', 'TOTSOMC', 'SMINN_30', 'SOLUTIONP_30',
                  'FPG_Spruce', 'FPG_Tamarack', 'FPG_Shrub',
                  'FPG_P_Spruce', 'FPG_P_Tamarack', 'FPG_P_Shrub', 
                  'FPI', 'FPI_P', 'SOMN_TO_SMINN_Spruce', 
                  'SOMN_TO_SMINN_Tamarack', 'SOMN_TO_SMINN_Shrub',
                  'SOMP_TO_SMINP_Spruce', 'SOMP_TO_SMINP_Tamarack',
                  'SOMP_TO_SMINP_Shrub']


prefix = "20240312"
#prefix  = "UQ_20240311_4" # UQ_20240311_3 = g03505, UQ_20240311_4 = g01540
extrafix = "" # "_alt_params"
growing_season = False
zwt_growing_season = True


year_range = range(2015, 2022)
if not "UQ" in prefix:
    runroot = os.path.join(os.environ['E3SM_ROOT'], 'output', 
                           f'{prefix}_US-SPR_ICB20TRCNPRDCTCBC', f'spruce_treatments{extrafix}')
else:
    runroot = os.path.join(os.environ['E3SM_ROOT'], 'output', 'UQ',
                           f'{prefix}_US-SPR_ICB20TRCNPRDCTCBC', 'g01438')

collect_part1 = get_sim_carbonfluxes(year_range, runroot, growing_season, 
                                     extra_pft_vars = ['LEAFC_ALLOC_TO_TOTVEGC_ABG'])

collect_part2 = pd.DataFrame(
    np.nan,
    index = pd.MultiIndex.from_product([['hummock', 'hollow', 'average'],
                                        plot_list,
                                        year_range], names = ['column', 'plot', 'year']),
    columns = var_list_extra)


for plot in plot_list:
    if not "UQ" in prefix:
        temp = plot.replace("P", '')
        rundir = os.path.join(runroot, f'plot{temp}_US-SPR_ICB20TRCNPRDCTCBC', 'run')
    else:
        temp = chamber_list_complete_dict[plot]
        rundir = os.path.join(runroot, temp)

    flist = sorted(glob(rundir + "/*.h2.*.nc"))
    flist = [f for f in flist if \
             int(f.split('/')[-1].split('.')[-2].split('-')[0]) in year_range]
    # if "UQ" in prefix:
    #     flist = flist[:-1]
    # if len(flist) == 0:
    #    continue
    hr = xr.open_mfdataset(flist)

    flist = sorted(glob(rundir + "/*.h1.*.nc"))
    flist = [f for f in flist if \
             int(f.split('/')[-1].split('.')[-2].split('-')[0]) in year_range]
    #if "UQ" in prefix:
    #    flist = flist[:-1]
    hr2 = xr.open_mfdataset(flist)

    filter = hr['time'].to_index().month *100 + hr['time'].to_index().day
    filter = (filter >= 515) & (filter <= 1015)

    if 'FPI' in var_list_extra:
        # NP limitation should always focus on growing season
        if 'FPG_PATCH' in hr.variables:
            # the FPG_PATCH variable can have benign extremely high values
            # because the current step's demand is too small compared to
            # available fie root uptake
            pot_smin_demand = hr['PLANT_NDEMAND_POT'][filter, :].resample({'time': '1Y'}).mean()
            smin_demand = hr['PLANT_NDEMAND'][filter, :].resample({'time': '1Y'}).mean()
            # som = hr['FUNGI_SOM_TO_NPOOL'][filter, :].resample({'time': '1Y'}).mean()
            for col, add in zip(['hummock', 'hollow'], [0, pft_stride]):
                collect_part2.loc[(col, plot), 'FPG_Spruce'] = (pot_smin_demand / smin_demand)[:, 2 + add]
                collect_part2.loc[(col, plot), 'FPG_Tamarack'] = (pot_smin_demand / smin_demand)[:, 3 + add]
                collect_part2.loc[(col, plot), 'FPG_Shrub'] = (pot_smin_demand / smin_demand)[:, 11 + add]
        else:
            temp = hr2['FPG'][filter, :].resample({'time': '1Y'}).mean()
            for num, col in enumerate(['hummock', 'hollow']):
                collect_part2.loc[(col, plot), 'FPG_Spruce'] = temp[:, num]
                collect_part2.loc[(col, plot), 'FPG_Tamarack'] = temp[:, num]
                collect_part2.loc[(col, plot), 'FPG_Shrub'] = temp[:, num]

        if 'FPG_P_PATCH' in hr.variables:
            pot_smin_demand = hr['PLANT_PDEMAND_POT'][filter, :].resample({'time': '1Y'}).mean()
            smin_demand = hr['PLANT_PDEMAND'][filter, :].resample({'time': '1Y'}).mean()
            som = hr['FUNGI_SOM_TO_PPOOL'][filter, :].resample({'time': '1Y'}).mean()
            for col, add in zip(['hummock', 'hollow'], [0, pft_stride]):
                collect_part2.loc[(col, plot), 'FPG_P_Spruce'] = (pot_smin_demand / smin_demand)[:, 2 + add]
                collect_part2.loc[(col, plot), 'FPG_P_Tamarack'] = (pot_smin_demand / smin_demand)[:, 3 + add]
                collect_part2.loc[(col, plot), 'FPG_P_Shrub'] = (pot_smin_demand / smin_demand)[:, 11 + add]
        else:
            temp = hr2['FPG_P'][filter, :].resample({'time': '1Y'}).mean()
            for num, col in enumerate(['hummock', 'hollow']):
                collect_part2.loc[(col, plot), 'FPG_P_Spruce'] = temp[:, num]
                collect_part2.loc[(col, plot), 'FPG_P_Tamarack'] = temp[:, num]
                collect_part2.loc[(col, plot), 'FPG_P_Shrub'] = temp[:, num]
    
    if 'FPG_PATCH' in hr.variables:
        # calculate the ratio of organic nutrient uptake to inorganic nutrient uptake
        smin_demand = hr['PLANT_NDEMAND'][filter, :].resample({'time': '1Y'}).mean()
        som = hr['FUNGI_SOM_TO_NPOOL'][filter, :].resample({'time': '1Y'}).mean()
        for col, add in zip(['hummock', 'hollow'], [0, pft_stride]):
            collect_part2.loc[(col, plot), 'SOMN_TO_SMINN_Spruce'] = (som / smin_demand)[:, 2 + add]
            collect_part2.loc[(col, plot), 'SOMN_TO_SMINN_Tamarack'] = (som / smin_demand)[:, 3 + add]
            collect_part2.loc[(col, plot), 'SOMN_TO_SMINN_Shrub'] = (som / smin_demand)[:, 11 + add]

    if 'FPG_P_PATCH' in hr.variables:
        smin_demand = hr['PLANT_PDEMAND'][filter, :].resample({'time': '1Y'}).mean()
        som = hr['FUNGI_SOM_TO_PPOOL'][filter, :].resample({'time': '1Y'}).mean()
        for col, add in zip(['hummock', 'hollow'], [0, pft_stride]):
            collect_part2.loc[(col, plot), 'SOMP_TO_SMINP_Spruce'] = (som / smin_demand)[:, 2 + add]
            collect_part2.loc[(col, plot), 'SOMP_TO_SMINP_Tamarack'] = (som / smin_demand)[:, 3 + add]
            collect_part2.loc[(col, plot), 'SOMP_TO_SMINP_Shrub'] = (som / smin_demand)[:, 11 + add]

    if 'FPI' in var_list_extra:
        col_list = ['TOTSOMC', 'FPI', 'FPI_P']
    else:
        col_list = ['TOTSOMC']
    for colvar in col_list: # 'FCH4'
        if growing_season:
            temp = hr2[colvar][filter, :].resample({'time': '1Y'}).mean().values
        else:
            temp = hr2[colvar][:-1, :].resample({'time': '1Y'}).mean().values
        for num, col in enumerate(['hummock', 'hollow']):
            collect_part2.loc[(col, plot), colvar] = temp[:, num]

    # soil N & P in plant accessible layers
    for colvar, thisvar in zip(['SMINN_30','SOLUTIONP_30'], ['SMINN_vr', 'SOLUTIONP_vr']):
        LEVGRND = np.array([0.007100635, 0.027925, 0.06225858, 0.1188651, 0.2121934,
                            0.3660658, 0.6197585, 1.038027, 1.727635, 2.864607, 4.739157,
                            7.829766, 12.92532, 21.32647, 35.17762])
        LEVGRND_I = np.append(np.insert(
            (LEVGRND[1:] + LEVGRND[:-1])*0.5, 0, 0
        ), LEVGRND[-1] + 0.5 * (LEVGRND[-1] - LEVGRND[-2]))
        THICKNESS = np.diff(LEVGRND_I)
        depth = 0.30 # 30 cm
        maxlayer = np.where(LEVGRND_I < depth)[0][-1]

        if growing_season:
            temp = hr2[thisvar][filter, :, :].resample({'time': '1Y'}).mean().values
        else:
            temp = hr2[thisvar][:-1, :, :].resample({'time': '1Y'}).mean().values
        data = np.zeros([temp.shape[0], 2])
        for i in range(maxlayer):
            data = data + temp[:, i, :] * THICKNESS[i]
        last_depth = min(THICKNESS[i], depth - LEVGRND_I[maxlayer])
        data = data + temp[:, maxlayer, :] * last_depth
        data = data / depth
        for num, col in enumerate(['hummock', 'hollow']):
            collect_part2.loc[(col, plot), colvar] = data[:, num]

    # ZWT needs special calculations
    if growing_season or zwt_growing_season:
        temp = hr2['ZWT'][filter, :].resample({'time': '1Y'}).mean()
        h2osfc = hr2['H2OSFC'][filter, :].resample({'time': '1Y'}).mean()
    else:
        temp = hr2['ZWT'][:-1, :].resample({'time': '1Y'}).mean()
        h2osfc = hr2['H2OSFC'][:-1, :].resample({'time': '1Y'}).mean()
    collect_part2.loc[('hummock', plot), 'ZWT'] = 0.3 - temp[:, 0]
    collect_part2.loc[('hollow', plot), 'ZWT'] = h2osfc[:, 1] / 1000 - temp[:, 1]

    hr.close()
    hr2.close()

collect = pd.concat([collect_part1, collect_part2], axis = 1)

temp = (collect.loc['hummock', :] * 0.64 + collect.loc['hollow' , :] * 0.36)
for ind, row in temp.iterrows():
    collect.loc[('average', *ind), :] = row.values

path_out = os.path.join(os.environ['PROJDIR'], 'ELM_Phenology', 'output', 'extract')
if len(extrafix) > 0:
    path_out = os.path.join(path_out, prefix + extrafix)
else:
    path_out = os.path.join(path_out, prefix)
if not os.path.exists(path_out):
    os.mkdir(path_out)

if growing_season:
    collect.to_csv(os.path.join(path_out, 'extract_ts_productivity_gs.csv'))
else:
    collect.to_csv(os.path.join(path_out, 'extract_ts_productivity.csv'))
