import os
import xarray as xr
import pandas as pd
import numpy as np
from glob import glob
from utils.constants import chamber_list_complete_dict

mossfrac = pd.read_excel(
    "Sphagnum_fraction.xlsx", index_col=0, skiprows=1, engine="openpyxl"
).drop(["plot", "Temp", "CO2"], axis=1)
mossfrac[2015] = mossfrac[2016]

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

#pft_stride = 32
pft_stride = 17
plot_list = ['P04', 'P06', 'P07', 'P08', 'P10', 'P11', 'P13', 'P16', 'P17', 'P19', 'P20']
var_list = ['TOTVEGC_ABG_pima', 'TOTVEGC_ABG_lala', 'TOTVEGC_ABG_shrub', 
            'AGNPP_shrub', 'AGNPP_pima', 'AGNPP_lala', 'AGNPP_tree', 'AGNPP_tree_shrub', 
            'BGNPP_tree_shrub', 'BGNPP_pima', 'BGNPP_lala', 'BGNPP_shrub', 
            'MR_pima', 'MR_lala', 'MR_shrub', 
            'HR', 'NEE', 'NPP_moss', 'TBOT', 'ZWT', 'TOTSOMC', 'SMINN_30', 'SOLUTIONP_30']
var_list_extra = ['FPG_pima', 'FPG_lala', 'FPG_shrub',
                  'FPG_P_pima', 'FPG_P_lala', 'FPG_P_shrub', 
                  'FPI', 'FPI_P']
var_list = var_list + var_list_extra

growing_season = False
zwt_growing_season = True


prefix = "20231113_3"

#prefix  = "UQ_20240315"
extrafix = "" # "_alt_params"

if "UQ" in prefix:
    year_range = range(2015, 2022)
else:
    year_range = range(2015, 2021)

collect = pd.DataFrame(
    np.nan,
    index = pd.MultiIndex.from_product([['hummock', 'hollow', 'average'],
                                        plot_list,
                                        year_range], names = ['column', 'plot', 'year']),
    columns = var_list)

for plot in plot_list:
    if not "UQ" in prefix:
        temp = plot.replace("P", '')    
        rundir = os.path.join(os.environ['PROJDIR'], 'E3SM', 'output', 
                              f'{prefix}_US-SPR_ICB20TRCNPRDCTCBC/spruce_treatments{extrafix}/plot{temp}_US-SPR_ICB20TRCNPRDCTCBC/run')
    else:
        temp = chamber_list_complete_dict[plot]
        rundir = os.path.join(os.environ["PROJDIR"], "E3SM", "output", "UQ", 
                              f"{prefix}_US-SPR_ICB20TRCNPRDCTCBC", "g01000", temp)

    flist = sorted(glob(rundir + "/*.h2.*.nc"))
    if "UQ" in prefix:
        flist = flist[:-1]
    # if len(flist) == 0:
    #    continue
    hr = xr.open_mfdataset(flist)

    flist = sorted(glob(rundir + "/*.h1.*.nc"))
    if "UQ" in prefix:
        flist = flist[:-1]
    hr2 = xr.open_mfdataset(flist)

    if growing_season or zwt_growing_season:
        filter = (hr['time'].to_index().month >= 5) & (hr['time'].to_index().month <= 10)

    # hummock: 0.64, hollow: 0.36
    # pima: 0.36, lala: 0.14
    for var in ['AGNPP', 'MR', 'TOTVEGC_ABG', 'BGNPP']:
        # temporary fix until better values become available
        if growing_season:
            temp = hr[var][filter, :].resample({'time': '1Y'}).mean()
        else:
            temp = hr[var][:-1, :].resample({'time': '1Y'}).mean()
        if not var == 'TOTVEGC_ABG':
            # convert gC/m2/s to gC/m2/year; otherwise gC m-2
            temp = temp.values * 365 * 86400
        else:
            temp = temp.values
        for col, add in zip(['hummock', 'hollow'], [0, pft_stride]):
            collect.loc[(col, plot), f'{var}_pima'] = temp[:, 2 + add] * 0.36
            collect.loc[(col, plot), f'{var}_lala'] = temp[:, 3 + add] * 0.14
            collect.loc[(col, plot), f'{var}_shrub'] = temp[:, 11 + add] * 0.25

    for var in ['AGNPP', 'BGNPP']:
        for col, add in zip(['hummock', 'hollow'], [0, pft_stride]):
            collect.loc[(col, plot), f'{var}_tree'] = \
                collect.loc[(col, plot), f'{var}_pima'].values + \
                collect.loc[(col, plot), f'{var}_lala'].values
            collect.loc[(col, plot), f'{var}_tree_shrub'] = \
                collect.loc[(col, plot), f'{var}_pima'].values + \
                collect.loc[(col, plot), f'{var}_lala'].values + \
                collect.loc[(col, plot), f'{var}_shrub'].values

    if growing_season:
        temp = hr['NPP'][filter, :].resample({'time': '1Y'}).mean()
    else:
        temp = hr['NPP'][:-1, :].resample({'time': '1Y'}).mean()
    # convert gC/m2/s to gC/m2/year
    temp = temp * 365 * 86400
    for col, add in zip(['hummock', 'hollow'], [0, pft_stride]):
        collect.loc[(col, plot), 'NPP_moss'] = temp[:, 12 + add] * \
            mossfrac.loc[plot_to_grid[plot], :].loc[year_range] / 100.

    #if growing_season:
    #    temp = hr['GPP'][filter, :].resample({'time': '1Y'}).mean()
    #else:
    #    temp = hr['GPP'][:-1, :].resample({'time': '1Y'}).mean()
    ## convert gC/m2/s to gC/m2/year
    #temp = temp * 365 * 86400
    #for col, add in zip(['hummock', 'hollow'], [0, pft_stride]):
    #    collect.loc[(col, plot), 'GPP_moss'] = temp[:, 12 + add] * \
    #        mossfrac.loc[plot_to_grid[plot], :].loc[year_range] / 100.

    if 'FPI' in var_list:
        # NP limitation should always focus on growing season
        if 'FPG_PATCH' in hr.variables:
            temp = hr['FPG_PATCH'][filter, :].resample({'time': '1Y'}).mean()
            for col, add in zip(['hummock', 'hollow'], [0, pft_stride]):
                collect.loc[(col, plot), 'FPG_pima'] = temp[:, 2 + add]
                collect.loc[(col, plot), 'FPG_lala'] = temp[:, 3 + add]
                collect.loc[(col, plot), 'FPG_shrub'] = temp[:, 11 + add]
        else:
            temp = hr2['FPG'][filter, :].resample({'time': '1Y'}).mean()
            for num, col in enumerate(['hummock', 'hollow']):
                collect.loc[(col, plot), 'FPG_pima'] = temp[:, num]
                collect.loc[(col, plot), 'FPG_lala'] = temp[:, num]
                collect.loc[(col, plot), 'FPG_shrub'] = temp[:, num]

        if 'FPG_PATCH' in hr.variables:
            temp = hr['FPG_P_PATCH'][filter, :].resample({'time': '1Y'}).mean()
            for col, add in zip(['hummock', 'hollow'], [0, pft_stride]):
                collect.loc[(col, plot), 'FPG_P_pima'] = temp[:, 2 + add]
                collect.loc[(col, plot), 'FPG_P_lala'] = temp[:, 3 + add]
                collect.loc[(col, plot), 'FPG_P_shrub'] = temp[:, 11 + add]
        else:
            temp = hr2['FPG_P'][filter, :].resample({'time': '1Y'}).mean()
            for num, col in enumerate(['hummock', 'hollow']):
                collect.loc[(col, plot), 'FPG_P_pima'] = temp[:, num]
                collect.loc[(col, plot), 'FPG_P_lala'] = temp[:, num]
                collect.loc[(col, plot), 'FPG_P_shrub'] = temp[:, num]

    if 'FPI' in var_list:
        col_list = ['TBOT', 'NEE', 'HR', 'TOTSOMC', 'FPI', 'FPI_P']
    else:
        col_list = ['TBOT', 'NEE', 'HR', 'TOTSOMC']
    for colvar in col_list: # 'FCH4'
        if growing_season:
            temp = hr2[colvar][filter, :].resample({'time': '1Y'}).mean().values
        else:
            temp = hr2[colvar][:-1, :].resample({'time': '1Y'}).mean().values
        if colvar == 'TBOT':
            temp = temp - 273.15
        elif colvar not in ["FPI", "FPI_P", "TOTSOMC"]:
            # convert gC/m2/s to gC/m2/year
            temp = temp * 365 * 86400
        for num, col in enumerate(['hummock', 'hollow']):
            if colvar == 'FCH4':
                collect.loc[(col, plot), 'CH4'] = temp[:, num]
            else:
                collect.loc[(col, plot), colvar] = temp[:, num]
    
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
            temp = hr2[thisvar][filter, :].resample({'time': '1Y'}).mean().values
        else:
            temp = hr2[thisvar][:-1, :].resample({'time': '1Y'}).mean().values
        data = 0.
        for i in range(maxlayer - 1):
            data = temp[:, i, :] * THICKNESS[i]
        last_depth = min(THICKNESS[i], depth - LEVGRND_I[maxlayer])
        data = temp[:, maxlayer, :] * last_depth
        data = data / depth
        for num, col in enumerate(['hummock', 'hollow']):
            collect.loc[(col, plot), colvar] = data[:, num]

    # ZWT needs special calculations
    if growing_season or zwt_growing_season:
        temp = hr2['ZWT'][filter, :].resample({'time': '1Y'}).mean()
        h2osfc = hr2['H2OSFC'][filter, :].resample({'time': '1Y'}).mean()
    else:
        temp = hr2['ZWT'][:-1, :].resample({'time': '1Y'}).mean()
        h2osfc = hr2['H2OSFC'][:-1, :].resample({'time': '1Y'}).mean()
    collect.loc[('hummock', plot), 'ZWT'] = 0.3 - temp[:, 0]
    collect.loc[('hollow', plot), 'ZWT'] = h2osfc[:, 1] / 1000 - temp[:, 1]

    hr.close()
    hr2.close()

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
