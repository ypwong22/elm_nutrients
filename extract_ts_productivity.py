import os
import xarray as xr
import pandas as pd
import numpy as np
from glob import glob


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

growing_season = False
zwt_growing_season = True

#prefix = "20231112"
prefix = "20240227"
collect = pd.DataFrame(
    np.nan,
    index = pd.MultiIndex.from_product([['hummock', 'hollow', 'average'],
                                        plot_list,
                                        range(2015, 2021)], names = ['column', 'plot', 'year']),
    columns = ['AGNPP_shrub', 'AGNPP_tree', 'BGNPP_tree_shrub', 'GPP_moss', 
                'HR', 'NEE', 'NPP_moss', 'TBOT', 'ZWT']
)

for plot in plot_list:
    temp = plot.replace("P", '')
    rundir = os.path.join(os.environ['PROJDIR'], 'E3SM', 'output', 
                          f'{prefix}_US-SPR_ICB20TRCNPRDCTCBC/spruce_treatments/plot{temp}_US-SPR_ICB20TRCNPRDCTCBC/run')

    flist = sorted(glob(rundir + "/*.h2.*.nc"))
    # if len(flist) == 0:
    #    continue

    hr = xr.open_mfdataset(flist)

    if growing_season or zwt_growing_season:
        filter = (hr['time'].to_index().month >= 5) & (hr['time'].to_index().month <= 10)

    # hummock: 0.64, hollow: 0.36
    # pima: 0.36, lala: 0.14
    if growing_season:
        temp = hr['AGNPP'][filter, :].resample({'time': '1Y'}).mean()
    else:
        temp = hr['AGNPP'][:-1, :].resample({'time': '1Y'}).mean()
    # convert gC/m2/s to gC/m2/year
    temp = temp * 365 * 86400
    for col, add in zip(['hummock', 'hollow'], [0, pft_stride]):
        collect.loc[(col, plot), 'AGNPP_tree'] = \
            temp[:, 2 + add] * 0.36 + temp[:, 3 + add] * 0.14
        collect.loc[(col, plot), 'AGNPP_shrub'] = temp[:, 11 + add] * 0.25

    if growing_season:
        temp = hr["BGNPP"][filter, :].resample({'time': '1Y'}).mean()
    else:
        temp = hr["BGNPP"][:-1, :].resample({'time': '1Y'}).mean() # hr['NPP'] - hr['AGNPP']
    # convert gC/m2/s to gC/m2/year
    temp = temp * 365 * 86400
    for col, add in zip(['hummock', 'hollow'], [0, pft_stride]):
        collect.loc[(col, plot), 'BGNPP_tree_shrub'] = \
            temp[:, 2 + add] * 0.36 + temp[:, 3 + add] * 0.14 + temp[:, 11 + add] * 0.25

    if growing_season:
        temp = hr['NPP'][filter, :].resample({'time': '1Y'}).mean()
    else:
        temp = hr['NPP'][:-1, :].resample({'time': '1Y'}).mean()
    # convert gC/m2/s to gC/m2/year
    temp = temp * 365 * 86400
    for col, add in zip(['hummock', 'hollow'], [0, pft_stride]):
        collect.loc[(col, plot), 'NPP_moss'] = temp[:, 12 + add] * \
            mossfrac.loc[plot_to_grid[plot], :].loc[range(2015,2021)] / 100.

    if growing_season:
        temp = hr['GPP'][filter, :].resample({'time': '1Y'}).mean()
    else:
        temp = hr['GPP'][:-1, :].resample({'time': '1Y'}).mean()
    # convert gC/m2/s to gC/m2/year
    temp = temp * 365 * 86400
    for col, add in zip(['hummock', 'hollow'], [0, pft_stride]):
        collect.loc[(col, plot), 'GPP_moss'] = temp[:, 12 + add] * \
            mossfrac.loc[plot_to_grid[plot], :].loc[range(2015,2021)] / 100.

    hr.close()

    flist = sorted(glob(rundir + "/*.h1.*.nc"))
    hr = xr.open_mfdataset(flist)

    if growing_season or zwt_growing_season:
        filter2 = (hr['time'].to_index().month >= 5) & (hr['time'].to_index().month <= 10)

    for colvar in ['TBOT', 'NEE', 'HR']: # 'FCH4'
        if growing_season:
            temp = hr[colvar][filter2, :].resample({'time': '1Y'}).mean()
        else:
            temp = hr[colvar][:-1, :].resample({'time': '1Y'}).mean()
        if colvar == 'TBOT':
            temp = temp - 273.15
        else:
            # convert gC/m2/s to gC/m2/year
            temp = temp * 365 * 86400
        for num, col in enumerate(['hummock', 'hollow']):
            if colvar == 'FCH4':
                collect.loc[(col, plot), 'CH4'] = temp[:, num]
            else:
                collect.loc[(col, plot), colvar] = temp[:, num]

    if growing_season or zwt_growing_season:
        temp = hr['ZWT'][filter2, :].resample({'time': '1Y'}).mean()
        h2osfc = hr['H2OSFC'][filter2, :].resample({'time': '1Y'}).mean()
    else:
        temp = hr['ZWT'][:-1, :].resample({'time': '1Y'}).mean()
        h2osfc = hr['H2OSFC'][:-1, :].resample({'time': '1Y'}).mean()
    collect.loc[('hummock', plot), 'ZWT'] = 0.3 - temp[:, 0]
    collect.loc[('hollow', plot), 'ZWT'] = h2osfc[:, 1] / 1000 - temp[:, 1]

    hr.close()

temp = (collect.loc['hummock', :] * 0.64 + collect.loc['hollow' , :] * 0.36)
for ind, row in temp.iterrows():
    collect.loc[('average', *ind), :] = row.values

if growing_season:
    collect.to_csv(os.path.join(os.environ['PROJDIR'], 'ELM_Phenology', 'output', 'extract', prefix, 'extract_ts_productivity_gs.csv'))
else:
    collect.to_csv(os.path.join(os.environ['PROJDIR'], 'ELM_Phenology', 'output', 'extract', prefix, 'extract_ts_productivity.csv'))
