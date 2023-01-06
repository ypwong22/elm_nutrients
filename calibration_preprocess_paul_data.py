import xarray as xr
from utils.paths import *
from utils.analysis import *
import os
import itertools as it
import pandas as pd


header = ['Variable', 'Startyear', 'endyear', 'Startday', 'endday', 'averaging period',
          'factor', 'add offset', 'pft', 'obs', 'obs_err', 'treatment']
pft_stride = 17


####################################################################################################################
hr = xr.open_dataset(os.path.join(path_intrim, 'spruce_validation_data.nc'))


# Edited syntax: use comma to separate the factors and pfts
# Simulated = Factor[1] x pft[1] + Factor[2] x pft[2] + ... + add offset


line = ''


# (1) hr['annual_anpp_tree']
for year, chamber in it.product([2016, 2017, 2018, 2019, 2020], [7, 21, 6, 19, 20, 11, 13, 4, 8, 16, 17, 10]):
    # pima: 0.36, lala: 0.14
    pairs = (
        (2, 0.36 * 86400 * 365),
        (3, 0.14 * 86400 * 365),
    )
    pft_list = ','.join([str(x[0]) for x in pairs])
    factor_list = ','.join([f'{x[1]:.2f}' for x in pairs])
    obs = float(hr['annual_anpp_tree'].loc[year, chamber])
    treatment = get_treatment_string(chamber)

    line += f'AGNPP\t{year}\t{year}\t1\t365\t365\t{factor_list}\t0\t{pft_list}\t{obs}\t{obs}\t{treatment}\n'


# (2) hr['annual_anpp_shrub']
for year, chamber in it.product([2016, 2017, 2018, 2019, 2020], [7, 21, 6, 19, 20, 11, 13, 4, 8, 16, 17, 10]):
    # hummock: 0.64, hollow: 0.36
    pairs = (
        (11, 0.25 * 86400 * 365),
    )
    pft_list = ','.join([str(x[0]) for x in pairs])
    factor_list = ','.join([f'{x[1]:.2f}' for x in pairs])
    obs = float(hr['annual_anpp_shrub'].loc[year, chamber])
    treatment = get_treatment_string(chamber)

    line += f'AGNPP\t{year}\t{year}\t1\t365\t365\t{factor_list}\t0\t{pft_list}\t{obs}\t{obs}\t{treatment}\n'


# (3) hr['annual_npp_moss']
for year, chamber in it.product([2016, 2017, 2018, 2019, 2020], [7, 21, 6, 19, 20, 11, 13, 4, 8, 16, 17, 10]):
    obs = float(hr['annual_npp_moss'].loc[year, chamber])
    treatment = get_treatment_string(chamber)
    mossfrac = get_mossfrac(year, treatment)

    # hummock: 0.64, hollow: 0.36
    pairs = (
        (12, mossfrac * 86400 * 365),
    )
    pft_list = ','.join([str(x[0]) for x in pairs])
    factor_list = ','.join([f'{x[1]:.2f}' for x in pairs])

    line += f'NPP\t{year}\t{year}\t1\t365\t365\t{factor_list}\t0\t{pft_list}\t{obs}\t{obs}\t{treatment}\n'


# (4) hr['annual_bnpp']
for year, chamber in it.product([2016, 2017, 2018, 2019, 2020], [7, 21, 6, 19, 20, 11, 13, 4, 8, 16, 17, 10]):
    pairs = (
        (2, 0.36 * 86400 * 365),
        (3, 0.14 * 86400 * 365),
        (11, 0.25 * 86400 * 365),
    )
    pft_list = ','.join([str(x[0]) for x in pairs])
    factor_list = ','.join([f'{x[1]:.2f}' for x in pairs])
    obs = float(hr['annual_bnpp'].loc[year, chamber])
    treatment = get_treatment_string(chamber)

    line += f'BGNPP\t{year}\t{year}\t1\t365\t365\t{factor_list}\t0\t{pft_list}\t{obs}\t{obs}\t{treatment}\n'


# This is irrelevant to the vegetation growth.
# (5) hr['annual_rh']
for year, chamber in it.product([2016, 2017, 2018, 2019, 2020], [7, 21, 6, 19, 20, 11, 13, 4, 8, 16, 17, 10]):
    # flip the simulation's sign
    pairs = (
        (0, - 86400 * 365),
    )
    pft_list = ','.join([str(x[0]) for x in pairs])
    factor_list = ','.join([f'{x[1]:.2f}' for x in pairs])
    obs = float(hr['annual_rh'].loc[year, chamber])
    treatment = get_treatment_string(chamber)

    line += f'HR\t{year}\t{year}\t1\t365\t365\t{factor_list}\t0\t{pft_list}\t{obs}\t{obs}\t{treatment}\n'

# (6) hr['annual_nee']
for year, chamber in it.product([2016, 2017, 2018, 2019, 2020], [7, 21, 6, 19, 20, 11, 13, 4, 8, 16, 17, 10]):
    # flip the simulation's sign
    pairs = (
        (0, - 86400 * 365),
    )
    pft_list = ','.join([str(x[0]) for x in pairs])
    factor_list = ','.join([f'{x[1]:.2f}' for x in pairs])
    obs = float(hr['annual_nee'].loc[year, chamber])
    treatment = get_treatment_string(chamber)

    line += f'NEE\t{year}\t{year}\t1\t365\t365\t{factor_list}\t0\t{pft_list}\t{obs}\t{obs}\t{treatment}\n'

hr.close()


f = open(os.path.join(path_out, 'postproc_vars_SPRUCE'), 'w')
f.write('#' + '\t'.join(header) + '\n')
f.write(line)
f.close()