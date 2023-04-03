import xarray as xr
from utils.paths import *
from utils.analysis import *
import os
import itertools as it
import pandas as pd
import numpy as np


header = ['Variable', 'Startyear', 'endyear', 'Startday', 'endday', 'averaging period',
          'factor', 'add offset', 'pft', 'obs', 'obs_err', 'treatment']
pft_stride = 17
chamber_list = [6, 19, 20, 11, 13, 4, 8, 16, 17, 10, 7]
chamber_list_names = ['T0.00', 'T0.00CO2', 'T2.25', 'T2.25CO2', 'T4.50', 'T4.50CO2', 'T6.75', 'T6.75CO2', 'T9.00', 'T9.00CO2', 'TAMB']


####################################################################################################################
hr = xr.open_dataset(os.path.join(path_intrim, 'spruce_validation_data.nc'))


# Edited syntax: use comma to separate the factors and pfts
# Simulated = Factor[1] x pft[1] + Factor[2] x pft[2] + ... + add offset


line = ''


# (1) hr['annual_anpp_tree'][7, 21, 6, 19, 20, 11, 13, 4, 8, 16, 17, 10]
for year, chamber in it.product([2016, 2017, 2018, 2019, 2020], chamber_list):
    # pima: 0.36, lala: 0.14
    pairs = (
        (2, 0.36 * 86400 * 365),
        (3, 0.14 * 86400 * 365),
    )
    pft_list = ','.join([str(x[0]) for x in pairs])
    factor_list = ','.join([f'{x[1]:.2f}' for x in pairs])
    obs = float(hr['annual_anpp_tree'].loc[year, chamber])
    treatment = get_treatment_string(chamber)

    line += f'AGNPP\t{year}\t{year}\t1\t365\t365\t{factor_list}\t0\t{pft_list}\t{obs}\t{obs*0.5}\t{treatment}\n'


# (2) hr['annual_anpp_shrub']
for year, chamber in it.product([2016, 2017, 2018, 2019, 2020], chamber_list):
    # hummock: 0.64, hollow: 0.36
    pairs = (
        (11, 0.25 * 86400 * 365),
    )
    pft_list = ','.join([str(x[0]) for x in pairs])
    factor_list = ','.join([f'{x[1]:.2f}' for x in pairs])
    obs = float(hr['annual_anpp_shrub'].loc[year, chamber])
    treatment = get_treatment_string(chamber)

    line += f'AGNPP\t{year}\t{year}\t1\t365\t365\t{factor_list}\t0\t{pft_list}\t{obs}\t{obs*0.5}\t{treatment}\n'


"""
# (3) hr['annual_npp_moss']
for year, chamber in it.product([2016, 2017, 2018, 2019, 2020], chamber_list):
    obs = float(hr['annual_npp_moss'].loc[year, chamber])
    treatment = get_treatment_string(chamber)
    mossfrac = get_mossfrac(year, treatment)

    # hummock: 0.64, hollow: 0.36
    pairs = (
        (12, mossfrac * 86400 * 365),
    )
    pft_list = ','.join([str(x[0]) for x in pairs])
    factor_list = ','.join([f'{x[1]:.2f}' for x in pairs])

    line += f'NPP\t{year}\t{year}\t1\t365\t365\t{factor_list}\t0\t{pft_list}\t{obs}\t{obs*0.5}\t{treatment}\n'
"""

# (4) hr['annual_bnpp']
for year, chamber in it.product([2016, 2017, 2018, 2019, 2020], chamber_list):
    pairs = (
        (2, 0.36 * 86400 * 365),
        (3, 0.14 * 86400 * 365),
        (11, 0.25 * 86400 * 365),
    )
    pft_list = ','.join([str(x[0]) for x in pairs])
    factor_list = ','.join([f'{x[1]:.2f}' for x in pairs])
    obs = float(hr['annual_bnpp'].loc[year, chamber])
    treatment = get_treatment_string(chamber)

    line += f'BGNPP\t{year}\t{year}\t1\t365\t365\t{factor_list}\t0\t{pft_list}\t{obs}\t{obs*0.5}\t{treatment}\n'


"""
# This is irrelevant to the vegetation growth.
# (5) hr['annual_rh']
for year, chamber in it.product([2016, 2017, 2018, 2019, 2020], chamber_list):
    # flip the simulation's sign
    pairs = (
        (0, - 86400 * 365),
    )
    pft_list = ','.join([str(x[0]) for x in pairs])
    factor_list = ','.join([f'{x[1]:.2f}' for x in pairs])
    obs = float(hr['annual_rh'].loc[year, chamber])
    treatment = get_treatment_string(chamber)

    line += f'HR\t{year}\t{year}\t1\t365\t365\t{factor_list}\t0\t{pft_list}\t{obs}\t{obs*0.5}\t{treatment}\n'

# (6) hr['annual_nee']
for year, chamber in it.product([2016, 2017, 2018, 2019, 2020], chamber_list):
    # flip the simulation's sign
    pairs = (
        (0, - 86400 * 365),
    )
    pft_list = ','.join([str(x[0]) for x in pairs])
    factor_list = ','.join([f'{x[1]:.2f}' for x in pairs])
    obs = float(hr['annual_nee'].loc[year, chamber])
    treatment = get_treatment_string(chamber)

    line += f'NEE\t{year}\t{year}\t1\t365\t365\t{factor_list}\t0\t{pft_list}\t{obs}\t{obs*0.5}\t{treatment}\n'

hr.close()
"""


# (7) hr['annual_lai'][7, 21, 6, 19, 20, 11, 13, 4, 8, 16, 17, 10]
for year, chamber in it.product([2016, 2017, 2018, 2019, 2020], chamber_list):
    # pima: 0.36, lala: 0.14, shrub: 0.25, sphagnum: 0.25
    # but do not calibrate the LAI of Sphagnum
    for pft, pft_name, frac in zip([2, 3, 11], ['EN', 'DN', 'SH'], [0.36, 0.14, 0.25, 0.25]):
        factor = frac * 66.4 / 114.8 # make the modeled result compatible with observation

        obs = float(hr['annual_lai'].loc[year, chamber, pft_name])
        
        if np.isnan(obs):
            continue

        treatment = get_treatment_string(chamber)

        if pft == 2 or pft == 12:
            # max LAI is in Oct
            line += f'TLAI\t{year}\t{year}\t274\t304\t31\t{factor}\t0\t{pft}\t{obs}\t{obs*0.25}\t{treatment}\n'
        else:
            # chose Aug
            line += f'TLAI\t{year}\t{year}\t213\t243\t31\t{factor}\t0\t{pft}\t{obs}\t{obs*0.5}\t{treatment}\n'


f = open(os.path.join('./temp/postproc_vars_SPRUCE'), 'w')
f.write('#' + '\t'.join(header) + '\n')
f.write(line)
f.close()