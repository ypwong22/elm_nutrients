from tokenize import PlainToken
import xarray as xr
import os
import pandas as pd
import numpy as np
from glob import glob
from datetime import datetime, timedelta
from utils.analysis import *
import matplotlib.pyplot as plt


####################################################################
# The vegetation patch
####################################################################
p = 4 # [3, 20, 4, 12, 21, 29] [evergreen grid 1, evergreen grid 2]
      #                        [deciduous tree, deciduous shrub][grid 1]
      #                        [deciduous tree, deciduous shrub][grid 2]
if p in [3,4,12]:
    g = 0
else:
    g = 1

####################################################################
# Constants
####################################################################
secspday  = 86400
secspstep = 3600
stop_date = 365
fracday   = 1/24

####################################################################
# Read the input dataset - soil temperature time series; use another simulations, 
# since the switch2 simulation codes are problematic
####################################################################
path_root = '/lustre/haven/user/ywang254/E3SM/output/20220407_stopdate365_switch1_plot17_US-SPR_ICB20TRCNPRDCTCBC/'
hr = xr.open_mfdataset(glob(os.path.join(path_root, 'run', '*.h0.*.nc')))
tsoi = hr['TSOI'][:, 3, g].values
tvec = pd.DatetimeIndex([datetime(i.year, i.month, i.day) for i in hr['time'].to_index()])
hr.close()

# speed things up by calculating these beforehand
# non-simple daylength calculation is too slow
day_of_year = pd.DatetimeIndex([datetime(1995, i,j) for i,j in zip(tvec.month, tvec.day)]).dayofyear
dayl        = np.array([daylength_simple(i, 47.563) for i in day_of_year]) * 3600
prev_dayl   = np.array(
    [daylength_simple(365, 47.563)] * 24 + 
    [daylength_simple(i, 47.563) for i in day_of_year[:-24]]) * 3600

# Offset days; use evergreen/evergreen +/- 15 days for simplicity
offset_flag = np.zeros(len(tsoi))
offset_counter = np.zeros(len(tsoi))
temp = 0
for i in range(len(tvec)):
    if p in [3,20]:
        crit_dayl = 39454.7
        ndays_off_evergreen = 30
    else:
        crit_dayl = 655
        ndays_off_evergreen = 30

    if (dayl[i] <= crit_dayl) & (prev_dayl[i] > crit_dayl):
        offset_flag[i] = 1
        offset_counter[i] = ndays_off_evergreen * secspday
        temp = offset_counter[i]

    if temp > 0:
        offset_flag[i] = 1
        offset_counter[i] = temp - secspstep
        temp = temp - secspstep
    else:
        offset_flag[i] = 0

# check
fig, axes = plt.subplots(1, 2, figsize = (8,6))
ax = axes.flat[0]
ax.plot(tvec, offset_flag, '-k')
ax.set_title('offset_flag')
ax = axes.flat[1]
ax.plot(tvec, offset_counter, '-k')
ax.set_title('offset_counter')
fig.savefig(os.path.join(f'/lustre/haven/proj/UTK0134/Phenology_ELM/check_leaf_offset_{p}.png'),
            dpi = 600.)
plt.close(fig)

####################################################################
# Read the parameter dataset
####################################################################
hr = xr.open_dataset(f'/lustre/haven/user/ywang254/E3SM/inputdata/atm/datm7/CLM1PT_data/SPRUCE_data/clm_params.nc_yang_dmr_yw_04052022_stopdate{stop_date}_switch2')
ndays_off_evergreen = hr['ndays_off_evergreen'].values[0]
proot_onset_tbase = hr['proot_onset_tbase'].values[p]
proot_onset_intercept = hr['proot_onset_intercept'].values[p]
proot_onset_slope = hr['proot_onset_slope'].values[p]
proot_a = hr['proot_a'].values[p]
proot_b = hr['proot_b'].values[p]
hr.close()

scaling_factor = 1 / \
    (np.exp(proot_a * stop_date) / (proot_b + np.exp(proot_a * stop_date)) - \
     1 / (proot_b + 1))
print(scaling_factor)

####################################################################
# Pre-allocate the variables following Fortran format
####################################################################
proot_gdd1flag_save        = np.full(len(tsoi), np.nan)
proot_gdd1_save            = np.full(len(tsoi), np.nan)
proot_onset_offset_save    = np.full(len(tsoi), np.nan)
proot_onset_flag_save      = np.full(len(tsoi), np.nan)
proot_onset_counter_save   = np.full(len(tsoi), np.nan)
proot_onset_counter2_save  = np.full(len(tsoi), np.nan)
xfer_percentage_save       = np.full(len(tsoi), np.nan)
proot_sumfrac_save         = np.full(len(tsoi), np.nan)
prev_ws_flag_save          = np.full(len(tsoi), np.nan)


proot_gdd1flag             = np.nan
proot_gdd1                 = np.nan
proot_onset_offset         = np.nan
proot_onset_flag           = np.nan
proot_onset_counter        = np.nan
proot_onset_counter2       = np.nan
xfer_percentage            = np.nan
proot_sumfrac              = np.nan
prev_ws_flag               = np.nan

####################################################################
# The algorithm
####################################################################
for i in range(len(tvec)):
    # This equation is wrong for evergreen!!! ndays_off_evergreen does not indicate
    # the total length of litterfall
    if (offset_flag[i] == 1) & \
       (offset_counter[i] == ndays_off_evergreen*secspday):
        proot_gdd1flag = 1
        proot_gdd1 = 0.

    if (dayl[i] >= prev_dayl[i]):
        ws_flag = 1
    else:
        ws_flag = 0

    if (ws_flag == 0) & (proot_gdd1flag == 1):
        proot_gdd1 = proot_gdd1 + max(tsoi[i] - proot_onset_tbase - 273.15, 0) * fracday

    if (prev_ws_flag == 0) & (ws_flag == 1) & (proot_gdd1flag == 1):
        proot_onset_offset = int(proot_onset_intercept + proot_onset_slope * proot_gdd1)
        proot_onset_flag   = 1 # indicate that the onset is pending switch

        # if the result is earlier than winter solstice, set to winter solstice
        if (day_of_year[i] > (proot_onset_offset + 365)):
            proot_onset_offset = day_of_year[i] - 365
        proot_gdd1 = 0
        proot_gdd1flag = 0

    if (proot_onset_flag == 1) and \
       (((day_of_year[i] - 365) == proot_onset_offset) or \
        (day_of_year[i] == proot_onset_offset)):
        if (proot_onset_counter > 0):
            proot_onset_counter2 = proot_onset_counter
        proot_onset_counter = stop_date * secspday
        proot_sumfrac = 0
        proot_onset_flag = 0 # ensure the switch is called only once

    if (proot_onset_counter > 0):
        if (proot_onset_counter2 > 0):
            xfer_percentage = proot_a * proot_b * \
                np.exp((stop_date * secspday - proot_onset_counter2)/secspday*proot_a) / \
                np.power(np.exp((stop_date * secspday - proot_onset_counter2)/secspday*proot_a) + \
                         proot_b, 2) * scaling_factor * fracday + \
                                 proot_a * proot_b * \
                np.exp((stop_date * secspday - proot_onset_counter)/secspday*proot_a) / \
                np.power(np.exp((stop_date * secspday - proot_onset_counter)/secspday*proot_a) + \
                         proot_b, 2) * scaling_factor * fracday
            proot_onset_counter2 = proot_onset_counter2 - secspstep
        else:
            xfer_percentage = proot_a * proot_b * \
                np.exp((stop_date * secspday - proot_onset_counter)/secspday*proot_a) / \
                np.power(np.exp((stop_date * secspday - proot_onset_counter)/secspday*proot_a) + \
                         proot_b, 2) * scaling_factor * fracday

        proot_onset_counter = proot_onset_counter - secspstep
    else:
        xfer_percentage = 0.

    proot_sumfrac = proot_sumfrac + xfer_percentage
    prev_ws_flag = ws_flag

    proot_gdd1flag_save[i]        = proot_gdd1flag
    proot_gdd1_save[i]            = proot_gdd1
    proot_onset_flag_save[i]      = proot_onset_flag
    proot_onset_offset_save[i]    = proot_onset_offset
    proot_onset_counter_save[i]   = proot_onset_counter
    proot_onset_counter2_save[i]  = proot_onset_counter2
    xfer_percentage_save[i]       = xfer_percentage
    proot_sumfrac_save[i]         = proot_sumfrac
    prev_ws_flag_save[i]          = prev_ws_flag


####################################################################
# Display the outcome
####################################################################
varnames = ['proot_gdd1flag', 'proot_gdd1', 'proot_onset_offset', 
            'proot_onset_counter', 'proot_onset_counter2', 'xfer_percentage', 
            'proot_sumfrac', 'prev_ws_flag']
fig, axes = plt.subplots(3, 3, figsize = (12, 12))
for i, var in enumerate(
    [proot_gdd1flag_save, proot_gdd1_save, proot_onset_offset_save, 
     proot_onset_counter_save, proot_onset_counter2_save, xfer_percentage_save, 
     proot_sumfrac_save, prev_ws_flag_save]):

    ax = axes.flat[i]
    ax.plot(tvec, var, '-')
    ax.set_title(varnames[i])

    if varnames[i] == 'proot_sumfrac':
        ax.set_ylim([0,1])
fig.savefig(os.path.join(f'/lustre/haven/proj/UTK0134/Phenology_ELM/check_result_{p}.png'), dpi = 600.)
plt.close(fig)