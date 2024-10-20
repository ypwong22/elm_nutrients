"""
20241014

We want to understand how to improve the response of 
 - AGNPP's mean and slope
 - HR's mean and slope
 to warming (CO2 elevation is an add-on), under reasonable
 levels of N saturation in the plants (FPG & FPG_P) and 
 pore water N. 

Knowledge from UQ_20240311_1 & UQ_20240312_testXXX
1. Higher shrub AGNPP translates to higher HR, when you just vary
   shrub vmax_froot_n/vmax_fungi_son + km_froot_n + fungi_cost_n. 
2. vmax_fungi_son seems to do nothing. This might be because fungi
   uptake fraction is too low in UQ_20240311_1. But UQ_20240312_testXXX
   bore this out. Very strange! We may need to vary vmax_fungi_son's 
   ratio to vmax_froot_n by at least one order of magnitude to see effect.

Knowledge from UQ_20240311_1
1. Shrub ANPP's sensitivity to temperature increases with shrub's
   vmax_froot_n, but then levels after vmax_froot_n hits 1.2e-11. 
   Shrub ANPP mean continus to increase a bit, with the observation being
   ~80 (ACO2)/~110 (ECO2) gC m-2, with the sweet spot more like 2e-11
   (depending on how high you set your vmax_froot_n and fungi_cost_n,
    but you probably want those to be high for HR reasons; this means 
    2e-11 to 3e-11). 

   Using tuned parameters from the default model, we are unlikely to need
   as high shrub ANPP. 

   Therefore, we probably only need to test up to 1e-11, with high
   km_froot_n (>10) & fungi_cost_n (>=100) values. 

   If we want to use a sample generator based on ratio of km_froot_n/vmax_froot_n
   and fungi_cost_n/vmax_froot_n, log-ratio works better than just ratio.
   
   log(km_froot_n / vmax_froot_n) between 25-30 seems okay.
   log(fungi_cost_n / vmax_froot_n) we don't seem to want those to drop below 26. 
   26-32 seems very reasonable. 

Knowledge from UQ_20240311_1 & UQ_20240312_test20241013
1. Higher km_froot_n decreases shrub AGNPP (bad), but does not affect shrub 
   AGNPP's sensitivity to temperature (good). It decreases HR and its 
   sensitivity to warming (good). Over 4-15 range, the effects are linear. 
2. Higher fungi_cost_n decreases shrub AGNPP (bad), increases shrub ANPP's
   sensitivity to temperature (good). It dcreases HR and it's sensitivity
   to temperature (good). Over 10-100 range, the effects are linear. 

Knowledge from UQ_20240311_2
"""
import numpy as np
import os
import xarray as xr
import itertools as it


nsamples = 3125
n_parameters = 1 + 1 + 3 + 3 + 3

# km_froot_n: 4 - 12
# fungi_cost_n: 80 - 160
# vmax_froot_n, spruce: 5e-12 to 2.5e-10
# vmax_froot_n, tamarack: 4e-11 to 2e-9
# vmax_froot_n, shrub: 8e-12 to 4e-10
# ratio of vmax_fungi_son to vmax_froot_n: keep constant at current file's levels

ratio_fungi_din = [224, 521, 557]
ratio_fungi_son = [703, 105, 3224228]

samples = np.zeros((nsamples, n_parameters), dtype = float)

count = 0
for km_froot_n, fungi_cost_n, vfn_spr, vfn_tama, vfn_shrub in \
   it.product(np.linspace(4, 12, 5), np.linspace(80, 160, 5), 
              np.linspace(5e-12, 2.5e-10, 5), np.linspace(4e-11, 2e-9, 5), 
              np.linspace(8e-12, 4e-10, 5)):

   samples[count, 0] = km_froot_n
   samples[count, 1] = fungi_cost_n
   samples[count, 2] = vfn_spr
   samples[count, 3] = vfn_tama
   samples[count, 4] = vfn_shrub
   samples[count, 5] = vfn_spr * ratio_fungi_din[0]
   samples[count, 6] = vfn_tama * ratio_fungi_din[1]
   samples[count, 7] = vfn_shrub * ratio_fungi_din[2]
   samples[count, 8] = vfn_spr * ratio_fungi_son[0]
   samples[count, 9] = vfn_tama * ratio_fungi_son[1]
   samples[count, 10] = vfn_shrub * ratio_fungi_son[2]

   count = count + 1

np.savetxt('./calibration_files/mcsamples_20240312_test20241017.txt', samples)
