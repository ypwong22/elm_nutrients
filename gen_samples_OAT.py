"""Generate one at a time perturbation in all of the parameters
+/-50%, 10 samples per parameter
"""
import os
import xarray as xr
import numpy as np

prefix = '20260224' # '20240106'
basic_file = os.path.join(
    os.environ['E3SM_ROOT'], 'inputdata', 'atm', 'datm7', 'CLM1PT_data', 'SPRUCE_data',
    'clm_params_SPRUCE_20231120_spruceroot.nc_npcompet' # _slatop'
)
hr = xr.open_dataset(basic_file)


f = open(os.path.join('calibration_files', f'parm_file_{prefix}_OAT'), 'w')

########################################################################
# PFT-specific parameters
########################################################################
for parname in ['vmax_froot_n', 'vmax_froot_p', 'vmax_fungi_din', 'vmax_fungi_dip', 
                'vmax_fungi_son', 'vmax_fungi_sop']:
    for pft in [2,3,11]:
        f.write(parname + ' ' + str(pft) + ' ' + str(float(hr[parname][pft]*0.1)) + ' ' \
                + str(float(hr[parname][pft]*10)) + '\n')
for parname in ['km_froot_n', 'km_froot_p']:
    for pft in [2,3,11]:
        f.write(parname + ' ' + str(pft) + ' ' + str(float(hr[parname][pft]*0.5)) + ' ' \
                + str(float(hr[parname][pft]*2)) + '\n')
f.write(f'inh_fungi_a 2 0.2 0.8\n')
f.write(f'inh_fungi_a 3 0.2 0.8\n')
f.write(f'inh_fungi_a 11 0.4 1\n')
f.write(f'inh_fungi_b 2 -0.1 0.1\n')
f.write(f'inh_fungi_b 3 -0.1 0.1\n')
f.write(f'inh_fungi_b 11 -0.1 0\n')
f.write(f'zwt_froot_a 11 0 2\n')
f.write(f'fungi_cost_n 0 0 20\n')
f.write(f'fungi_cost_p 0 0 200\n')
f.write(f'q10_upt 0 1 3\n')
f.write(f'swc_opt 0 0.3 0.9\n')
f.write(f'alpha_fpg 0 0.75 2.25\n')
f.write(f'km_nsc 0 1 3\n')

f.close()

########################################################################
# Put everything into the same mc samples file
# Need to make sure the order is correct!!!
########################################################################
n_samples = 50
npars     = 37

# default values
samples_default = np.full(npars, np.nan)
count = 0
for parname in ['vmax_froot_n', 'vmax_froot_p', 'vmax_fungi_din', 'vmax_fungi_dip', 
                'vmax_fungi_son', 'vmax_fungi_sop', 'km_froot_n', 'km_froot_p',
                'inh_fungi_a', 'inh_fungi_b']:
    for pft in [2,3,11]:
        samples_default[count] = hr[parname][pft]
        count = count + 1
samples_default[count] = hr['zwt_froot_a'][pft]
count = count + 1
for parname in ['fungi_cost_n', 'fungi_cost_p', 'q10_upt', 'swc_opt', 'alpha_fpg', 'km_nsc']:
    samples_default[count] = hr[parname]
    count = count + 1

# samples
samples = (np.broadcast_to(samples_default.reshape(1, -1), (npars*n_samples, npars))).copy()
count = 0
for parname in ['vmax_froot_n', 'vmax_froot_p', 'vmax_fungi_din', 'vmax_fungi_dip', 
                'vmax_fungi_son', 'vmax_fungi_sop']:
    for pft in [2,3,11]:
        samples[(count*n_samples):((count+1)*n_samples), count] = \
            np.linspace(hr[parname][pft]*0.1, hr[parname][pft]*10, n_samples)
        count = count + 1
for parname in ['km_froot_n', 'km_froot_p']:
    for pft in [2,3,11]:
        samples[(count*n_samples):((count+1)*n_samples), count] = \
            np.linspace(hr[parname][pft]*0.5, hr[parname][pft]*2, n_samples)
        count = count + 1
samples[(count*n_samples):((count+1)*n_samples), count] = np.linspace(0.2, 0.8, n_samples)
count = count + 1
samples[(count*n_samples):((count+1)*n_samples), count] = np.linspace(0.2, 0.8, n_samples)
count = count + 1
samples[(count*n_samples):((count+1)*n_samples), count] = np.linspace(0.4, 1, n_samples)
count = count + 1
samples[(count*n_samples):((count+1)*n_samples), count] = np.linspace(-0.1, 0.1, n_samples)
count = count + 1
samples[(count*n_samples):((count+1)*n_samples), count] = np.linspace(-0.1, 0.1, n_samples)
count = count + 1
samples[(count*n_samples):((count+1)*n_samples), count] = np.linspace(-0.1, 0, n_samples)
count = count + 1
samples[(count*n_samples):((count+1)*n_samples), count] = np.linspace(0, 2, n_samples)
count = count + 1
samples[(count*n_samples):((count+1)*n_samples), count] = np.linspace(10, 100, n_samples)
count = count + 1
samples[(count*n_samples):((count+1)*n_samples), count] = np.linspace(100, 1000, n_samples)
count = count + 1
samples[(count*n_samples):((count+1)*n_samples), count] = np.linspace(1, 3, n_samples)
count = count + 1
samples[(count*n_samples):((count+1)*n_samples), count] = np.linspace(0.3, 0.9, n_samples)
count = count + 1
samples[(count*n_samples):((count+1)*n_samples), count] = np.linspace(0.75, 2.25, n_samples)
count = count + 1
samples[(count*n_samples):((count+1)*n_samples), count] = np.linspace(1, 3, n_samples)
count = count + 1

np.savetxt(os.path.join('calibration_files', f'mcsamples_{prefix}_OAT.txt'), samples)


f.close()


hr.close()