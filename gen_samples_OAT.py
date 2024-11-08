"""Generate one at a time perturbation in all of the parameters
+/-50%, 10 samples per parameter
"""
import os
import xarray as xr
import numpy as np

prefix = '20240407'
basic_file = os.path.join(
    os.environ['E3SM_ROOT'], 'output', f'{prefix}_US-SPR_ICB1850CNRDCTCBC_ad_spinup', 'run', 
    'clm_params.nc'
)
hr = xr.open_dataset(basic_file)


f = open(os.path.join('calibration_files', f'parm_file_{prefix}_OAT'), 'w')

########################################################################
# PFT-specific parameters
########################################################################
for parname in ['vmax_froot_n', 'vmax_froot_p', 'vmax_fungi_din', 'vmax_fungi_dip', 
                'vmax_fungi_son', 'vmax_fungi_sop']:
    for pft in [2,3,11]:
        f.write(parname + ' ' + str(pft) + ' ' + str(float(hr[parname][pft]*0.5)) + ' ' \
                + str(float(hr[parname][pft]*1.5)) + '\n')
f.write(f'zwt_fungi_a 2 0.3 0.9\n')
f.write(f'zwt_fungi_a 3 0.3 0.9\n')
f.write(f'zwt_fungi_a 11 0.4 1\n')
f.write(f'zwt_fungi_b 2 -0.1 0.1\n')
f.write(f'zwt_fungi_b 3 -0.1 0.1\n')
f.write(f'zwt_fungi_b 11 0.4 1\n')
f.write(f'zwt_froot_a 11 0 2\n')

########################################################################
# Non-PFT parameters
########################################################################
for parname in ['km_froot_n', 'km_froot_p', 'alpha_fpg', 'km_nsc', 'fungi_cost_n', 'fungi_cost_p',
                'q10_upt', 'swc_opt']:
    f.write(parname + ' ' + str(pft) + ' ' + str(float(hr[parname]*0.5)) + ' ' \
            + str(float(hr[parname]*1.5)) + '\n')

f.close()

########################################################################
# Put everything into the same mc samples file
########################################################################
n_samples = 15

# default values
samples_default = np.full(33, np.nan)
count = 0
for parname in ['vmax_froot_n', 'vmax_froot_p', 'vmax_fungi_din', 'vmax_fungi_dip', 
                'vmax_fungi_son', 'vmax_fungi_sop', 'zwt_fungi_a', 'zwt_fungi_b']:
    for pft in [2,3,11]:
        samples_default[count] = hr[parname][pft]
        count = count + 1
samples_default[count] = hr['zwt_froot_a'][pft]
count = count + 1
for parname in ['km_froot_n', 'km_froot_p', 'alpha_fpg', 'km_nsc', 'fungi_cost_n', 'fungi_cost_p',
                'q10_upt', 'swc_opt']:
    samples_default[count] = hr[parname]
    count = count + 1

# samples
samples = (np.broadcast_to(samples_default.reshape(1, -1), (33*n_samples, 33))).copy()
count = 0
for parname in ['vmax_froot_n', 'vmax_froot_p', 'vmax_fungi_din', 'vmax_fungi_dip', 
                'vmax_fungi_son', 'vmax_fungi_sop']:
    for pft in [2,3,11]:
        samples[(count*n_samples):((count+1)*n_samples), count] = np.linspace(hr[parname][pft]*0.5, 
                                                            hr[parname][pft]*1.5, n_samples)
        count = count + 1
samples[(count*n_samples):((count+1)*n_samples), count] = np.linspace(0.3, 0.9, n_samples)
count = count + 1
samples[(count*n_samples):((count+1)*n_samples), count] = np.linspace(0.3, 0.9, n_samples)
count = count + 1
samples[(count*n_samples):((count+1)*n_samples), count] = np.linspace(0.4, 1, n_samples)
count = count + 1
samples[(count*n_samples):((count+1)*n_samples), count] = np.linspace(-0.1, 0.1, n_samples)
count = count + 1
samples[(count*n_samples):((count+1)*n_samples), count] = np.linspace(-0.1, 0.1, n_samples)
count = count + 1
samples[(count*n_samples):((count+1)*n_samples), count] = np.linspace(0.4, 1, n_samples)
count = count + 1
samples[(count*n_samples):((count+1)*n_samples), count] = np.linspace(0, 2, n_samples)
count = count + 1
for parname in ['km_froot_n', 'km_froot_p', 'alpha_fpg', 'km_nsc', 'fungi_cost_n', 'fungi_cost_p',
                'q10_upt', 'swc_opt']:
    samples[(count*n_samples):((count+1)*n_samples), count] = \
        np.linspace(float(hr[parname]*0.5), float(hr[parname]*1.5), n_samples)
    count = count + 1

np.savetxt(os.path.join('calibration_files', f'mcsamples_{prefix}_OAT.txt'), samples)


f.close()


hr.close()