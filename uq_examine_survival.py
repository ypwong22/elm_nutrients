"""Check the PFT survival situation compared to parameter values during spin-up"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

PARM_SAMPLES = os.environ['HOME'] + '/Git/elm_nutrients/calibration_files/mcsamples_20240101_OAT.txt'
UQ_DIR       = os.environ['E3SM_ROOT'] + '/output/UQ/20260224_US-SPR_ICB1850CNMYCICTCBC_ad_spinup_MYCI_OAT'

# 1. Load parameter samples
# Each column = one parameter, each row = one sample
params = np.loadtxt(PARM_SAMPLES)
if params.ndim == 1:
    params = params[:, np.newaxis]

n_samples = params.shape[0]
print(f"Loaded {n_samples} samples, {params.shape[1]} parameters")

# 2. For each row, check subfolder g{row:05d} and evaluate TLAI survival criterion.
#    Pairs to sum: (2,34), (3,35), (11,43), (12,44) — indices into TLAI's last dimension.
#    If the minimum of all pair-sums is non-zero → 1 (survived), else → 0 (died).
#    Missing subfolder or <2 matching files → NaN.

SURVIVAL_PAIRS = [(2, 34), (3, 35), (11, 43), (12, 44)]

survival = np.full(n_samples, np.nan)

for row in range(n_samples):
    subfolder = os.path.join(UQ_DIR, f'g{row:05d}')
    if not os.path.isdir(subfolder):
        continue  # leave as NaN

    files = sorted(glob.glob(os.path.join(subfolder, '*.h1.*-00000.nc')))
    if len(files) < 2:
        continue  # skip, leave as NaN

    target_file = files[-2]  # second-to-last file

    with nc.Dataset(target_file) as ds:
        tlai = ds.variables['TLAI'][:]  # shape: (time, n) or (n,)

    # Fill masked values with NaN
    if np.ma.is_masked(tlai):
        tlai = tlai.filled(np.nan)

    # Ensure 2D: (time, n)
    if tlai.ndim == 1:
        tlai = tlai[np.newaxis, :]

    # For each pair, compute element-wise sum across the last axis and find the minimum over time
    pair_min_sums = [
        (tlai[:, a]*0.64 + tlai[:, b]*0.36).min()
        for a, b in SURVIVAL_PAIRS
    ]

    survival[row] = 1 if min(pair_min_sums) >1e-3 else 0

n_survived = int(np.nansum(survival == 1))
n_died     = int(np.nansum(survival == 0))
n_nan      = int(np.sum(np.isnan(survival)))
print(f"Survival — survived: {n_survived}, died: {n_died}, NaN: {n_nan}")

# 3. Visualize: boxplot of each parameter's values split by outcome (0 vs 1)
n_params = params.shape[1]
nrows = int(np.sqrt(n_params))
ncols = int(np.ceil(n_params/nrows))

fig, axes = plt.subplots(nrows, ncols, figsize=(3*ncols, 3*nrows), squeeze=False)
axes = axes.flat

for i, ax in enumerate(axes):
    if i >= n_params:
        continue

    plot_data   = []
    plot_labels = []

    data_died = params[survival == 0, i]
    if len(data_died) > 0:
        plot_data.append(data_died)
        plot_labels.append(f'Died\n(n={len(data_died)})')

    data_survived = params[survival == 1, i]
    if len(data_survived) > 0:
        plot_data.append(data_survived)
        plot_labels.append(f'Survived\n(n={len(data_survived)})')

    ax.boxplot(plot_data, labels=plot_labels)
    ax.set_title(f'Parameter {i + 1}')
    ax.set_ylabel('Parameter value')

plt.suptitle('Parameter distributions by PFT survival outcome')
plt.tight_layout()
plt.savefig('survival_boxplot.png', dpi=150, bbox_inches='tight')
plt.show()
print("Saved survival_boxplot.png")
