""" Pull from the clm_params_g{XXXXX}.nc in each UQ directory, and save back into a mcsamples.txt file
in the exact old OLMT format.
"""

import os
import netCDF4

# ── Paths ──────────────────────────────────────────────────────────────────────
PARM_FILE = os.environ['HOME'] + "/Git/elm_nutrients/calibration_files/parm_file_20231118_compact"
UQ_DIR    = os.environ['E3SM_ROOT'] + "/output/UQ/UQ_20231118_backup/UQ_20231118_US-SPR_ICB1850CNRDCTCBC_ad_spinup"
OUTPUT    = os.environ['HOME'] + "/Git/elm_nutrients/calibration_files/mcsamples_UQ_20231118_4000.txt"

N_SAMPLES = 4000


# ── 1. Read the parameter file ─────────────────────────────────────────────────
params = []  # list of (name, index)
with open(PARM_FILE) as f:
    for line in f:
        parts = line.split()
        if len(parts) >= 2:
            name  = parts[0]
            index = int(parts[1])
            params.append((name, index))

n_params = len(params)
print(f"Read {n_params} parameters from parm_file")


def read_nc_variable(path, varname, index):
    with netCDF4.Dataset(path) as ds:
        var = ds.variables[varname]
        if var.ndim == 0:
            return float(var[()])
        return float(var[index])


# ── 2. Loop through sample directories ────────────────────────────────────────
rows = []
for i in range(1, N_SAMPLES + 1):
    subdir  = f"g{i:05d}"
    nc_file = os.path.join(UQ_DIR, subdir, f"clm_params_{i:05d}.nc")

    if not os.path.isfile(nc_file):
        print(f"WARNING: missing file {nc_file}, filling row with NaN")
        rows.append([float("nan")] * n_params)
        continue

    row = []
    for name, idx in params:
        val = read_nc_variable(nc_file, name, idx)
        row.append(val)

    rows.append(row)

    if i % 500 == 0:
        print(f"  processed {i}/{N_SAMPLES}")

print(f"Finished reading {len(rows)} samples.")


# ── 3. Write mcsamples.txt ─────────────────────────────────────────────────────
# Each column: 24 chars wide, scientific notation, 1 digit left of decimal.
# Columns separated by a single space.
# Format: 20 chars for the float mantissa + 4 chars for e+XX exponent = 24 total.
# Using %.19e gives 1+1+19+1+1+2 = 25 chars for positive values
def fmt_val(v):
    s = f"{v:.18e}"   # e.g. "1.088246959962290e+01" (21 chars for positive)
    return f"{s:>24}"

with open(OUTPUT, "w") as out:
    for row in rows:
        line = " ".join(fmt_val(v) for v in row)
        out.write(line + "\n")

print(f"Wrote {OUTPUT}")
