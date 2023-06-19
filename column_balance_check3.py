import os
from datetime import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


casename = "test_US-SPR_ICB1850CNRDCTCBC_ad_spinup"
lnd_log = "lnd.log.4056315.230616-192217"
first_year = 28

filename = os.path.join(
    os.environ["PROJDIR"], "E3SM", "output", casename, "run", lnd_log
)

delta_state = pd.DataFrame(
    np.nan,
    columns=pd.MultiIndex.from_product(
        [
            ["hollow", "hummock"],
            [2, 3, 11, 12],
            ["totpftc", "xsmrpool", "ctrunc", "dispvegc", "storvegc"],
        ]
    ),
    index=pd.date_range(f"20{first_year}-01-01", "2029-05-20"),
)
flux = pd.DataFrame(
    np.nan,
    columns=pd.MultiIndex.from_product(
        [
            ["hollow", "hummock"],
            [2, 3, 11, 12],
            [
                "gpp",
                "ar",
                "litfall",
                "mr",
                "gr",
                "xr",
                "cpool_gr",
                "transfer_gr",
                "stor_gr",
            ],
        ]
    ),
    index=pd.date_range(f"20{first_year}-01-01", "2029-05-20"),
)


keeptrack = 0
f = open(filename, encoding="latin-1")
while True:
    line = f.readline()
    if len(line) == 0:
        break

    if line.startswith(" Beginning timestep   :"):
        lsplit = line.split()[-1].split("-")
        if int(lsplit[0]) < first_year:
            continue

        lsplit = line.split()[-1]

        if keeptrack % (24 * 100) == 0:
            print(lsplit)

        lsplit = datetime.strptime(lsplit, "%Y-%m-%d_%H:%M:%S")
        date = datetime(
            lsplit.year + 2000,
            lsplit.month,
            lsplit.day,
            lsplit.hour,
            lsplit.minute,
            lsplit.second,
        )

        line = f.readline()
        start = {}
        while line.startswith(" Start State"):
            lsplit = line.split()
            start[lsplit[2]] = np.array(lsplit[3:]).astype(float)
            line = f.readline()

        for count in range(6):
            line = f.readline()

        end = {}
        while line.startswith(" End State"):
            lsplit = line.split()
            end[lsplit[2]] = np.array(lsplit[3:]).astype(float)
            line = f.readline()

        for k in start.keys():
            diff = end[k] - start[k]
            pft = np.mod(int(k), 17) - 1
            if int(k) > 17:
                patch = "hollow"
            else:
                patch = "hummock"
            delta_state.loc[date, (patch, pft)] = diff

        while line.startswith(" Flux"):
            lsplit = line.split()
            pft = np.mod(int(lsplit[1]), 17) - 1
            if int(lsplit[1]) > 17:
                patch = "hollow"
            else:
                patch = "hummock"
            flux.loc[date, (patch, pft)] = np.array(lsplit[2:]).astype(float)
            line = f.readline()

        keeptrack = keeptrack + 1

f.close()


#
fig, axes = plt.subplots(4, 2, figsize=(15, 15))
for i, pft in enumerate([2, 3, 11, 12]):
    for j, patch in enumerate(["hummock", "hollow"]):
        ax = axes[i, j]
        ax.plot(
            flux.index,
            delta_state[(patch, pft, "totpftc")]
            - (
                flux[(patch, pft, "gpp")]
                - flux[(patch, pft, "ar")]
                - flux[(patch, pft, "litfall")]
            )
            * 3600,
            "-",
        )
        ax.set_title("end C - start C - flux C")
fig.savefig("test.png", dpi=600.0)
plt.close(fig)


fig, axes = plt.subplots(4, 2, figsize=(15, 15))
for i, pft in enumerate([2, 3, 11, 12]):
    for j, patch in enumerate(["hummock", "hollow"]):
        ax = axes[i, j]
        ax.plot(
            flux.index,
            flux[(patch, pft, "transfer_gr")],
            "-",
        )
        ax.set_title("end C - start C - flux C")
fig.savefig("test2.png", dpi=600.0)
plt.close(fig)
