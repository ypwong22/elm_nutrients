import os
from datetime import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


casename = "test_US-SPR_ICB1850CNRDCTCBC_ad_spinup"
lnd_log = "lnd.log.4060371.230619-114723"
first_year = 28
last_year = 35

filename = os.path.join(
    os.environ["PROJDIR"], "E3SM", "output", casename, "run", lnd_log
)


colerr = pd.DataFrame(
    np.nan,
    columns=["hollow", "hummock"],
    index=pd.date_range(f"20{first_year}-01-01", f"20{last_year}-04-18"),
)


keeptrack = 0
f = open(filename, encoding="latin-1")
while True:
    line = f.readline()
    if len(line) == 0:
        break

    if line.startswith(" Beginning timestep   :"):
        lsplit = int(line.split()[-1].split("-")[0])
        if lsplit < first_year:
            continue
        if lsplit > last_year:
            break

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

        count = 0
        while count < 2:
            found = False
            while not found:
                line = f.readline()
                if line.startswith(" column cbalance error "):
                    lerr = float(line.split()[-2])

                    if count == 0:
                        patch = "hummock"
                    else:
                        patch = "hollow"

                    colerr.loc[date, patch] = lerr
                    found = True
            count = count + 1

        keeptrack = keeptrack + 1

f.close()


#
fig, ax = plt.subplots()
ax.plot(colerr.index, colerr.iloc[:, 0], "-o")
ax.set_ylabel("column cbalance error")
ax.set_xlabel("mm-dd hh")
fig.savefig("test4.png")
