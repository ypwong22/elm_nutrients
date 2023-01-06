import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np


gcc = pd.read_csv('austincary_EN_1000_1day.csv', skiprows = 24, index_col = 0,
                   parse_dates = True)['smooth_gcc_90']

lai = pd.read_csv('landsat_lai_austincary.csv', parse_dates = True, index_col = 0)
lai.loc[lai['QA']!=0, 'LAI'] = np.nan

fig, ax = plt.subplots(figsize = (10, 6))
ax.plot(gcc.index, gcc, '-b')
ax2 = ax.twinx()
ax2.plot(lai.index, lai['LAI'], '-r')
ax.set_xlabel('Date')
ax.set_ylabel('GCC')
ax2.set_ylabel('LAI')
fig.savefig('compare_ts.png', dpi = 600, bbox_inches = 'tight')
plt.close(fig)


fig, ax = plt.subplots(figsize = (10, 6))

gcc_cycle = gcc.groupby(gcc.index.dayofyear).mean()
lai_cycle = lai['LAI'].groupby(lai.index.dayofyear).mean()

ax.plot(gcc_cycle.index, gcc_cycle, '-b')
ax2 = ax.twinx()
ax2.plot(lai_cycle.index, lai_cycle, '-r')
ax.set_xlabel('Day of the year')
ax.set_ylabel('GCC')
ax2.set_ylabel('LAI')
fig.savefig('compare_cycle.png', dpi = 600, bbox_inches = 'tight')
plt.close(fig)