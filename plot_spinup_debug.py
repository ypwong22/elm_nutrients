""" Check if the PFT died and the nutrients balance during the spinup
    Requires h0, h1, h2 outputs like spruce treatment runs.
"""
import xarray as xr
import matplotlib.pyplot as plt
from glob import glob
import os


prefix = '20240304_5'
nutrient = 'RD'

# Path to the folder containing the NetCDF files
folder_path = os.path.join(os.environ['PROJDIR'], 'E3SM', 'output',             
                           f'{prefix}_US-SPR_ICB1850CN{nutrient}CTCBC_ad_spinup', 'run')
print(folder_path)

#####################################################################
# PFT growth
#####################################################################
nc_files = sorted(glob(os.path.join(folder_path, "*.h2.*.nc")))[:-1]

ds = xr.open_mfdataset(nc_files, combine='by_coords')

# (1) Patch level nutrient limitation
fpg_monthly = ds['FPG_PATCH'].resample(time='M').mean()
fpg_p_monthly = ds['FPG_P_PATCH'].resample(time='M').mean()

for pft in [2,3,11,12]:
   plt.figure(figsize=(10, 6))
   fpg_monthly[:,pft].plot(label='FPG_PATCH')
   fpg_p_monthly[:,pft].plot(label='FPG_P_PATCH')
   plt.title(pft)
   plt.xlabel('Time')
   plt.ylabel('gC/m2')
   plt.legend()
   plt.show()
   plt.savefig(os.path.join(folder_path, f'fpg_{pft}.png'), dpi = 600.)

# (2) Maintenace respiration
frootc_monthly = ds['FROOT_MR'].resample(time='M').mean()
leafc_monthly = ds['LEAF_MR'].resample(time='M').mean()

for pft in [2,3,11,12]:
   plt.figure(figsize=(10, 6))
   frootc_monthly[:,pft].plot(label='FROOT_MR')
   leafc_monthly[:,pft].plot(label='LEAF_MR')
   plt.title(pft)
   plt.xlabel('Time')
   plt.ylabel('gC/m2')
   plt.legend()
   plt.show()
   plt.savefig(os.path.join(folder_path, f'mr_{pft}.png'), dpi = 600.)

# (3) Vegetation biomass
frootc_monthly = ds['FROOTC'].resample(time='M').mean()
leafc_monthly = ds['LEAFC'].resample(time='M').mean()

for pft in [2,3,11,12]:
   plt.figure(figsize=(10, 6))
   frootc_monthly[:,pft].plot(label='FROOTC')
   leafc_monthly[:,pft].plot(label='LEAFC')
   plt.title(pft)
   plt.xlabel('Time')
   plt.ylabel('gC/m2')
   plt.legend()
   plt.show()
   plt.savefig(os.path.join(folder_path, f'pft_{pft}.png'), dpi = 600.)

ds.close()

#####################################################################
# Microbial nutrients limitation
#####################################################################
nc_files = sorted(glob(os.path.join(folder_path, "*.h1.*.nc")))[:-1]
ds = xr.open_mfdataset(nc_files, combine='by_coords')

# Need to plot this on annual level to understand the constraints
# The monthly variations are too huge
fpi_monthly = ds['FPI'].resample(time='1Y').mean()
fpi_p_monthly = ds['FPI_P'].resample(time='1Y').mean()

plt.figure(figsize=(10, 6))
fpi_monthly[:,0].plot(label='FPI hum')
fpi_p_monthly[:,0].plot(label='FPI_P hum')
fpi_monthly[:,1].plot(label='FPI hol')
fpi_p_monthly[:,1].plot(label='FPI_P hol')
plt.xlabel('Time')
plt.ylabel('gC/m2')
plt.legend()
plt.show()
plt.savefig(os.path.join(folder_path, f'fpi.png'), dpi = 600.)

ds.close()