import xarray as xr
import numpy as np
import os


hr = xr.open_dataset('clm_params.nc_yang_dmr_02242021', decode_times = False)

hr['fcur'][2] = 0.5
hr['fcur'][3] = 0.5
hr['fcur'][11] = 0.5

hr['crit_gdd1'][2] = 5.92009 # parameter 1 in the critical gdd function
hr['crit_gdd2'][2] = 0.051719 # parameter 2 in the critical gdd function
hr['crit_dayl'] = xr.DataArray([np.nan, np.nan, 39454.7, 39300, np.nan, np.nan, 36000, 39300, 39300, np.nan, 36000, 39300, 
                                 np.nan, 36000, 36000, 36000, 36000, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                                coords = {'pft': hr['pft']},
                                dims = ['pft'], 
                                attrs = {'units': 'seconds',
                                        'long_name': 'Critical day length for senescence'})
hr['ndays_on'] = xr.DataArray([np.nan, np.nan, 51, 30, np.nan, np.nan, 30, 30, 30, np.nan, 30, 30, np.nan, 30, 30, 30, 30, 
                               np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                               coords = {'pft': hr['pft']},
                               dims = ['pft'],
                               attrs = {'units': 'days',
                                        'long_name': 'Number of days to complete leaf onset'})
hr['ndays_off'] = xr.DataArray([np.nan, np.nan, 14, 15, np.nan, np.nan, 15, 15, 15, np.nan, 15, 15, np.nan, 15, 15, 15, 15, 
                                np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                                coords = {'pft': hr['pft']},
                                dims = ['pft'], 
                                attrs = {'units': 'days', 
                                         'long_name': 'Number of days to complete leaf litterfall (for evergreen this is only the relevant parameter not exact)'})
hr['maxday_off'] = xr.DataArray([320.], coords = {'allpfts': hr['allpfts']},
                                dims = ['allpfts'], 
                                attrs = {'units': 'days', 
                                         'long_name': 'Day of maximum litterfall'})
hr['rstor2tran'] = xr.DataArray([180.], coords = {'allpfts': hr['allpfts']},
                                dims = ['allpfts'], 
                                attrs = {'units': 'days', 
                                         'long_name': 'Number of days it takes to deplete the storage pool.'})

encoding = {}
for data_var in hr.data_vars:
    if '_FillValue' in hr[data_var].encoding.keys():
        continue
    else:
        encoding[data_var] = {'_FillValue': None}
encoding['crit_dayl'] = {'_FillValue': -1e20}
encoding['ndays_on']  = {'_FillValue': -1e20}
encoding['ndays_off'] = {'_FillValue': -1e20}

hr.to_netcdf(f'clm_params.nc_yang_dmr_yw_20230212', encoding = encoding, format = 'NETCDF3_CLASSIC')