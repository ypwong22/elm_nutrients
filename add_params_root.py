import xarray as xr
import numpy as np
import os


hr = xr.open_dataset('clm_params.nc_yang_dmr_02242021', decode_times = False)
# hr = xr.open_dataset(os.environ['PROJDIR'] + '/E3SM/output/UQ/UQ_default_US-SPR_ICB20TRCNPRDCTCBC/g00012/clm_params_00012.nc', decode_times = False)

hr['crit_gdd1'][2] = 5.92009 # parameter 1 in the critical gdd function
hr['crit_gdd2'][2] = 0.051719 # parameter 2 in the critical gdd function
hr['crit_dayl_evergreen' ] = xr.DataArray([39454.7], coords = {'allpfts': hr['allpfts']},
                                          dims = ['allpfts'], 
                                          attrs = {'units': 'seconds',
                                                   'long_name': 'Critical day length for senescence'})
hr['ndays_on_evergreen'  ] = xr.DataArray([51.], coords = {'allpfts': hr['allpfts']},
                                          dims = ['allpfts'], 
                                          attrs = {'units': 'days', 
                                                   'long_name': 'Number of days to complete leaf onset'})
hr['ndays_off_evergreen' ] = xr.DataArray([14.], coords = {'allpfts': hr['allpfts']},
                                          dims = ['allpfts'], 
                                          attrs = {'units': 'none', 
                                                   'long_name': 'Parameter controlling the number of days of intense litterfall'})
hr['maxday_off_evergreen'] = xr.DataArray([305.], coords = {'allpfts': hr['allpfts']},
                                          dims = ['allpfts'], 
                                          attrs = {'units': 'days', 
                                                   'long_name': 'Day of maximum litterfall'})
hr['ndays_on_root'       ] = xr.DataArray([131.], coords = {'allpfts': hr['allpfts']},
                                          dims = ['allpfts'], 
                                          attrs = {'units': 'days', 
                                                   'long_name': 'Number of days to complete root onset'})


hr['proot_a'] = xr.DataArray([ \
    np.nan, np.nan, 0.027987847634371753, 0.027987847634371753, np.nan, 
    np.nan, np.nan, np.nan, np.nan, np.nan,
    np.nan, -0.013590746753969844, np.nan, np.nan, np.nan, 
    np.nan, np.nan, np.nan, np.nan, np.nan, 
    np.nan, np.nan, np.nan, np.nan, np.nan],
    coords = {'pft': hr['pft']}, dims = ['pft'],
    attrs = {'units': 'none', 'long_name': 'first parameter in the root growth phenology'})

hr['proot_b'] = xr.DataArray([ \
    np.nan, np.nan, 7.582252148060458, 7.582252148060458, np.nan, 
    np.nan, np.nan, np.nan, np.nan, np.nan,
    np.nan, 0.8284980816958448, np.nan, np.nan, np.nan, 
    np.nan, np.nan, np.nan, np.nan, np.nan, 
    np.nan, np.nan, np.nan, np.nan, np.nan],
    coords = {'pft': hr['pft']}, dims = ['pft'],
    attrs = {'units': 'none', 'long_name': 'second parameter in the root growth phenology'})

hr['proot_a'].encoding['_FillValue'] = -999.99
hr['proot_b'].encoding['_FillValue'] = -999.99

encoding = {}
for data_var in hr.data_vars:
    if '_FillValue' in hr[data_var].encoding.keys():
        continue
    else:
        encoding[data_var] = {'_FillValue': None}

hr.to_netcdf(f'clm_params.nc_yang_dmr_yw_20230101', encoding = encoding, format = 'NETCDF3_CLASSIC')
# hr.to_netcdf('clm_params.nc_yang_dmr_yw_20230120_root', encoding = encoding, format = 'NETCDF3_CLASSIC')