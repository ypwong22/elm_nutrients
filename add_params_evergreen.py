import xarray as xr
import numpy as np
import os


hr = xr.open_dataset('clm_params.nc_yang_dmr_02242021', decode_times = False)
# hr = xr.open_dataset(os.environ['PROJDIR'] + '/E3SM/output/UQ/UQ_default_US-SPR_ICB20TRCNPRDCTCBC/g00012/clm_params_00012.nc', decode_times = False)

# hr['fcur'][2] = 0.5 // transferring litterfall "storage" means old leaves dying to give to new leaves?
hr['crit_gdd1'][2] = 5.92009
hr['crit_gdd2'][2] = 0.051719
hr['crit_dayl_evergreen'] = xr.DataArray([39454.7], coords = {'allpfts': hr['allpfts']},
                                        dims = ['allpfts'], 
                                        attrs = {'units': 'seconds', 
                                                 'long_name': 'Critical day length for senescence'})
hr['ndays_on_evergreen' ] = xr.DataArray([51.], coords = {'allpfts': hr['allpfts']},
                                        dims = ['allpfts'], 
                                        attrs = {'units': 'days', 
                                                 'long_name': 'Number of days to complete leaf onset'})
hr['ndays_off_evergreen'] = xr.DataArray([14.], coords = {'allpfts': hr['allpfts']},
                                        dims = ['allpfts'], 
                                        attrs = {'units': 'none', 
                                                 'long_name': 'Parameter controlling the number of days of intense litterfall'})
hr['maxday_off_evergreen'] = xr.DataArray([305.], coords = {'allpfts': hr['allpfts']},
                                        dims = ['allpfts'], 
                                        attrs = {'units': 'days', 
                                                 'long_name': 'Day of maximum litterfall'})
encoding = {}
for data_var in hr.data_vars:
    if '_FillValue' in hr[data_var].encoding.keys():
        continue
    else:
        encoding[data_var] = {'_FillValue': None}

hr.to_netcdf('clm_params.nc_yang_dmr_yw_12312022', encoding = encoding, format = 'NETCDF3_CLASSIC')
# hr.to_netcdf('clm_params.nc_yang_dmr_yw_01202023_leaf', encoding = encoding, format = 'NETCDF3_CLASSIC')