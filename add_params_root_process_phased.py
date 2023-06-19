import xarray as xr
import numpy as np
import os
from utils.paths import *

path_parameter = os.path.join(os.environ['SCRATCHDIR'], 'E3SM', 'inputdata', 'atm', 'datm7', 'CLM1PT_data', 'SPRUCE_data')

orgfile = 'clm_params.nc_yang_dmr_20230518_UQ_default2_optimized_XGBClassifier_MLPRegressor'
newfile = 'clm_params.nc_yang_dmr_20230518_UQ_default2_optimized_XGBClassifier_MLPRegressor_root_phased'

hr = xr.open_dataset(os.path.join(path_parameter, orgfile), decode_times = False)

hr['fcur'][2] = 0.8
hr['fcur'][3] = 0.
hr['fcur'][11] = 0.


hr['gdd_tbase'] = xr.DataArray([np.nan, np.nan, 278.15, 279.5, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 279.05, 
                                 np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                                coords = {'pft': hr['pft']},
                                dims = ['pft'],
                                attrs = {'units': 'K',
                                        'long_name': 'Base temperature for GDD accumulation'})
hr['crit_chil1'] = xr.DataArray([np.nan, np.nan, 199.30927021132482, 9., np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 33, 
                                 np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                                coords = {'pft': hr['pft']},
                                dims = ['pft'], 
                                attrs = {'units': '',
                                        'long_name': 'Parameter in alternating model'})
hr['crit_chil2'] = xr.DataArray([np.nan, np.nan, 3751.8473351447456, 2112, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 1388, 
                                 np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                                coords = {'pft': hr['pft']},
                                dims = ['pft'], 
                                attrs = {'units': '',
                                        'long_name': 'Parameter in alternating model'})
hr['crit_chil3'] = xr.DataArray([np.nan, np.nan, -0.03728095348592797, -0.04, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, -0.02, 
                                 np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                                coords = {'pft': hr['pft']},
                                dims = ['pft'], 
                                attrs = {'units': '',
                                        'long_name': 'Parameter in alternating model'})
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
hr['ndays_on_root'] = xr.DataArray([np.nan, np.nan, 90, 90, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 90, np.nan, np.nan,
                                    np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                                   coords = {'pft': hr['pft']},
                                   dims = ['pft'],
                                   attrs = {'units': 'days',
                                            'long_name': 'Number of days to deplete root storage by onset growth'})
hr['ndays_off'] = xr.DataArray([np.nan, np.nan, 14, 15, np.nan, np.nan, 15, 15, 15, np.nan, 15, 15, np.nan, 15, 15, 15, 15, 
                                np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                                coords = {'pft': hr['pft']},
                                dims = ['pft'], 
                                attrs = {'units': 'days', 
                                         'long_name': 'Number of days to complete leaf litterfall (for evergreen this is only the relevant parameter not exact)'})
# might need an ndays_off_root here
hr['maxday_off'] = xr.DataArray([320.], coords = {'allpfts': hr['allpfts']},
                                dims = ['allpfts'], 
                                attrs = {'units': 'days', 
                                         'long_name': 'Day of maximum litterfall'})

hr['crit_h2osoi'] = xr.DataArray([0.6], dims = ['allpft'], attrs = {'units': '', 'long_name': 'Parameter in soil moisture control on root growth and mortality'})
hr['crit_psi'] = xr.DataArray([35], dims = ['allpft'], attrs = {'units': '', 'long_name': 'Parameter in soil moisture control on root growth and mortality'})
hr['crit_tsoi'] = xr.DataArray([10], dims = ['allpft'], attrs = {'units': 'degC', 'long_name': 'Parameter in soil moisture control on root growth and mortality'})
hr['nmin_beta'] = xr.DataArray([1], dims = ['allpft'], attrs = {'units': '', 'long_name': 'Parameter in soil moisture control on root growth and mortality'})
hr['sw_scale'] = xr.DataArray([100], dims = ['allpft'], attrs = {'units': 'mm', 'long_name': 'Parameter in soil moisture control on root growth and mortality'})


# additional decidous phenology parameters
hr['off_pstart'] = xr.DataArray([np.nan, np.nan, np.nan, 46800., np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 54600., 
                                 np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                                coords = {'pft': hr['pft']},
                                dims = ['pft'], 
                                attrs = {'units': '',
                                        'long_name': 'Paramter for deciduous offset'})
hr['off_pend'] = xr.DataArray([np.nan, np.nan, np.nan, 1750., np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 1600., 
                                 np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                                coords = {'pft': hr['pft']},
                                dims = ['pft'], 
                                attrs = {'units': '',
                                        'long_name': 'Paramter for deciduous offset'})
hr['off_tbase'] = xr.DataArray([np.nan, np.nan, np.nan, 294.5, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 290.15, 
                                np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                                coords = {'pft': hr['pft']},
                                dims = ['pft'], 
                                attrs = {'units': '',
                                         'long_name': 'Paramter for deciduous offset'})


encoding = {}
for data_var in hr.data_vars:
    if '_FillValue' in hr[data_var].encoding.keys():
        continue
    elif np.any(np.isnan(hr[data_var].values)):
        encoding[data_var] = {'_FillValue': -1e20}
    else:
        encoding[data_var] = {'_FillValue': None}
hr.to_netcdf(os.path.join(path_parameter, newfile), encoding = encoding, format = 'NETCDF3_CLASSIC')