import xarray as xr
import numpy as np

hr = xr.open_dataset('clm_params.nc_yang_dmr_02242021', decode_times = False)
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

for i, stopdate in enumerate([365, 171]):
    hr['stop_date'] = xr.DataArray([stopdate],
        coords = {'allpfts': hr['allpfts']}, dims = ['allpfts'], 
        attrs =  {'units': 'none', 
                  'long_name': 'Day of the year when storage-driven root growth stops'})

    hr['proot_onset_tbase'] = xr.DataArray([ \
        np.nan, np.nan, 0., 0., np.nan, 
        np.nan, np.nan, np.nan, np.nan, np.nan,
        np.nan, 0., np.nan, np.nan, np.nan, 
        np.nan, np.nan, np.nan, np.nan, np.nan, 
        np.nan, np.nan, np.nan, np.nan, np.nan],
        coords = {'pft': hr['pft']}, dims = ['pft'],
        attrs = {'units': 'none', 'long_name': 'Base temperature for calculating the soil temperature controlled onset of growth'})
    hr['proot_onset_tbase'].encoding['_FillValue'] = -999.99

    hr['proot_onset_intercept'] = xr.DataArray([ \
        np.nan, np.nan, 25.348249556614952, 25.348249556614952, np.nan, 
        np.nan, np.nan, np.nan, np.nan, np.nan,
        np.nan, 21.854824068918475, np.nan, np.nan, np.nan, 
        np.nan, np.nan, np.nan, np.nan, np.nan, 
        np.nan, np.nan, np.nan, np.nan, np.nan],
        coords = {'pft': hr['pft']}, dims = ['pft'],
        attrs = {'units': 'none', 'long_name': 'Critical growing degree days for temperature-controlled root growth onset'})
    hr['proot_onset_intercept'].encoding['_FillValue'] = -999.99

    hr['proot_onset_slope'] = xr.DataArray([ \
        np.nan, np.nan, -0.08860982483861916, -0.08860982483861916, np.nan, 
        np.nan, np.nan, np.nan, np.nan, np.nan,
        np.nan, -0.12558838304033987, np.nan, np.nan, np.nan, 
        np.nan, np.nan, np.nan, np.nan, np.nan, 
        np.nan, np.nan, np.nan, np.nan, np.nan],
        coords = {'pft': hr['pft']}, dims = ['pft'],
        attrs = {'units': 'none', 'long_name': 'Critical growing degree days for temperature-controlled root growth onset'})
    hr['proot_onset_slope'].encoding['_FillValue'] = -999.99

    hr['proot_npp_frac'] = xr.DataArray([ \
        np.nan, np.nan, 0.99, 0.99, np.nan, 
        np.nan, np.nan, np.nan, np.nan, np.nan,
        np.nan, 0.99, np.nan, np.nan, np.nan, 
        np.nan, np.nan, np.nan, np.nan, np.nan, 
        np.nan, np.nan, np.nan, np.nan, np.nan],
        coords = {'pft': hr['pft']}, dims = ['pft'],
        attrs = {'units': 'none', 'long_name': 'Fraction of displayed root carbon out of all root carbon allocated from net primary productivity'})
    hr['proot_npp_frac'].encoding['_FillValue'] = -999.99

    if stopdate == 365:
        for switch in range(1,5):
            hr['switch']    = xr.DataArray([switch], 
                coords = {'allpfts': hr['allpfts']}, dims = ['allpfts'], 
                attrs  = {'units': 'none', 
                          'long_name': 'Control which root phenology scheme to use'})

            if switch == 1:
                hr['proot_a'] = xr.DataArray([ \
                    np.nan, np.nan, 0.01724373420691372, 0.01724373420691372, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan,
                    np.nan, -0.013919512505618499, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan],
                    coords = {'pft': hr['pft']}, dims = ['pft'],
                    attrs = {'units': 'none', 'long_name': 'first parameter in the root growth phenology'})

                hr['proot_b'] = xr.DataArray([ \
                    np.nan, np.nan, 1.7757576078230444, 1.7757576078230444, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan,
                    np.nan, 0.7989308814801414, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan],
                    coords = {'pft': hr['pft']}, dims = ['pft'],
                    attrs = {'units': 'none', 'long_name': 'second parameter in the root growth phenology'})

            elif switch == 2:
                hr['proot_a'] = xr.DataArray([ \
                    np.nan, np.nan, 0.015217772115630676, 0.015217772115630676, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan,
                    np.nan, -0.012735805399980716, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan],
                    coords = {'pft': hr['pft']}, dims = ['pft'],
                    attrs = {'units': 'none', 'long_name': 'first parameter in the root growth phenology'})
                hr['proot_b'] = xr.DataArray([ \
                    np.nan, np.nan, 0.7971664109443781, 0.7971664109443781, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan,
                    np.nan, 1.3319086767688797, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan],
                    coords = {'pft': hr['pft']}, dims = ['pft'],
                    attrs = {'units': 'none', 'long_name': 'second parameter in the root growth phenology'})

            elif switch == 3:
                hr['proot_a'] = xr.DataArray([ \
                    np.nan, np.nan, 0.04216671197151266, 0.04216671197151266, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan,
                    np.nan, 0.038769065865924386, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan],
                    coords = {'pft': hr['pft']}, dims = ['pft'],
                    attrs = {'units': 'none', 'long_name': 'first parameter in the root growth phenology'})
                hr['proot_b'] = xr.DataArray([ \
                    np.nan, np.nan, 1.2193508483196647, 1.2193508483196647, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan,
                    np.nan, 1.233444342100033, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan],
                    coords = {'pft': hr['pft']}, dims = ['pft'],
                    attrs = {'units': 'none', 'long_name': 'second parameter in the root growth phenology'})

            elif switch == 4:
                hr['proot_a'] = xr.DataArray([ \
                    np.nan, np.nan, 0.3488034767961882, 0.3488034767961882, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan,
                    np.nan, 0.41669626687671196, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan],
                    coords = {'pft': hr['pft']}, dims = ['pft'],
                    attrs = {'units': 'none', 'long_name': 'first parameter in the root growth phenology'})
                hr['proot_b'] = xr.DataArray([ \
                    np.nan, np.nan, 1.4142135623730951, 1.4142135623730951, np.nan,
                    np.nan, np.nan, np.nan, np.nan, np.nan,
                    np.nan, 1.4142135623730951, np.nan, np.nan, np.nan,
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
            hr.to_netcdf(f'clm_params.nc_yang_dmr_yw_04052022_stopdate{stopdate}_switch{switch}', encoding = encoding, format = 'NETCDF3_CLASSIC')

    else:
        for switch in range(1,5):
            hr['switch']     = xr.DataArray([switch], 
                coords = {'allpfts': hr['allpfts']}, dims = ['allpfts'], 
                attrs = {'units': 'none', 
                         'long_name': 'Control which root phenology scheme to use'})

            if switch == 1:
                hr['proot_a'] = xr.DataArray([ \
                    np.nan, np.nan, 0.027797606206855563, 0.027797606206855563, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan,
                    np.nan, -1.6435917288672965e-05, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan],
                    coords = {'pft': hr['pft']}, dims = ['pft'],
                    attrs = {'units': 'none', 'long_name': 'first parameter in the root growth phenology'})
                hr['proot_b'] = xr.DataArray([ \
                    np.nan, np.nan, 0.22074491604111696, 0.22074491604111696, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan,
                    np.nan, -1.0050318083152572, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan],
                    coords = {'pft': hr['pft']}, dims = ['pft'],
                    attrs = {'units': 'none', 'long_name': 'second parameter in the root growth phenology'})

            elif switch == 2:
                hr['proot_a'] = xr.DataArray([ \
                    np.nan, np.nan, -0.022382370109498437, -0.022382370109498437, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan,
                    np.nan, 218.35011544555346, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan],
                    coords = {'pft': hr['pft']}, dims = ['pft'],
                    attrs = {'units': 'none', 'long_name': 'first parameter in the root growth phenology'})
                hr['proot_b'] = xr.DataArray([ \
                    np.nan, np.nan, 0.5782882711084616, 0.5782882711084616, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan,
                    np.nan, np.nan, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan],
                    coords = {'pft': hr['pft']}, dims = ['pft'],
                    attrs = {'units': 'none', 'long_name': 'second parameter in the root growth phenology'})

            elif switch == 3:
                hr['proot_a'] = xr.DataArray([ \
                    np.nan, np.nan, 0.03931414498134363, 0.03931414498134363, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan,
                    np.nan, 0.012099852706207014, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan],
                    coords = {'pft': hr['pft']}, dims = ['pft'],
                    attrs = {'units': 'none', 'long_name': 'first parameter in the root growth phenology'})
                hr['proot_b'] = xr.DataArray([ \
                    np.nan, np.nan, 1.5290218443090728, 1.5290218443090728, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan,
                    np.nan, np.nan, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan],
                    coords = {'pft': hr['pft']}, dims = ['pft'],
                    attrs = {'units': 'none', 'long_name': 'second parameter in the root growth phenology'})

            elif switch == 4:
                hr['proot_a'] = xr.DataArray([ \
                    np.nan, np.nan, 0.31170893555680557, 0.31170893555680557, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan,
                    np.nan, 0.34624944092803683, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan],
                    coords = {'pft': hr['pft']}, dims = ['pft'],
                    attrs = {'units': 'none', 'long_name': 'first parameter in the root growth phenology'})
                hr['proot_b'] = xr.DataArray([ \
                    np.nan, np.nan, 0.8944271909999159, 0.8944271909999159, np.nan, 
                    np.nan, np.nan, np.nan, np.nan, np.nan,
                    np.nan, 0.8944271909999159, np.nan, np.nan, np.nan, 
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
            hr.to_netcdf(f'clm_params.nc_yang_dmr_yw_04052022_stopdate{stopdate}_switch{switch}', encoding = encoding, format = 'NETCDF3_CLASSIC')
