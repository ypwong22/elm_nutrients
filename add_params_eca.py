import xarray as xr
import os

inputfile = os.path.join(os.environ['E3SM_ROOT'], 'inputdata', 'lnd', 'clm2', 'paramdata',
                         'clm_params.eca.c190307.nc')
h_in = xr.open_dataset(inputfile)

hr = xr.open_dataset(os.path.join(os.environ['E3SM_ROOT'], 'inputdata',
                                  'atm', 'datm7', 'CLM1PT_data', 'SPRUCE_data',
                                  'clm_params_SPRUCE_20231120_spruceroot.nc'
                                  ), decode_times = False) # 'clm_params.nc_yang_dmr_02242021'
hr['alpha_nfix'] = h_in['alpha_nfix']
hr['alpha_ptase'] = h_in['alpha_ptase']
hr['VMAX_PTASE'] = h_in['VMAX_PTASE']

# Need to use SPRUCE-specific chemistry 
hr['leafcn_obs'] = h_in['leafcn'] # hr['leafcn']
hr['leafcn_obs_flex'] = h_in['leafcn_obs_flex'] # - h_in['leafcn_obs'] + hr['leafcn']
hr['frootcn_obs'] = h_in['frootcn_obs'] # hr['frootcn']
hr['frootcn_obs_flex'] = h_in['frootcn_obs_flex'] # / h_in['frootcn_obs'] * hr['frootcn']
hr['livewdcn_obs'] = h_in['livewdcn_obs'] # hr['livewdcn']
hr['livewdcn_obs_flex'] = h_in['livewdcn_obs_flex'] # / h_in['livewdcn_obs'] * hr['livewdcn']
hr['deadwdcn_obs'] = h_in['deadwdcn_obs'] # hr['deadwdcn']
hr['deadwdcn_obs_flex'] = h_in['deadwdcn_obs_flex'] # / h_in['deadwdcn_obs'] * hr['deadwdcn']

hr['leafcp_obs'] = hr['leafcp']
hr['leafcp_obs_flex'] = h_in['leafcp_obs_flex'] - h_in['leafcp_obs'] + hr['leafcp']
hr['frootcp_obs'] = h_in['frootcp_obs'] # hr['frootcp']
hr['frootcp_obs_flex'] = h_in['frootcp_obs_flex'] # / h_in['frootcp_obs'] * hr['frootcp']
hr['livewdcp_obs'] = h_in['livewdcp_obs'] # hr['livewdcp']
hr['livewdcp_obs_flex'] = h_in['livewdcp_obs_flex'] # / h_in['livewdcp_obs'] * hr['livewdcp']
hr['deadwdcp_obs'] = h_in['deadwdcp_obs'] # hr['deadwdcp']
hr['deadwdcp_obs_flex'] = h_in['deadwdcp_obs_flex'] # / h_in['deadwdcp_obs'] * hr['deadwdcp']

hr['vcmax_np1'] = h_in['vcmax_np1']
hr['vcmax_np2'] = h_in['vcmax_np2']
hr['vcmax_np3'] = h_in['vcmax_np3']
hr['vcmax_np4'] = h_in['vcmax_np4']

hr['jmax_np1'] = h_in['jmax_np1']
hr['jmax_np2'] = h_in['jmax_np2']
hr['jmax_np3'] = h_in['jmax_np3']

hr['laimax'] = h_in['laimax']
hr['laimx'] = h_in['laimx'] # red flag: those are NaN, perhaps not good?

encoding = {}
for data_var in hr.data_vars:
    if '_FillValue' in hr[data_var].encoding.keys():
        continue
    else:
        encoding[data_var] = {'_FillValue': None}

outputfile = os.path.join(os.environ['E3SM_ROOT'], 'inputdata',
                          'atm', 'datm7', 'CLM1PT_data', 'SPRUCE_data',
                          'clm_params_SPRUCE_20231120_spruceroot_eca.nc') 
                          # 'clm_params.nc_eca_yang_dmr_02242021')
hr.to_netcdf(outputfile, encoding = encoding, format = 'NETCDF3_CLASSIC')

hr.close()
h_in.close()