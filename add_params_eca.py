import xarray as xr
import os


inputfile = os.path.join(os.environ['PROJDIR'], 'E3SM', 'inputdata', 'lnd', 'clm2', 'paramdata', 
                         'clm_params.eca.c190307.nc')
h_in = xr.open_dataset(inputfile)

hr = xr.open_dataset(os.path.join(os.environ['PROJDIR'], 'E3SM', 'inputdata',
                                  'atm', 'datm7', 'CLM1PT_data', 'SPRUCE_data',
                                  'clm_params.nc_yang_dmr_02242021'), decode_times = False)
hr['alpha_nfix'] = h_in['alpha_nfix']
hr['alpha_ptase'] = h_in['alpha_ptase']
hr['VMAX_PTASE'] = h_in['VMAX_PTASE']

hr['leafcn_ob'] = h_in['leafcn_obs']
hr['leafcn_obs_flex'] = h_in['leafcn_obs_flex']
hr['frootcn_ob'] = h_in['frootcn_obs']
hr['frootcn_obs_flex'] = h_in['frootcn_obs_flex']
hr['livewdcn_obs'] = h_in['livewdcn_obs']
hr['livewdcn_obs_flex'] = h_in['livewdcn_obs_flex']
hr['deadwdcn_obs'] = h_in['deadwdcn_obs']
hr['deadwdcn_obs_flex'] = h_in['deadwdcn_obs_flex']

hr['leafcp_obs'] = h_in['leafcp_obs']
hr['leafcp_obs_flex'] = h_in['leafcp_obs_flex']
hr['frootcp_obs'] = h_in['frootcp_obs']
hr['frootcp_obs_flex'] = h_in['frootcp_obs_flex']
hr['livewdcp_obs'] = h_in['livewdcp_obs']
hr['livewdcp_obs_flex'] = h_in['livewdcp_obs_flex']
hr['deadwdcp_obs'] = h_in['deadwdcp_obs']
hr['deadwdcp_obs_flex'] = h_in['deadwdcp_obs_flex']

hr['vcmax_np1'] = h_in['vcmax_np1']
hr['vcmax_np2'] = h_in['vcmax_np2']
hr['vcmax_np3'] = h_in['vcmax_np3']
hr['vcmax_np4'] = h_in['vcmax_np4']

hr['jmax_np1'] = h_in['jmax_np1']
hr['jmax_np2'] = h_in['jmax_np2']
hr['jmax_np3'] = h_in['jmax_np3']

hr['laimax'] = h_in['laimax']
hr['laimx'] = h_in['laimx']

encoding = {}
for data_var in hr.data_vars:
    if '_FillValue' in hr[data_var].encoding.keys():
        continue
    else:
        encoding[data_var] = {'_FillValue': None}


outputfile = os.path.join(os.environ['PROJDIR'], 'E3SM', 'inputdata',
                          'atm', 'datm7', 'CLM1PT_data', 'SPRUCE_data',
                          'clm_params.nc_eca_yang_dmr_02242021')
hr.to_netcdf(outputfile, encoding = encoding, format = 'NETCDF3_CLASSIC')


hr.close()
h_in.close()
