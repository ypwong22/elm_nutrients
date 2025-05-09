from netCDF4 import Dataset
import numpy as np
import os

casename = 'UQ_20231118_US-SPR_ICB20TRCNPRDCTCBC'
ambcase=os.path.join(os.environ['E3SM_ROOT'], 'output', 'UQ', 
                     'UQ_20231118_US-SPR_ICB20TRCNPRDCTCBC', 'g03067')
cases = ['TAMB', 'T0.00','T2.25','T4.50','T6.75','T9.00', \
         'T0.00CO2','T2.25CO2','T4.50CO2','T6.75CO2','T9.00CO2']

parmfile = Dataset(ambcase+'/clm_params_03067.nc', 'r') #Dataset(ambcase+'/run/clm_params.nc','r')
surffile = Dataset(ambcase+'/surfdata_03067.nc', 'r') #Dataset(ambcase+'/run/surfdata.nc','r')

#Hummock/hollow parameters
hum_frac = parmfile['hum_frac'][0]
humhol_ht = parmfile['humhol_ht'][0]

ystart = 2016
yend = 2021 # 2023
myvars = ['GPP','NEE','NPP','HR','NPP_pft2','NPP_pft11','NPP_pft12','TOTVEGC_ABG_pft2', 'TOTVEGC_ABG_pft11', \
        'TOTVEGC_ABG_pft12','TOTSOMC','TOTLITC','TLAI_pft2','TLAI_pft11','TLAI_pft12',\
        'CPOOL', #,'NFIX_TO_SMINN','POTENTIAL_IMMOB',\
       'ACTUAL_IMMOB','FPG','FPI','FPG_P','FPI_P'] # 'FCH4','GROSS_NMIN','PLANT_NDEMAND','SMINN_TO_PLANT',
ncol = 17 #32

for v in myvars:
  fac = 24*3600*365
  if 'CH4' in v:
      fac=fac*1000
  if 'TOT' in v or 'TLAI' in v or 'LEAFC' in v or 'FP' in v or 'POOL' in v:
      fac = 1
  var_annual = np.zeros([yend-ystart+1,len(cases)], float)
  cnum=0
  for c in cases:
    thiscase = ambcase + '/' + c # ambcase.replace('TAMB',c)
    lnd_in = open(ambcase+'/lnd_in','r')
    #Get dynamic PFT info
    for line in lnd_in:
        if ('landuse_timeseries' in line):
            lufile = line.split('=')[1][:-1].replace("'","")
    lnd_in.close()
    pftdyn = Dataset(lufile,'r')
    #Process model output
    for y in range(ystart,yend+1):
        pct_pft = pftdyn['PCT_NAT_PFT'][min(y-1850,170),:,0]
        htype = 'h1'
        if ('_pft' in v):
            htype = 'h2'
        myfile = thiscase+'/'+casename+'.elm.'+htype+'.'+str(y)+'-01-01-00000.nc'
        mydata = Dataset(myfile,'r')
        if '_pft' in v:
            thispft = int(v.split('_pft')[1])
            myvar = mydata[v.split('_pft')[0]][:,thispft]*hum_frac*pct_pft[thispft]/100. + \
                    mydata[v.split('_pft')[0]][:,thispft+ncol]*(1.0-hum_frac)*pct_pft[thispft]/100.
            if '_pft2' in v:
                #add larch, both tree types
                thispft = thispft+1
                myvar = myvar+ mydata[v.split('_pft')[0]][:,thispft]*hum_frac*pct_pft[thispft]/100. + \
                    mydata[v.split('_pft')[0]][:,thispft+ncol]*(1.0-hum_frac)*pct_pft[thispft]/100.
        else:
            myvar = mydata[v][:,0]*hum_frac + mydata[v][:,1]*(1.0-hum_frac)
        var_annual[y-ystart,cnum] = np.mean(myvar)*fac
    cnum=cnum+1
    pftdyn.close()
  x = [0,2.25,4.5,6.75,9]    
  slope, intercept = np.polyfit(x, np.mean(var_annual[:,0:5],axis=0), 1)
  print(v, np.mean(var_annual[:,0]), np.mean(var_annual[:,4]), slope, \
          np.mean(var_annual[:,5]), \
          np.mean(var_annual[:,9]))
