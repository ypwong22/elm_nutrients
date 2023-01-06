""" test if this works """
import pandas as pd
import numpy as np
import os


def getvar(fname, varname):
    from netCDF4 import Dataset
    nffile = Dataset(fname,"r")
    if varname in nffile.variables:
      varvals = nffile.variables[varname][:]
    else:
      # print('Warning: '+varname+' not in '+fname)
      raise ValueError('"%s" not in %s'%(varname,fname))
      varvals=[-1]
    nffile.close()
    return varvals

def putvar(fname, varname, varvals):
    from netCDF4 import Dataset
    import numpy as np
    nffile = Dataset(fname,"a")
    if (varname in nffile.variables):
      nffile.variables[varname][...] = varvals
    else:
      print('Warning: '+varname+' not in '+fname)
    nffile.close()
    ierr = 0
    return ierr


def postproc(myvars, myyear_start, myyear_end, myday_start, myday_end, myavg, \
             myfactor, myoffset, mypft, mytreatment, data):
    index=0
    ierr = 0
    thiscol = 0

    # dummy
    case = '20211008_plot04_US-SPR_ICB20TRCNPRDCTCBC'
    rundir = os.path.join(os.environ['PROJDIR'], 'E3SM/output/20211008_plot04_US-SPR_ICB20TRCNPRDCTCBC/run/')
    model_name = 'clm2'

    index=0
    ierr = 0
    thiscol = 0
    for v in myvars:
        if 'US-SPR' in case:
          # split the pft and multiplication factors
          pft_list = [int(t) for t in mypft[index].split(',')]
          factor_list = [float(f) for f in myfactor[index].split(',')]

        ndays_total = 0
        output = []
        n_years = myyear_end[index]-myyear_start[index]+1
        npy=1
        for y in range(myyear_start[index],myyear_end[index]+1):
            if (pft_list[0] <= 0 or 'PFT' in v):
              fname = rundir+case+'.'+model_name+'.h0.'+str(10000+y)[1:]+'-01-01-00000.nc'
              myindex_list = [max(0, pft_list[0])]
              hol_add = 1
            else:
              fname = rundir+case+'.'+model_name+'.h1.'+str(10000+y)[1:]+'-01-01-00000.nc'
              myindex_list = pft_list
              hol_add = 17
            if (os.path.exists(fname)):
              mydata = getvar(fname,v) 
              if ('ZWT' in v):
                mydata2 = getvar(fname,'H2OSFC')
              if (len(mydata) < 10):
                npy = 1 
              elif (len(mydata) >= 365):    #does not currently allow hourly
                npy = 365
            else:
              #print(fname)
              mydata = np.zeros([npy,34], float)+np.NaN
            #get output and average over days/years
            n_days = myday_end[index]-myday_start[index]+1
            ndays_total = ndays_total + n_days
            #get number of timesteps per output file
            #print(v, n_days, ndays_total)

            if (npy == 365):
                for d in range(myday_start[index]-1,myday_end[index]):
                    if ('US-SPR' in case and 'ZWT' in v):
                      #Use hollows for water table height
                      output.append(mydata[d][myindex_list[0]+hol_add]*factor_list[0] \
                             +myoffset[index]+mydata2[d][myindex_list[0]+hol_add]/1000.)
                    elif ('US-SPR' in case):
                      temp = 0.
                      for m, myindex in enumerate(myindex_list):
                        temp = temp + (0.36 * mydata[d][myindex] + 0.64 * mydata[d][myindex+hol_add]) * factor_list[m]
                      temp = temp + myoffset[index]
                      output.append(temp)
                    else:
                      output.append(mydata[d][myindex]*myfactor[index] + myoffset[index])
            elif (npy == 1):                    #Assume annual output (ignore days)
               for d in range(myday_start[index]-1,myday_end[index]):    #28-38 was myindex
                 if ('SCPF' in v):
                   output.append(sum(mydata[0,28:38])/10.0*myfactor[index]+myoffset[index])
                 elif ('NPLANT_SCLS' in v):
                   output.append(sum(mydata[0,1:])*myfactor[index]+myoffset[index])
                 elif ('SCLS' in v):
                    output.append(sum(mydata[0,:])*myfactor[index]+myoffset[index])
                 else:
                   try:
                     output.append(mydata[0,myindex]*myfactor[index]+myoffset[index])
                   except:
                     output.append(np.NaN)
        for i in range(0,int(ndays_total/myavg[index])):
            data[thiscol] = sum(output[(i*myavg[index]):((i+1)*myavg[index])])/myavg[index]
            thiscol=thiscol+1
        index=index+1


if __name__ == '__main__':

    postproc_file = os.path.join(os.environ['PROJDIR'], 'Phenology_ELM', 'output_analysis', 'postproc_vars_SPRUCE')

    myvars=[]
    myyear_start=[]
    myyear_end=[]
    myday_start=[]
    myday_end=[]
    myavg_pd=[]
    myfactor=[]
    myoffset=[]
    mypft=[]
    myobs=[]
    myobs_err=[]
    mytreatment=[]

    postproc_input = open(postproc_file,'r')
    data_cols = 0
    for s in postproc_input:
        if (s[0:1] != '#'):
            myvars.append(s.split()[0])
            myyear_start.append(int(s.split()[1]))
            myyear_end.append(int(s.split()[2]))
            myday_start.append(int(s.split()[3]))
            myday_end.append(int(s.split()[4]))
            myavg_pd.append(int(s.split()[5]))
            myfactor.append(s.split()[6])
            myoffset.append(float(s.split()[7]))
            if (len(s.split()) >= 9):
                mypft.append(s.split()[8])
            else:
                mypft.append(-1)
            if (len(s.split()) >= 11):
                myobs.append(float(s.split()[9]))
                myobs_err.append(float(s.split()[10]))
            else: 
                myobs.append(-9999)
                myobs_err.append(-9999)
            if (len(s.split()) == 12):        
                mytreatment.append(s.split()[11])     
            else:
                mytreatment.append('NA')
            days_total = (int(s.split()[2]) - int(s.split()[1])+1)*(int(s.split()[4]) - int(s.split()[3])+1)        
            data_cols = int(round(data_cols + days_total / int(s.split()[5])))
            print('DATA_COLS',data_cols)
    print(mytreatment)
    data_row = np.zeros([data_cols], float)-999
    postproc_input.close()

    postproc(myvars, myyear_start, myyear_end, myday_start, myday_end, myavg_pd, \
             myfactor, myoffset, mypft, mytreatment, data_row)

    data = pd.DataFrame([list(data_row), myobs, myobs_err], index = ['data', 'myobs', 'myobs_err']).T
    data.to_csv(os.path.join(os.environ['PROJDIR'], 'Phenology_ELM', 'output_analysis', 'test_postproc.csv'))