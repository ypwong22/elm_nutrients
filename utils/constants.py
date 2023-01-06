import pandas as pd

chamber_levels = {'04': [4.5,500], '06': [0   ,  0], '08': [6.75,  0], 
                  '10': [9  ,500], '11': [2.25,500], '13': [4.5 ,  0], '16': [6.75,500], 
                  '17': [9  ,  0], '19': [0   ,500], '20': [2.25,  0],
                  '07': [0  ,  0], '21': [0  ,  0]}

# T: 0, 2.25, 4.5, 6.75, 9; first ambient, then elevated CO2
chamber_list = ['06', '20', '13', '08', '17', '19', '11', '04', '16', '10']

sims_prefix = ['20211008', '20220103',
               '20220407_stopdate171_switch1', '20220407_stopdate365_switch1',
               '20220407_stopdate171_switch2', '20220407_stopdate365_switch2',
               '20220407_stopdate171_switch3', '20220407_stopdate365_switch3',
               '20220407_stopdate171_switch4', '20220407_stopdate365_switch4']
sims_names =  ['Original', 'AltEvgrPheno',
               'AltRoot_171_1', 'AltRoot_365_1',
               'AltRoot_171_2', 'AltRoot_365_2',
               'AltRoot_171_3', 'AltRoot_365_3',
               'AltRoot_171_4', 'AltRoot_365_4']
sims_colors = ['#08519c', '#54278f', # blue, purple
               '#a50f15', '#006d2c',
               '#de2d26', '#31a354',
               '#fb6a4a', '#74c476', 
               '#fc9272', '#a1d99b'] # red - 171, green - 365

sim_tvec = pd.date_range('2015-01-01', '2020-12-31', freq = '1D')
sim_tvec = sim_tvec[(sim_tvec.month != 2) | (sim_tvec.day != 29)]

soil_interfaces = [0, 0.0175, 0.0451, 0.0906, 0.1655, 0.2891, 0.4929, 0.8289, 1.3828,
                   2.2961, 3.8019, 6.2845, 10.3775, 17.1259, 28.2520, 42.1032]

# AR     : autotrophic respiration
# FSDS   : atmospheric incident solar radiation
# FSA    : absorbed solar radiation
# QVEGT  : canopy transpiration
# FSH    : sensible heat
# TBOT   : atmospheric air temperature
# H2OSOI : vertically resolved soil water; limit to 40 cm
# CH4PROD: gridcell total production of CH4
var_list = ['TLAI_lala'   , 'TLAI_pima'   , 'TLAI_shrub'    , 'TLAI_moss'   , 'TLAI'   ,
            'QVEGT_lala'  , 'QVEGT_pima'  , 'QVEGT_shrub'   , 'QVEGT_moss'  , 'QVEGT'  ,
            'GPP_lala'    , 'GPP_pima'    , 'GPP_shrub'     , 'GPP_moss'    , 'GPP'    ,
            'NPP_lala'    , 'NPP_pima'    , 'NPP_shrub'     , 'NPP_moss'    , 'NPP'    ,
            'AGNPP_lala'  , 'AGNPP_pima'  , 'AGNPP_shrub'   , 'AGNPP_moss'  , 'AGNPP'  ,
            'BGNPP_lala'  , 'BGNPP_pima'  , 'BGNPP_shrub'   , 'BGNPP_moss'  , 'BGNPP'  ,
            'FROOTC_lala' , 'FROOTC_pima' , 'FROOTC_shrub'  , 'FROOTC_moss' , 'FROOTC' ,
            'AR_lala'     , 'AR_pima'     , 'AR_shrub'      , 'AR_moss'     , 'AR'     ,
            'FSDS', 'FSA', 'FSH', 'TBOT', 'NEE', 'HR',
            'H2OSOI', 'QOVER', 'ZWT'] # , 'CH4PROD', 'LITFALL_lala', 'LITFALL_pima', 'LITFALL_shrub' , 'LITFALL_moss', 'LITFALL',
unit_list = ['none'        , 'none'        , 'none'      , 'none'      , 'none'      ,
             'mm day-1'    , 'mm day-1'    , 'mm day-1'  , 'mm day-1'  , 'mm day-1'  ,
             'gC m-2 day-1', 'gC m-2 day-1', 'gC m-2 day-1', 'gC m-2 day-1', 'gC m-2 day-1',
             'gC m-2 day-1', 'gC m-2 day-1', 'gC m-2 day-1', 'gC m-2 day-1', 'gC m-2 day-1',
             'gC m-2 day-1', 'gC m-2 day-1', 'gC m-2 day-1', 'gC m-2 day-1', 'gC m-2 day-1',
             'gC m-2 day-1', 'gC m-2 day-1', 'gC m-2 day-1', 'gC m-2 day-1', 'gC m-2 day-1',
             'gC m-2'      , 'gC m-2'      , 'gC m-2'    , 'gC m-2'    , 'gC m-2'    ,
             'gC m-2 day-1','gC m-2 day-1','gC m-2 day-1','gC m-2 day-1','gC m-2 day-1',
             'W m-2', 'W m-2', 'W m-2', 'degC', 'gC m-2 day-1', 'gC m-2 day-1',
             'mm3 mm-3', 'mm day-1', 'mm'] # 'gC m-2 day-1', 'gC m-2 day-1','gC m-2 day-1','gC m-2 day-1','gC m-2 day-1','gC m-2 day-1',