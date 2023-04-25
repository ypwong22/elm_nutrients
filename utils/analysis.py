# https://gist.github.com/anttilipp/ed3ab35258c7636d87de6499475301ce
import numpy as np
import pandas as pd
import os
from glob import glob
import xarray as xr
from .constants import *


def get_treatment_string(chamber):
    chamber_levels = {'07': ['TAMB',  0], '14': ['TAMB',  0], '21': ['TAMB',  0],
                      '06': [0     ,  0], '19': [0     ,500], 
                      '11': [2.25  ,500], '20': [2.25  ,  0],
                      '04': [4.5   ,500], '13': [4.5   ,  0],
                      '08': [6.75  ,  0], '16': [6.75  ,500],
                      '10': [9     ,500], '17': [9     ,  0]}
    temperature, co2 = chamber_levels[f'{chamber:02g}']
    if temperature == 'TAMB':
        treatment = 'TAMB'
    else:
        treatment = f'T{temperature:.2f}'
    if co2 > 0:
        treatment = treatment + 'CO2'
    return treatment


def get_mossfrac(year, treatment):
    mossfrac = pd.read_excel(os.path.join(os.environ['HOME'], 'Git', 'phenology_elm', 'Sphagnum_fraction.xlsx'), index_col = 0, skiprows = 1,
                             engine = 'openpyxl').drop(['plot','Temp','CO2'], axis = 1)
    mossfrac[2015] = mossfrac[2016]
    return mossfrac.loc[treatment, year]


def read_mortality():
    """ Return the data from SPRUCE_S1_Minirhizotron data_2012_For Yaoping.xls """
    dates = pd.DatetimeIndex([datetime(2012, 5,16), datetime(2012, 5,21), datetime(2012, 5,28),
                              datetime(2012, 6, 5), datetime(2012, 6,12), datetime(2012, 6,18),
                              datetime(2012, 7, 3), datetime(2012, 7, 9), datetime(2012, 7,17),
                              datetime(2012, 7,24), datetime(2012, 7,30), datetime(2012, 8,17),
                              datetime(2012, 9, 1), datetime(2012, 9,15)])
    # 
    mortality = np.array([0, 33.4986, 104.1712, 27.1749, 45.6811, 30.8191, 93.8463, 49.2897, 61.8272, 69.4396, 50.9765, 158.3639, 176.7845, 99.4316])

    days = (dates[1:] - dates[:-1]).days

    mortality = mortality[1:] / days

    date_list = []
    mort_list = []
    for i in range(len(days)):
        date_list.extend(list(pd.date_range(dates[i], dates[i+1] - timedelta(days = 1))))
        mort_list.extend([mortality[i]]*days[i])

    mort = pd.Series(mort_list, index = date_list)
    mort_eos = mort.index[np.where(mort.cumsum() >= mort.cumsum()[-1] * 0.75)[0]][0]

    return mort, mort_eos


def daylength_simple(dayOfYear, lat):
    """
    Computes the length of the day (the time between sunrise and
    sunset) given the day of the year and latitude of the location.
    Function uses the Brock model for the computations.
    For more information see, for example,
    Forsythe et al., "A model comparison for daylength as a
    function of latitude and day of year", Ecological Modelling,
    1995.
    Parameters
    ----------
    dayOfYear : int
        The day of the year. 1 corresponds to 1st of January
        and 365 to 31st December (on a non-leap year).
    lat : float
        Latitude of the location in degrees. Positive values
        for north and negative for south.
    Returns
    -------
    d : float
        Daylength in hours.
    """
    latInRad = np.deg2rad(lat)
    declinationOfEarth = 23.45*np.sin(np.deg2rad(360.0*(283.0+dayOfYear)/365.0))
    if -np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth)) <= -1.0:
        return 24.0
    elif -np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth)) >= 1.0:
        return 0.0
    else:
        hourAngle = np.rad2deg(np.arccos(-np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth))))
        return 2.0*hourAngle/15.0


# https://stackoverflow.com/questions/54485777/finding-twilight-times-with-skyfield
from skyfield import api, almanac
from datetime import datetime, timedelta
import pytz
from skyfield.nutationlib import iau2000b

def daylength(year, month, day, lat, lon, elv, tzn):
    """Build a function of time that returns the daylength.

    The function that this returns will expect a single argument that is a 
    :class:`~skyfield.timelib.Time` and will return ``True`` if the sun is up
    or twilight has started, else ``False``.
    """
    DAYLENGTH_CENTER_HORIZON = 0.0
    DAYLENGTH_TOP_HORIZON = 0.26667
    DAYLENGTH_TOP_HORIZON_APPARENTLY = 0.8333
    DAYLENGTH_CIVIL_TWILIGHT = 6.0
    DAYLENGTH_NAUTICAL_TWILIGHT = 12.0
    DAYLENGTH_ASTRONOMICAL_TWILIGHT = 18.0

    planets = api.load('de421.bsp')
    topos   = api.Topos(lat, lon, elevation_m=elv)
    ts = api.load.timescale()
    tz = pytz.timezone(tzn)

    sun = planets['sun']
    topos_at = (planets['earth'] + topos).at

    t0 = ts.utc(datetime(year, month, day, tzinfo = tz))
    t1 = ts.utc(tz.normalize(datetime(year, month, day, tzinfo = tz) + timedelta(1)))

    def is_sun_up_at(t):
        """Return `True` if the sun has risen by time `t`."""
        t._nutation_angles = iau2000b(t.tt)
        return topos_at(t).observe(sun).apparent().altaz()[0].degrees > -DAYLENGTH_TOP_HORIZON_APPARENTLY

    is_sun_up_at.rough_period = 0.5  # twice a day

    center_time, _ = almanac.find_discrete(t0, t1, is_sun_up_at)
    up, down       = center_time.utc_iso()
    daylength = datetime.strptime(down, '%Y-%m-%dT%H:%M:%SZ') - \
                datetime.strptime(up, '%Y-%m-%dT%H:%M:%SZ') 
    return daylength.total_seconds()


def kge(simulations, evaluation):
    """
    Stolen from the hydroeval package
    
    Original Kling-Gupta Efficiency (KGE) and its three components
    (r, α, β) as per `Gupta et al., 2009
    <https://doi.org/10.1016/j.jhydrol.2009.08.003>`_.

    Note, all four values KGE, r, α, β are returned, in this order.

    :Calculation Details:
        .. math::
           E_{\\text{KGE}} = 1 - \\sqrt{[r - 1]^2 + [\\alpha - 1]^2
           + [\\beta - 1]^2}
        .. math::
           r = \\frac{\\text{cov}(e, s)}{\\sigma({e}) \\cdot \\sigma(s)}
        .. math::
           \\alpha = \\frac{\\sigma(s)}{\\sigma(e)}
        .. math::
           \\beta = \\frac{\\mu(s)}{\\mu(e)}

        where *e* is the *evaluation* series, *s* is (one of) the
        *simulations* series, *cov* is the covariance, *σ* is the
        standard deviation, and *μ* is the arithmetic mean.

    """
    # calculate error in timing and dynamics r
    # (Pearson's correlation coefficient)
    sim_mean = np.mean(simulations, axis=0, dtype=np.float64)
    obs_mean = np.mean(evaluation, dtype=np.float64)

    r_num = np.sum((simulations - sim_mean) * (evaluation - obs_mean),
                   axis=0, dtype=np.float64)
    r_den = np.sqrt(np.sum((simulations - sim_mean) ** 2,
                           axis=0, dtype=np.float64)
                    * np.sum((evaluation - obs_mean) ** 2,
                             dtype=np.float64))
    r = r_num / r_den
    # calculate error in spread of flow alpha
    alpha = np.std(simulations, axis=0) / np.std(evaluation, dtype=np.float64)
    # calculate error in volume beta (bias of mean discharge)
    beta = (np.sum(simulations, axis=0, dtype=np.float64)
            / np.sum(evaluation, dtype=np.float64))
    # calculate the Kling-Gupta Efficiency KGE
    kge_ = 1 - np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2)

    return np.vstack((kge_, r, alpha, beta))


def extract_sim_one(prefix, var_list = {'pft': [], 'col': []}):
    tvec = pd.date_range('2015-01-01', '2020-12-31', freq = '1D')
    tvec = tvec[(tvec.month != 2) | (tvec.day != 29)]

    pft_list = [2, 3, 11, 12]
    hol_add = 17

    var_tuples = [(i, j) for i in var_list['pft'] for j in pft_list]
    var_tuples = var_tuples + [(i, 0) for i in var_list['col']]

    colnames = pd.MultiIndex.from_tuples([(i, j[0], j[1], k) for i in chamber_list_complete for j in var_tuples for k in ['hummock', 'hollow']])
    collection = pd.DataFrame(np.nan, index = tvec, columns = colnames)

    for plot in chamber_list_complete:
        path_data = os.path.join(os.environ['PROJDIR'], 'E3SM', 'output', f'{prefix}_plot{plot:02d}_US-SPR_ICB20TRCNPRDCTCBC', 'run')

        var_list_pft = var_list['pft']
        flist = sorted(glob(os.path.join(path_data, '*.h2.*.nc')))[:-1]
        hr = xr.open_mfdataset(flist, decode_times = False)
        for var in var_list_pft:
            for pft in pft_list:
                collection.loc[:, (plot, var, pft, 'hummock')] = hr[var][:, pft].values
                collection.loc[:, (plot, var, pft, 'hollow')] = hr[var][:, pft + hol_add].values
        hr.close()

        var_list_col  = var_list['pft']
        flist = sorted(glob(os.path.join(path_data, '*.h1.*.nc')))[:-1]
        hr = xr.open_mfdataset(flist, decode_times = False)
        for var in var_list_col:
            if var == 'TSOI_3':
                pass
            elif var == 'TSOI_AVG':
                pass
            elif var == 'SWC_3':
                pass
            elif var == 'SWC_AVG':
                pass
            elif var == 'SOIPSI_3':
                pass
            elif var == 'SOIPSI_AVG':
                pass
            else:
                collection.loc[:, (plot, var, 0, 'hummock')] = hr[var][:, 0].values
                collection.loc[:, (plot, var, 0, 'hollow')] = hr[var][:, 1].values
        hr.close()

    return collection