# https://gist.github.com/anttilipp/ed3ab35258c7636d87de6499475301ce
import numpy as np
import pandas as pd
import os


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