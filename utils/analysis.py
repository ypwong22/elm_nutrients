# https://gist.github.com/anttilipp/ed3ab35258c7636d87de6499475301ce
import numpy as np
def daylength(dayOfYear, lat):
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