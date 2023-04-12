import numpy as np
import os
import pandas as pd
from datetime import datetime, timedelta
from .analysis import *
from scipy.optimize import curve_fit


lat, lon, tzn, elv = '47.563 N', '93.438 W', 'US/Central', 430


def _find_solstice(year_list):
    # find last winter solstice
    solstice  = pd.Series(index = year_list)
    for year in year_list: 
        last_year = pd.DatetimeIndex([datetime(year-1, 12, 1) + timedelta(days = i) \
            for i in range((datetime(year,1,1)-datetime(year-1, 1, 1)).days)])

        dl_prev = 1e20
        for ind, yy, mm, dd in zip(range(len(last_year)), 
                                        last_year.year, last_year.month, last_year.day):
            dl = daylength(yy, mm, dd, lat, lon, elv, tzn)
            if dl >= dl_prev:
                break
            else:
                dl_prev = dl
        solstice[year] = last_year[ind]
    return solstice


class GDD():
    """ Vanilla growing degree days model implemented in the original ELM 

    dayofyear : pandas.DataFrame
        Start of season dates. The index is year. The columns are different chambers. 

    tair      : pandas.DataFrame
        Annual mean air temperature in degree celsius. The index is date time. The columns are different chambers. 

    tsoi      : pandas.DataFrame
        Soil temperature in degrees celsius. The index is date time. The columns are different chambers. 
    """
    def __init__(self):
        self.params    = None

    def _crit_gdd(self, x, a, b):
        # crit_gdd  = pd.DataFrame(np.exp(4.8 + 0.13 * annt2m.values), index = annt2m.index + 1)
        return np.exp(a + b * x)

    def fit(self, dayofyear, tair, tsoi):
        self.dayofyear = dayofyear
        self.tair      = tair
        self.tsoi      = tsoi
        self.solstice = _find_solstice(self.dayofyear.index)

        print(self.solstice)

        gdds = pd.DataFrame(np.nan, index = self.dayofyear.index, columns = self.dayofyear.columns)
        for year in self.dayofyear.index:
            for c in self.dayofyear.columns:
                tsoi_accu = self.tsoi.loc[(self.tsoi.index >= self.solstice.loc[year]) & \
                                          (self.tsoi.index <= (datetime(year-1, 12, 31) \
                                                               + timedelta(days = float(self.dayofyear.loc[year, c])))), c]
                tsoi_accu = (tsoi_accu * (tsoi_accu > 0.)).sum(axis = 0) * \
                            (self.tsoi.index[1] - self.tsoi.index[0]).total_seconds() / 86400
                gdds.loc[year, c] = tsoi_accu

        print(gdds)

        x = self.tair.loc[gdds.index - 1, :].values.reshape(-1)
        y = gdds.values.reshape(-1)
        popt, pcov = curve_fit(self._crit_gdd, x, y)
        ypred = self._crit_gdd(x, *popt)
        self.params = popt

        return popt, ypred, x, y

    def predict(self, x):
        return self._crit_gdd(x, *self.params)


class ALT():
    """ Alternate model 

    Fit parameters

    dayofyear : pandas.DataFrame
        Start of season dates. The index is year. The columns are different chambers. 

    tsoi      : pandas.DataFrame
        Soil temperature in degrees celsius. The index is date time. The columns are different chambers. 
    """
    def __init__(self):
        # a, b, c, t_base
        self.params    = None

    def _leafout(self, tsoi, a, b, c, t_base):
        """ tsoi is a dataframe """
        year_list = tsoi.index.year.unique()
        solstice  = _find_solstice(year_list)

        if tsoi.index[0] > solstice[0]:
            year_list = year_list[1:]

        leafout_dates = pd.DataFrame(np.nan, index = year_list, columns = tsoi.columns)
        for year in year_list:
            for c in tsoi.columns:
                tsoi_accu = tsoi.loc[(tsoi.index >= solstice.loc[year]) & \
                                     (tsoi.index <= datetime(year, 6, 30)), c]

                gdd_accu = (tsoi_accu * (tsoi_accu > t_base)).cumsum(axis = 0) * \
                           (tsoi.index[1] - tsoi.index[0]).total_seconds() / 86400
                cdd_accu = (tsoi_accu * (tsoi_accu < t_base)).cumsum(axis = 0) * \
                           (tsoi.index[1] - tsoi.index[0]).total_seconds() / 86400

                diff = (gdd_accu - a - b * np.exp(c * cdd_accu)).values
                ind = np.where((diff[1:] >= 0) & (diff[:-1] < 0))[0][0] + 1
                if tsoi_accu.index[ind].year == year:
                    leafout_dates.loc[year, c] = tsoi_accu.index[ind].dayofyear
                elif tsoi_accu.index[ind].year == (year - 1):
                    leafout_dates.loc[year, c] = tsoi_accu.index[ind].dayofyear - 365
                else:
                    raise 'Wrong year'
        return leafout_dates.values.reshape(-1)


    def _leafout_fit(self, x, a, b, c, t_base):
        # Growth happens when Sf(t) >= a + b * exp(c * Sc(t))
        # Obtain the date when the above equation is satisfied
        leafout_dates = pd.DataFrame(np.nan, index = self.dayofyear.index, 
                                     columns = self.dayofyear.columns)
        for year in self.dayofyear.index:
            for c in self.dayofyear.columns:
                tsoi_accu = self.tsoi.loc[(self.tsoi.index >= self.solstice.loc[year]) & \
                                          (self.tsoi.index <= datetime(year, 6, 30)), c]

                gdd_accu = (tsoi_accu * (tsoi_accu > t_base)).cumsum(axis = 0) * \
                           (self.tsoi.index[1] - self.tsoi.index[0]).total_seconds() / 86400
                cdd_accu = (tsoi_accu * (tsoi_accu < t_base)).cumsum(axis = 0) * \
                           (self.tsoi.index[1] - self.tsoi.index[0]).total_seconds() / 86400

                diff = (gdd_accu - a - b * np.exp(c * cdd_accu)).values
                ind = np.where((diff[1:] >= 0) & (diff[:-1] < 0))[0][0] + 1
                if tsoi_accu.index[ind].year == year:
                    leafout_dates.loc[year, c] = tsoi_accu.index[ind].dayofyear
                elif tsoi_accu.index[ind].year == (year - 1):
                    leafout_dates.loc[year, c] = tsoi_accu.index[ind].dayofyear - 365
                else:
                    raise 'Wrong year'
        return leafout_dates.values.reshape(-1)


    def fit(self, dayofyear, tsoi):
        self.dayofyear = dayofyear
        self.tsoi      = tsoi
        self.solstice  = _find_solstice(self.dayofyear.index)

        y = self.dayofyear.values.reshape(-1)
        popt, pcov = curve_fit(self._leafout_wrap, 0, y)
        ypred = self._leafout_wrap(*popt)
        self.params = popt

        return popt, ypred, y

    def predict(self, x):
        return self._leafout(x, *self.params)


class SWC():
    """ soil moisture based """
    pass
