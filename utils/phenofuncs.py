import numpy as np
import os
import pandas as pd
from datetime import datetime, timedelta
from .analysis import *
from scipy.optimize import minimize, basinhopping, differential_evolution
import itertools as it


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
        Annual mean air temperature in degree celsius. The index is year. The columns are different chambers. 

    tsoi      : pandas.DataFrame
        Soil temperature in degrees celsius. The index is date time. The columns are different chambers. 
    """
    def __init__(self, dayofyear, tsoi, tair):
        self.params    = None
        self.dayofyear = dayofyear # day of year
        self.tair      = tair # annual air temperature
        self.tsoi      = tsoi # daily soil temperature

    def _budbreak(self, a, b):
        year_list = np.unique(self.tsoi.index.year)
        solstice = _find_solstice(year_list)
        budbreak = pd.DataFrame(np.nan, index = year_list, columns = self.tsoi.columns)
        for year in year_list:
            crit_gdd = np.exp(a + b * self.tair.loc[year - 1, :])
            temp = self.tsoi.loc[(self.tsoi.index >= solstice.loc[year]) & (self.tsoi.index.year <= year), :]
            for c in self.tsoi.columns:
                temp_temp = np.where(np.cumsum(temp[c] * (temp[c] > 0.)) >= crit_gdd.loc[c])[0]
                if len(temp_temp) == 0:
                    budbreak.loc[year, c] = 999
                else:
                    temp_temp = temp_temp[0]
                    budbreak.loc[year, c] = temp.index[temp_temp].dayofyear
        return budbreak

    def _calc(self, popt):
        # varied amount by python is too small
        temp = self.dayofyear - self._budbreak(*popt)
        print(*popt, np.sqrt(np.mean(np.power(temp.values.reshape(-1), 2))))
        return np.sqrt(np.mean(np.power(temp.values.reshape(-1), 2)))

    def fit(self):
        res = differential_evolution(self._calc, [(2, 8), (0, 0.1)])
        params = res.x
        ypred = self._budbreak(*params)
        temp = self.dayofyear - ypred
        rmse = np.sqrt(np.mean(np.power(temp.values.reshape(-1), 2)))
        return params, rmse, ypred