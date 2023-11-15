# https://gist.github.com/anttilipp/ed3ab35258c7636d87de6499475301ce
import numpy as np
import pandas as pd
import os
from glob import glob
import xarray as xr
from .constants import *
from .paths import *


def get_treatment_string(chamber):
    chamber_levels = {
        "07": ["TAMB", 0],
        "14": ["TAMB", 0],
        "21": ["TAMB", 0],
        "06": [0, 0],
        "19": [0, 500],
        "11": [2.25, 500],
        "20": [2.25, 0],
        "04": [4.5, 500],
        "13": [4.5, 0],
        "08": [6.75, 0],
        "16": [6.75, 500],
        "10": [9, 500],
        "17": [9, 0],
    }
    temperature, co2 = chamber_levels[f"{chamber:02g}"]
    if temperature == "TAMB":
        treatment = "TAMB"
    else:
        treatment = f"T{temperature:.2f}"
    if co2 > 0:
        treatment = treatment + "CO2"
    return treatment


def get_mossfrac(year, treatment):
    mossfrac = pd.read_excel(
        os.path.join(
            os.environ["HOME"], "Git", "phenology_elm", "Sphagnum_fraction.xlsx"
        ),
        index_col=0,
        skiprows=1,
        engine="openpyxl",
    ).drop(["plot", "Temp", "CO2"], axis=1)
    mossfrac[2015] = mossfrac[2016]
    return mossfrac.loc[treatment, year]


def read_mortality():
    """Return the data from SPRUCE_S1_Minirhizotron data_2012_For Yaoping.xls"""
    dates = pd.DatetimeIndex(
        [
            datetime(2012, 5, 16),
            datetime(2012, 5, 21),
            datetime(2012, 5, 28),
            datetime(2012, 6, 5),
            datetime(2012, 6, 12),
            datetime(2012, 6, 18),
            datetime(2012, 7, 3),
            datetime(2012, 7, 9),
            datetime(2012, 7, 17),
            datetime(2012, 7, 24),
            datetime(2012, 7, 30),
            datetime(2012, 8, 17),
            datetime(2012, 9, 1),
            datetime(2012, 9, 15),
        ]
    )
    #
    mortality = np.array(
        [
            0,
            33.4986,
            104.1712,
            27.1749,
            45.6811,
            30.8191,
            93.8463,
            49.2897,
            61.8272,
            69.4396,
            50.9765,
            158.3639,
            176.7845,
            99.4316,
        ]
    )

    days = (dates[1:] - dates[:-1]).days

    mortality = mortality[1:] / days

    date_list = []
    mort_list = []
    for i in range(len(days)):
        date_list.extend(
            list(pd.date_range(dates[i], dates[i + 1] - timedelta(days=1)))
        )
        mort_list.extend([mortality[i]] * days[i])

    mort = pd.Series(mort_list, index=date_list)
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
    declinationOfEarth = 23.45 * np.sin(np.deg2rad(360.0 * (283.0 + dayOfYear) / 365.0))
    if -np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth)) <= -1.0:
        return 24.0
    elif -np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth)) >= 1.0:
        return 0.0
    else:
        hourAngle = np.rad2deg(
            np.arccos(-np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth)))
        )
        return 2.0 * hourAngle / 15.0


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

    planets = api.load("de421.bsp")
    topos = api.Topos(lat, lon, elevation_m=elv)
    ts = api.load.timescale()
    tz = pytz.timezone(tzn)

    sun = planets["sun"]
    topos_at = (planets["earth"] + topos).at

    t0 = ts.utc(datetime(year, month, day, tzinfo=tz))
    t1 = ts.utc(tz.normalize(datetime(year, month, day, tzinfo=tz) + timedelta(1)))

    def is_sun_up_at(t):
        """Return `True` if the sun has risen by time `t`."""
        t._nutation_angles = iau2000b(t.tt)
        return (
            topos_at(t).observe(sun).apparent().altaz()[0].degrees
            > -DAYLENGTH_TOP_HORIZON_APPARENTLY
        )

    is_sun_up_at.rough_period = 0.5  # twice a day

    center_time, _ = almanac.find_discrete(t0, t1, is_sun_up_at)
    up, down = center_time.utc_iso()
    daylength = datetime.strptime(down, "%Y-%m-%dT%H:%M:%SZ") - datetime.strptime(
        up, "%Y-%m-%dT%H:%M:%SZ"
    )
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

    r_num = np.sum(
        (simulations - sim_mean) * (evaluation - obs_mean), axis=0, dtype=np.float64
    )
    r_den = np.sqrt(
        np.sum((simulations - sim_mean) ** 2, axis=0, dtype=np.float64)
        * np.sum((evaluation - obs_mean) ** 2, dtype=np.float64)
    )
    r = r_num / r_den
    # calculate error in spread of flow alpha
    alpha = np.std(simulations, axis=0) / np.std(evaluation, dtype=np.float64)
    # calculate error in volume beta (bias of mean discharge)
    beta = np.sum(simulations, axis=0, dtype=np.float64) / np.sum(
        evaluation, dtype=np.float64
    )
    # calculate the Kling-Gupta Efficiency KGE
    kge_ = 1 - np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2)

    return np.vstack((kge_, r, alpha, beta))


def extract_sims(prefix, var_list={"pft": [], "col": [], "const": []}):
    tvec = pd.date_range("2015-01-01", "2020-12-31", freq="1D")
    tvec = tvec[(tvec.month != 2) | (tvec.day != 29)]

    pft_list = [2, 3, 11, 12]
    hol_add = 17

    collection_ts = {}
    collection_const = pd.DataFrame(
        np.nan, index=["hummock", "hollow"], columns=var_list["const"]
    )

    for plot in chamber_list_complete:
        print(plot)
        path_data = os.path.join(
            os.environ["PROJDIR"],
            "E3SM",
            "output",
            f"{prefix}_plot{plot:02d}_US-SPR_ICB20TRCNPRDCTCBC",
            "run",
        )
        flist_pft = sorted(glob(os.path.join(path_data, "*.h2.*.nc")))[:-1]
        flist_col = sorted(glob(os.path.join(path_data, "*.h1.*.nc")))[:-1]
        flist_const = sorted(glob(os.path.join(path_data, "*.h0.*.nc")))[:-1]

        var_list_pft = var_list["pft"]
        hr = xr.open_mfdataset(flist_pft, decode_times=False)
        for var in var_list_pft:
            for pft in pft_list:
                if var == "TSOI_AVG":
                    rootfr = hr["ROOTFR"][:, :, pft].values
                    hr2 = xr.open_mfdataset(flist_col)
                    collection_ts[(plot, var, pft, "hummock")] = np.sum(
                        (hr2["TSOI"][:, :, 0].values - 273.15) * rootfr, axis=1
                    )
                    collection_ts[(plot, var, pft, "hollow")] = np.sum(
                        (hr2["TSOI"][:, :, 1].values - 273.15) * rootfr, axis=1
                    )
                    hr2.close()
                elif var == "SWC_AVG":
                    rootfr = hr["ROOTFR"][:, :, pft].values
                    hr2 = xr.open_dataset(flist_const[0])
                    dzsoi = hr2["DZSOI"].values
                    hr2.close()
                    hr2 = xr.open_mfdataset(flist_col)
                    collection_ts[(plot, var, pft, "hummock")] = np.sum(
                        hr2["SOILLIQ"][:, :, 0].values
                        / dzsoi[:, 0].reshape(1, -1)
                        / 1000
                        * rootfr,
                        axis=1,
                    )
                    collection_ts[(plot, var, pft, "hollow")] = np.sum(
                        hr2["SOILLIQ"][:, :, 1].values
                        / dzsoi[:, 1].reshape(1, -1)
                        / 1000
                        * rootfr,
                        axis=1,
                    )
                    hr2.close()
                elif var == "ZWT":
                    hr2 = xr.open_mfdataset(flist_col)
                    collection_ts[(plot, var, pft, "hummock")] = (
                        0.3 - hr2["ZWT"][:, 0].values
                    )
                    collection_ts[(plot, var, pft, "hollow")] = (
                        hr2["H2OSFC"][:, 1] / 1000.0 - hr2["ZWT"][:, 1]
                    )
                    hr2.close()
                else:
                    collection_ts[(plot, var, pft, "hummock")] = hr[var][:, pft].values
                    collection_ts[(plot, var, pft, "hollow")] = hr[var][
                        :, pft + hol_add
                    ].values
        hr.close()

        var_list_col = var_list["col"]
        hr = xr.open_mfdataset(flist_col, decode_times=False)
        for var in var_list_col:
            if var.startswith("TSOI_"):
                layer = int(var.split("_")[1])
                collection_ts[(plot, var, 0, "hummock")] = (
                    hr["TSOI"][:, layer - 1, 0].values - 273.15
                )
                collection_ts[(plot, var, 0, "hollow")] = (
                    hr["TSOI"][:, layer - 1, 1].values - 273.15
                )
            elif var.startswith("SWC_"):
                layer = int(var.split("_")[1])
                hr2 = xr.open_dataset(flist_const[0])
                dzsoi = hr2["DZSOI"].values
                hr2.close()
                collection_ts[(plot, var, 0, "hummock")] = (
                    hr["SOILLIQ"][:, layer - 1, 0].values / dzsoi[layer - 1, 0] / 1000
                )
                collection_ts[(plot, var, 0, "hollow")] = (
                    hr["SOILLIQ"][:, layer - 1, 1].values / dzsoi[layer - 1, 1] / 1000
                )
            elif var.startswith("H2OSOI_"):
                layer = int(var.split("_")[1])
                collection_ts[(plot, var, 0, "hummock")] = hr["H2OSOI"][
                    :, layer - 1, 0
                ].values
                collection_ts[(plot, var, 0, "hollow")] = hr["H2OSOI"][
                    :, layer - 1, 1
                ].values
            elif var == "SMP_MAX":
                collection_ts[(plot, var, 0, "hummock")] = np.nanmax(
                    hr["SMP"][:, :, 0].values, axis=1
                )  # mm
                collection_ts[(plot, var, 0, "hollow")] = np.nanmax(
                    hr["SMP"][:, :, 1].values, axis=1
                )
            else:
                collection_ts[(plot, var, 0, "hummock")] = hr[var][:, 0].values
                collection_ts[(plot, var, 0, "hollow")] = hr[var][:, 1].values
        hr.close()

        # some time-constant variables are in the h0 files
        var_list_const = var_list["const"]
        hr = xr.open_mfdataset(flist_const, decode_times=False)
        for var in var_list_const:
            if var == "SUCSAT":
                collection_const.loc[:, "SUCSAT"] = pd.Series(
                    hr["SUCSAT"][2, :].values, index=["hummock", "hollow"]
                )
            elif var == "WATSAT":
                collection_const.loc[:, "WATSAT"] = pd.Series(
                    hr["WATSAT"][2, :].values, index=["hummock", "hollow"]
                )
            elif var == "BSW":
                collection_const.loc[:, "BSW"] = pd.Series(
                    hr["BSW"][2, :].values, index=["hummock", "hollow"]
                )
            else:
                raise Exception("Not implemented")
        hr.close()

    collection_ts = pd.DataFrame(collection_ts)
    collection_ts.index = tvec
    collection_ts.columns.names = ["plot", "variable", "pft", "topo"]

    return collection_ts, collection_const


def extract_xys(var_list={"pft": [], "col": [], "const": []}):
    tvec = pd.date_range("2015-01-01", "2020-12-31", freq="1D")
    tvec = tvec[(tvec.month != 2) | (tvec.day != 29)]

    pft_list = [2, 3, 11, 12]
    hol_add = 17

    collection_ts = {}
    collection_const = pd.DataFrame(
        np.nan, index=["hummock", "hollow"], columns=var_list["const"]
    )

    count = 0
    for plot, plot_name in zip(chamber_list_complete, chamber_list_names_complete):
        count = count + 1

        if count >= 2:
            continue

        path_data = os.path.join(
            os.environ["SCRATCH"],
            "../y9s",
            f"20230414_spruceroot_{plot_name}_US-SPR_ICB20TRCNPRDCTCBC",
            "run",
        )
        print(path_data)
        flist_pft = sorted(glob(os.path.join(path_data, "*.h3.*.nc")))[:-1, :]
        flist_col = sorted(glob(os.path.join(path_data, "*.h2.*.nc")))[:-1, :]
        flist_const = sorted(glob(os.path.join(path_data, "*.h0.*.nc")))

        # skip 2021
        flist_pft = [p for p in flist_pft if "2021" not in p]
        flist_col = [p for p in flist_col if "2021" not in p]
        flist_const = [p for p in flist_const if "2021" not in p]

        var_list_pft = var_list["pft"]
        hr = xr.open_mfdataset(flist_pft, decode_times=False)
        for var in var_list_pft:
            for pft in pft_list:
                if var == "TSOI_AVG":
                    rootfr = hr["ROOTFR"][:, :, pft].values
                    hr2 = xr.open_mfdataset(flist_col)
                    collection_ts[(plot, var, pft, "hummock")] = np.sum(
                        (hr2["TSOI"][:, :, 0].values - 273.15) * rootfr, axis=1
                    )
                    collection_ts[(plot, var, pft, "hollow")] = np.sum(
                        (hr2["TSOI"][:, :, 1].values - 273.15) * rootfr, axis=1
                    )
                    hr2.close()
                elif var == "SWC_AVG":
                    rootfr = hr["ROOTFR"][:, :, pft].values
                    hr2 = xr.open_dataset(flist_const[0])
                    dzsoi = hr2["DZSOI"].values
                    hr2.close()
                    hr2 = xr.open_mfdataset(flist_col)
                    collection_ts[(plot, var, pft, "hummock")] = np.sum(
                        hr2["SOILLIQ"][:, :, 0].values
                        / dzsoi[:, 0].reshape(1, -1)
                        / 1000
                        * rootfr,
                        axis=1,
                    )
                    collection_ts[(plot, var, pft, "hollow")] = np.sum(
                        hr2["SOILLIQ"][:, :, 1].values
                        / dzsoi[:, 1].reshape(1, -1)
                        / 1000
                        * rootfr,
                        axis=1,
                    )
                    hr2.close()
                elif var == "ZWT":
                    hr2 = xr.open_mfdataset(flist_col)
                    collection_ts[(plot, var, pft, "hummock")] = (
                        0.3 - hr2["ZWT"][:, 0].values
                    )
                    collection_ts[(plot, var, pft, "hollow")] = (
                        hr2["H2OSFC"][:, 1] / 1000.0 - hr2["ZWT"][:, 1]
                    )
                    hr2.close()
                else:
                    collection_ts[(plot, var, pft, "hummock")] = hr[var][:, pft].values
                    collection_ts[(plot, var, pft, "hollow")] = hr[var][
                        :, pft + hol_add
                    ].values
        hr.close()

        var_list_col = var_list["col"]
        hr = xr.open_mfdataset(flist_col, decode_times=False)
        for var in var_list_col:
            if var.startswith("TSOI_"):
                layer = int(var.split("_")[1])
                collection_ts[(plot, var, 0, "hummock")] = (
                    hr["TSOI"][:, layer - 1, 0].values - 273.15
                )
                collection_ts[(plot, var, 0, "hollow")] = (
                    hr["TSOI"][:, layer - 1, 1].values - 273.15
                )
            elif var.startswith("SWC_"):
                layer = int(var.split("_")[1])
                hr2 = xr.open_dataset(flist_const[0])
                dzsoi = hr2["DZSOI"].values
                hr2.close()
                collection_ts[(plot, var, 0, "hummock")] = (
                    hr["SOILLIQ"][:, layer - 1, 0].values / dzsoi[layer - 1, 0] / 1000
                )
                collection_ts[(plot, var, 0, "hollow")] = (
                    hr["SOILLIQ"][:, layer - 1, 1].values / dzsoi[layer - 1, 1] / 1000
                )
            elif var.startswith("H2OSOI_"):
                layer = int(var.split("_")[1])
                collection_ts[(plot, var, 0, "hummock")] = hr["H2OSOI"][
                    :, layer - 1, 0
                ].values
                collection_ts[(plot, var, 0, "hollow")] = hr["H2OSOI"][
                    :, layer - 1, 1
                ].values
            elif var == "SMP_MAX":
                collection_ts[(plot, var, 0, "hummock")] = np.nanmax(
                    hr["SMP"][:, :, 0].values, axis=1
                )  # mm
                collection_ts[(plot, var, 0, "hollow")] = np.nanmax(
                    hr["SMP"][:, :, 1].values, axis=1
                )
            else:
                collection_ts[(plot, var, 0, "hummock")] = hr[var][:, 0].values
                collection_ts[(plot, var, 0, "hollow")] = hr[var][:, 1].values
        hr.close()

        # some time-constant variables are in the h0 files
        var_list_const = var_list["const"]
        hr = xr.open_mfdataset(flist_const, decode_times=False)
        for var in var_list_const:
            if var == "SUCSAT":
                collection_const.loc[:, "SUCSAT"] = pd.Series(
                    hr["SUCSAT"][2, :].values, index=["hummock", "hollow"]
                )
            elif var == "WATSAT":
                collection_const.loc[:, "WATSAT"] = pd.Series(
                    hr["WATSAT"][2, :].values, index=["hummock", "hollow"]
                )
            elif var == "BSW":
                collection_const.loc[:, "BSW"] = pd.Series(
                    hr["BSW"][2, :].values, index=["hummock", "hollow"]
                )
            else:
                raise Exception("Not implemented")
        hr.close()

    collection_ts = pd.DataFrame(collection_ts)
    collection_ts.index = tvec
    collection_ts.columns.names = ["plot", "variable", "pft", "topo"]

    return collection_ts, collection_const


def read_extract_sims_ts(prefix):
    collection_ts = pd.read_csv(
        os.path.join(path_out, "extract", prefix, "analysis_ts.csv"),
        index_col=0,
        header=[0, 1, 2, 3],
        parse_dates=True,
    )
    indices = collection_ts.columns.to_list()
    indices = [(int(i), j, int(k), l) for i, j, k, l in indices]
    collection_ts.columns = pd.MultiIndex.from_tuples(
        indices, names=collection_ts.columns.names
    )
    return collection_ts


def read_sims_tair_daily():
    prefix = "20221212"  # identical for any
    collection_ts = read_extract_sims_ts(prefix)
    temperature = (
        0.64 * collection_ts.loc[:, (slice(None), "TBOT", 0, "hummock")]
        + 0.36 * collection_ts.loc[:, (slice(None), "TBOT", 0, "hollow")].values
    )
    temperature.columns = temperature.columns.droplevel([1, 2, 3])
    return temperature - 273.15


def read_sims_tair_annual():
    temperature = read_sims_tair_daily()
    temperature = temperature.groupby(temperature.index.year).mean()
    return temperature


def read_obs_tair_annual():
    obs_data = pd.read_excel(
        os.path.join(path_input, "SPRUCE C Budget Summary 28Apr2022EXP.xlsx"),
        sheet_name="DataForPythonRead",
        skiprows=1,
        engine="openpyxl",
    )
    obs_data = obs_data.set_index(["Year", "Plot"]).sort_index(axis=0)
    obs_data = obs_data.sort_index()
    t2m_obs = obs_data.loc[:, "Mean Annual Temp. at 2 m"].unstack()
    t2m_obs.columns = [int(p[1:]) for p in t2m_obs.columns]
    t2m_obs = t2m_obs.loc[:, chamber_list]
    t2m_obs.loc[2015, :] = np.nan
    t2m_obs = t2m_obs.drop(2021, axis=0)
    t2m_obs = t2m_obs.sort_index(axis=0)
    return t2m_obs


def read_leaf_sos():
    pheno_obs = {}
    for var in ["EN", "DN", "SH"]:
        hr = xr.open_dataset(os.path.join(path_intrim, "spruce_validation_data.nc"))
        temp = hr["pheno_dates_lai"].loc[:, "SOS", :, var]
        temp = pd.DataFrame(
            temp.values,
            index=temp["year"],
            columns=["%02d" % i for i in temp["chamber"]],
        )
        pheno_obs[var] = (
            temp.dropna(axis=1, how="all").dropna(axis=0, how="all").sort_index(axis=1)
        )
        hr.close()
    return pheno_obs


def read_leaf_eos():
    pheno_obs = {}
    for var in ["EN", "DN", "SH"]:
        hr = xr.open_dataset(os.path.join(path_intrim, "spruce_validation_data.nc"))
        temp = hr["pheno_dates_lai"].loc[:, "EOS", :, var]
        temp = pd.DataFrame(
            temp.values,
            index=temp["year"],
            columns=["%02d" % i for i in temp["chamber"]],
        )
        pheno_obs[var] = (
            temp.dropna(axis=1, how="all").dropna(axis=0, how="all").sort_index(axis=1)
        )
        hr.close()
    return pheno_obs


def read_obs_tsoi_daily():
    tvec = pd.date_range("2015-01-01", "2021-12-31", freq="1D")
    tvec = tvec[~((tvec.month == 2) & (tvec.day == 29))]

    annt2m = {}
    tsoi_collect = {}  # '2m', '10cm'
    for fid in chamber_levels.keys():
        env = pd.read_csv(
            os.path.join(
                path_input,
                "WEW_Complete_Environ_20220518",
                "WEW PLOT_{}_Complete_Environ_20220518.csv".format(fid),
            )
        )
        env = env.loc[(env["Year"] >= 2015) & (env["Year"] <= 2021), :]
        env.index = pd.DatetimeIndex(env["TIMESTAMP"])
        env = env.replace(to_replace="^\s+", value=np.nan, regex=True).sort_index()
        env = env.loc[~env.index.duplicated(keep="first"), :]

        # These do not match ELM input, use ELM input
        tseries = (env["TA_2_0__1"].astype(float) + env["TA_2_0__2"].astype(float)) / 2
        ##tmax      = tseries.groupby(env.index.year * 1000 + env.index.month * 100 + env.index.day).max()
        ##tmin      = tseries.groupby(env.index.year * 1000 + env.index.month * 100 + env.index.day).min()
        ##tavg      = tseries.groupby(env.index.year * 1000 + env.index.month * 100 + env.index.day).mean()
        annt2m[fid] = tseries.groupby(env.index.year).mean()
        tseries = tseries.resample("1d").mean()
        tseries = tseries.loc[~((tseries.index.month == 2) & (tseries.index.day == 29))]

        # Average hollow & hummock @ 10cm below surface
        temp = (
            env.loc[:, ["TS_ 10__A3", "TS_ 10__B3", "TS_ 10__C3"]].mean(axis=1) * 0.64
            + env.loc[:, ["TS_Hummock_A2", "TS_Hummock_B2"]].mean(axis=1).values * 0.36
        )
        temp = temp.astype(float)
        # has NaNs if keep hourly level
        temp = temp.resample("1d").mean()
        temp = temp.loc[~((temp.index.month == 2) & (temp.index.day == 29))]

        tsoi = pd.DataFrame(np.nan, index=tvec, columns=["2m", "10cm"])
        tsoi.loc[temp.index, "10cm"] = temp.values

        # check: tseries.loc[(tsoi.index[0] - timedelta(days = 119, hours = 23, minutes = 30)):(tsoi.index[0])].mean()
        tair = tseries.rolling("21d").mean()
        tsoi.loc[tair.index, "2m"] = tair.values

        # fill NaNs
        for col in tsoi.columns:
            tsoi.loc[:, col] = tsoi.loc[:, col].interpolate(
                limit=5, limit_direction="both"
            )

        # skip NaNs in the beginning
        tsoi = tsoi.loc[tsoi.index >= datetime(2015, 11, 1), :]

        narows = np.where(tsoi.isna().any(axis=1))[0]
        if len(narows) > 0:
            print(fid, tsoi["10cm"][narows])
            import pdb

            pdb.set_trace()
            raise Exception("check")

        tsoi = dict(tsoi.dropna(axis=0, how="any"))

        for k in tsoi.keys():
            if k not in tsoi_collect:
                tsoi_collect[k] = pd.DataFrame({fid: tsoi[k]})
            else:
                tsoi_collect[k][fid] = tsoi[k]
    annt2m = pd.DataFrame(annt2m)
    annt2m = annt2m.sort_index(axis=1)
    for k in tsoi_collect.keys():
        tsoi_collect[k] = tsoi_collect[k].sort_index(axis=1)

    return tsoi_collect, annt2m
