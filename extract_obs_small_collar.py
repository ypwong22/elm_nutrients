"""
Flux partitioning: Q10-based Reco estimation and GPP model fitting.

Workflow:
  1. Separate nighttime (PAR < 50) and daytime (PAR >= 50) records.
  2. Fit Q10 model to nighttime NEE (= Reco at night) vs Tsoil.
  3. Predict daytime Reco from Q10 model + daytime Tsoil.
  4. Derive daytime GPP = Reco_predicted - NEE  (since NEE = GPP - Reco → GPP = Reco - NEE).
     NOTE: sign convention — NEE is negative when ecosystem is a net sink.
           If your NEE is defined as GPP - Reco (positive = uptake), then
           GPP = NEE + Reco.  Adjust SIGN_CONVENTION below.
  5. Fit four GPP models (PAR; PAR+Wh; PAR+Tair; PAR+Wh+Tair).
  6. Plot observed vs fitted Reco and GPP with RMSE and AIC.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os

from utils.nee_flux_partitioning import (
    fit_q10,
    fit_gpp,
)


# Column names for training (adjust to match your DataFrame)
COL_DATETIME = "Timestamp"
COL_NEE      = "CO2_flux"
COL_PAR      = "PAR"
COL_TAIR     = "TA_0_5"
COL_TSOIL    = "TS_5_B2" # "collar_soil_temp"
COL_WH       = "WTD (m)"

# Column names for prediction (adjust to match your DataFrame)
PRED_PAR     = "PAR_2"
PRED_TAIR    = "TA_0_5"
PRED_TSOIL   = "TS_5_B2"
PRED_WH      = "WTD (m)"


# "ecology":  NEE = Reco - GPP  (negative NEE = net uptake)
# "atmosphere": NEE = GPP - Reco
SIGN_CONVENTION = "ecology"

PAR_THRESHOLD = 50.0
Q10_T_REF = 15.0  # reference temperature for Q10 fit

# =====================================================================
# HELPER FUNCTIONS
# =====================================================================

def calc_rmse(obs, pred):
    return np.sqrt(np.mean((obs - pred) ** 2))


def calc_aic(obs, pred, k):
    n = len(obs)
    ss_res = np.sum((obs - pred) ** 2)
    return n * np.log(ss_res / n) + 2 * k


# =====================================================================
# Collect the small collar measurements into well-formatted 
# single dataframe
# =====================================================================
# It turns out 2023 data is much more complete than 2022; 2022 only has a few months
# Focus on 2023 then
# Aggregate to half-hourly to be compatible with water table observation frequency
nee_data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ELM_Nutrients', 'input', 'Data and Guide for SPRUCE.104', 
                                    'SPRUCE_GHG_Flux_15_min_2023.csv'), 
                       index_col = 0, na_values = -9999, parse_dates=True)
nee_data = nee_data.groupby(['date_time_30_min_CST', 'plot', 'warming_treatment', 'CO2_treatment', 'collar_number']).mean()
nee_data = nee_data[['CO2_flux','PAR','collar_soil_temp','collar_VWC']].reset_index()
nee_data.rename({'date_time_30_min_CST': 'Timestamp'}, axis = 1, inplace = True)
nee_data['Timestamp'] = pd.to_datetime(nee_data['Timestamp'])

#
warming_list = nee_data['warming_treatment'].drop_duplicates().sort_values().values
co2_list = nee_data['CO2_treatment'].drop_duplicates().values
plot_list = nee_data['plot'].drop_duplicates().values

# Environmental soil temperature, air temperature, and PAR data
env_data = []
for plot in plot_list:
    warming = nee_data.loc[nee_data['plot'] == plot, 'warming_treatment'].values[0]
    co2 = nee_data.loc[nee_data['plot'] == plot, 'CO2_treatment'].values[0]
    # Obtain the PAR and 5-10cm soil temperature from complete environmental data
    temp = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ELM_Nutrients', 'input', 'WEW_Complete_Environ_20250311',
                                        f'WEW PLOT_{plot[1:]}_Complete_Environ_20250311.csv'),
                           na_values = -9999)
    temp = temp[['TIMESTAMP','TA_0_5','TS_0_B1','TS_5_B2','TS_10_B3','PAR_2']].dropna(axis = 0, how = 'any')
    temp['TIMESTAMP'] = pd.to_datetime(temp['TIMESTAMP'])
    temp = temp.loc[temp['TIMESTAMP'].dt.year == 2023, :]
    temp.rename({'TIMESTAMP': 'Timestamp'}, axis = 1, inplace=True)
    temp['plot'] = plot
    temp['warming_treatment'] = warming
    temp['CO2_treatment'] = co2
    env_data.append(temp)
env_data = pd.concat(env_data, axis = 0)

# Water table depth data: relative to hollow height, in meters
temp = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ELM_Nutrients', 'input', 'Data and Guide for SPRUCE.079',
                                'WEW_WT_HalfHour_Updated_20250408.csv'), na_values = -9999)
wtd_data = []
for plot in plot_list:
    warming = nee_data.loc[nee_data['plot'] == plot, 'warming_treatment'].values[0]
    co2 = nee_data.loc[nee_data['plot'] == plot, 'CO2_treatment'].values[0]

    if f'WTNorm_{plot[1:]}' in temp.columns:
        subset = temp[['TIMESTAMP', f'WTNorm_{plot[1:]}']].copy()
    else:
        subset = temp[['TIMESTAMP', f'WTNorm{plot[1:]}']].copy()
    subset.columns = ['Timestamp', 'WTD (m)']
    subset['Timestamp'] = pd.to_datetime(subset['Timestamp'])
    subset = subset.loc[subset['Timestamp'].dt.year == 2023, :]

    subset['plot'] = plot
    subset['warming_treatment'] = warming
    subset['CO2_treatment'] = co2
    wtd_data.append(subset)
wtd_data = pd.concat(wtd_data, axis = 0)

# 
env_data.set_index(['plot', 'warming_treatment', 'CO2_treatment', 'Timestamp'], inplace = True)
wtd_data.set_index(['plot', 'warming_treatment', 'CO2_treatment', 'Timestamp'], inplace = True)
env_wtd = pd.concat([env_data, wtd_data], axis=1, join='outer')

#
nee_data.set_index(['plot', 'warming_treatment', 'CO2_treatment', 'Timestamp', 'collar_number'], inplace = True)
# Join to nee_data on the 4 common levels — env/wtd values broadcast to all collar_numbers
df_all_plots = nee_data.join(env_wtd, on=['plot', 'warming_treatment', 'CO2_treatment', 'Timestamp'], 
                             how='outer')


# -------------------------------------------------------------------
# Loop through all the plots individually
# -------------------------------------------------------------------
nighttime_params = {
    "Q10 only": pd.DataFrame(np.nan, index = plot_list, 
                             columns = ["Base","Q10","Base_std","Q10_std","RMSE","AIC"]),
    "Q10+Wh": pd.DataFrame(np.nan, index = plot_list, 
                           columns = ["Base","rw","Q10","Base_std","rw_std","Q10_std","RMSE","AIC"])
}

daytime_params = {
    "PAR only": pd.DataFrame(np.nan, index = plot_list,
                             columns = ["alpha","gmax","alpha_std","gmax_std","RMSE","AIC"]),
    "PAR+Wh": pd.DataFrame(np.nan, index = plot_list, 
                           columns = ["alpha","gmax","gw","alpha_std","gmax_std","gw_std","RMSE","AIC"]),
    "PAR+Tair": pd.DataFrame(np.nan, index = plot_list,
                             columns = ["alpha","gmax","Tmin","Tmax","Topt","alpha_std","gmax_std","Tmin_std","Tmax_std","Topt_std","RMSE","AIC"]),
    "PAR+Wh+Tair": pd.DataFrame(np.nan, index = plot_list, 
                                columns = ["alpha","gmax","gw","Tmin","Tmax","Topt","alpha_std","gmax_std","gw_std","Tmin_std","Tmax_std","Topt_std","RMSE","AIC"])
}
save_fit = []
save_pred = []


plot_treatment_co2 = df_all_plots.index.to_frame().iloc[:,:3].drop_duplicates()
for _, (plot, warming, co2) in plot_treatment_co2.iterrows():

    df = df_all_plots.loc[plot,:]

    cols_needed = [COL_NEE, COL_PAR, COL_TAIR, COL_TSOIL, COL_WH]
    fit_df = df.dropna(subset=cols_needed)
    if SIGN_CONVENTION != "ecology":
        fit_df[COL_NEE] *= -1

    cols_preds = [PRED_PAR, PRED_TAIR, PRED_TSOIL, PRED_WH]
    pred_df = env_wtd.loc[plot,:]

    # use the ambient 6.75's PAR to gapfill the elevated 6.75's PAR, because
    # measurements are nonexistant
    if plot == "P16":
        index_1 = env_wtd.loc["P16","PAR_2"].index.get_level_values('Timestamp')
        index_2 = env_wtd.loc["P08","PAR_2"].index.get_level_values('Timestamp')
        overlap = index_1.intersection(index_2)
        pred_df.loc[(slice(None),slice(None),overlap),"PAR_2"] = env_wtd.loc["P08","PAR_2"].loc[(slice(None),slice(None),overlap)].values

    pred_df = pred_df.dropna(subset=cols_preds[1:])

    # =====================================================================
    # STEP 1: FIT Q10 ON NIGHTTIME DATA
    # =====================================================================

    # Nighttime: GPP ≈ 0, so Reco ≈ NEE (ecology) or -NEE (atmosphere)
    # Keep only positive Reco for Q10 fitting 
    fit_night_mask = (fit_df[COL_PAR] < PAR_THRESHOLD) & \
                     (fit_df[COL_NEE].values > 0)
    fit_df_night = fit_df.loc[fit_night_mask].copy()

    print(f"{plot} night records: {len(fit_df_night)}")

    # --- Q10 only ---
    popt, pcov, names, model = fit_q10(fit_df_night[COL_TSOIL], fit_df_night[COL_NEE], Q10_T_REF)
    fitted = model(fit_df_night[COL_TSOIL], *popt)
    k = len(popt)
    rmse = calc_rmse(fit_df_night[COL_NEE], fitted)
    aic = calc_aic(fit_df_night[COL_NEE], fitted, k)
    for nm,v,e in zip(names, popt, np.sqrt(np.diag(pcov))):
        nighttime_params["Q10 only"].loc[plot, nm] = v
        nighttime_params["Q10 only"].loc[plot, f"{nm}_std"] = e
    nighttime_params["Q10 only"].loc[plot, "RMSE"] = rmse
    nighttime_params["Q10 only"].loc[plot, "AIC"] = aic

    # make gapfilled data
    pred_df["Reco_pred Q10 only"] = model(pred_df[PRED_TSOIL].values, *popt)
    fit_df["Reco_pred Q10 only"] = model(fit_df[COL_TSOIL].values, *popt)

    # --- Q10 + water table height ---
    popt, pcov, names, model = fit_q10(fit_df_night[COL_TSOIL], fit_df_night[COL_NEE], Q10_T_REF, fit_df_night[COL_WH])
    fitted = model(np.vstack([fit_df_night[COL_TSOIL], fit_df_night[COL_WH]]), *popt)
    k = len(popt)
    rmse = calc_rmse(fit_df_night[COL_NEE], fitted)
    aic = calc_aic(fit_df_night[COL_NEE], fitted, k)
    for nm,v,e in zip(names, popt, np.sqrt(np.diag(pcov))):
        nighttime_params["Q10+Wh"].loc[plot, nm] = v
        nighttime_params["Q10+Wh"].loc[plot, f"{nm}_std"] = e
    nighttime_params["Q10+Wh"].loc[plot, "RMSE"] = rmse
    nighttime_params["Q10+Wh"].loc[plot, "AIC"] = aic

    # make gapfilled data
    pred_df["Reco_pred Q10+Wh"] = model(np.vstack([pred_df[PRED_TSOIL].values, pred_df[PRED_WH].values]), *popt)
    fit_df["Reco_pred Q10+Wh"] = model(np.vstack([fit_df[COL_TSOIL].values, fit_df[COL_WH].values]), *popt)

    # =====================================================================
    # STEP 2: PREDICT RECO EVERYWHERE, DERIVE DAYTIME GPP
    # =====================================================================
    fit_df["GPP_obs"] = fit_df["Reco_pred Q10+Wh"] - df[COL_NEE]

    fit_day_mask = (fit_df[COL_PAR] > PAR_THRESHOLD) & (fit_df["GPP_obs"] > 0)
    fit_df_day = fit_df.loc[fit_day_mask].copy()
    print(f"{plot} day records:   {len(fit_df_day)}")

    PAR_d = fit_df_day[COL_PAR].values
    GPP_d = fit_df_day["GPP_obs"].values
    Wh_d = fit_df_day[COL_WH].values
    Tair_d = fit_df_day[COL_TAIR].values
    time_day_valid = fit_df_day.index

    pred_day_mask = pred_df[PRED_PAR] > PAR_THRESHOLD

    # =====================================================================
    # STEP 3: FIT FOUR GPP MODELS
    # =====================================================================
    # --- PAR only ---
    popt, pcov, names, model = fit_gpp(PAR_d, GPP_d)
    fitted = model(PAR_d, *popt)
    k = len(popt)
    rmse = calc_rmse(GPP_d, fitted)
    aic = calc_aic(GPP_d, fitted, k)
    for nm,v,e in zip(names, popt, np.sqrt(np.diag(pcov))):
        daytime_params["PAR only"].loc[plot, nm] = v
        daytime_params["PAR only"].loc[plot, f"{nm}_std"] = e
    daytime_params["PAR only"].loc[plot, "RMSE"] = rmse
    daytime_params["PAR only"].loc[plot, "AIC"] = aic
    
    # ------ make gapfilled data ---
    pred_df.loc[:, "GPP_pred PAR only"] = np.where(pred_day_mask, model(pred_df[PRED_PAR].values, *popt), 0)

    # --- PAR + Wh ---
    popt, pcov, names, model = fit_gpp(PAR_d, GPP_d, Wh=Wh_d)
    fitted = model(np.vstack([PAR_d, Wh_d]), *popt)
    k = len(popt)
    rmse = calc_rmse(GPP_d, fitted)
    aic = calc_aic(GPP_d, fitted, k)
    for nm,v,e in zip(names, popt, np.sqrt(np.diag(pcov))):
        daytime_params["PAR+Wh"].loc[plot, nm] = v
        daytime_params["PAR+Wh"].loc[plot, f"{nm}_std"] = e
    daytime_params["PAR+Wh"].loc[plot, "RMSE"] = rmse
    daytime_params["PAR+Wh"].loc[plot, "AIC"] = aic

    # ------ make gapfilled data ---
    pred_df.loc[:, "GPP_pred PAR+Wh"] = np.where(pred_day_mask, model(np.vstack([pred_df[PRED_PAR].values, pred_df[PRED_WH].values]), *popt), 0)

   # --- PAR + Tair ---
    popt, pcov, names, model = fit_gpp(PAR_d, GPP_d, Tair=Tair_d)
    fitted = model(np.vstack([PAR_d, Tair_d]), *popt)
    k = len(popt)
    rmse = calc_rmse(GPP_d, fitted)
    aic = calc_aic(GPP_d, fitted, k)
    for nm,v,e in zip(names, popt, np.sqrt(np.diag(pcov))):
        daytime_params["PAR+Tair"].loc[plot, nm] = v
        daytime_params["PAR+Tair"].loc[plot, f"{nm}_std"] = e
    daytime_params["PAR+Tair"].loc[plot, "RMSE"] = rmse
    daytime_params["PAR+Tair"].loc[plot, "AIC"] = aic

    # ------ make gapfilled data ---
    pred_df.loc[:, "GPP_pred PAR+Tair"] = np.where(pred_day_mask, model(np.vstack([pred_df[PRED_PAR].values, pred_df[PRED_TAIR]]), *popt), 0)

    # --- PAR + Wh + Tair ---
    popt, pcov, names, model = fit_gpp(PAR_d, GPP_d, Wh=Wh_d, Tair=Tair_d)
    fitted = model(np.vstack([PAR_d, Wh_d, Tair_d]), *popt)
    k = len(popt)
    rmse = calc_rmse(GPP_d, fitted)
    aic = calc_aic(GPP_d, fitted, k)
    for nm,v,e in zip(names, popt, np.sqrt(np.diag(pcov))):
        daytime_params["PAR+Wh+Tair"].loc[plot, nm] = v
        daytime_params["PAR+Wh+Tair"].loc[plot, f"{nm}_std"] = e
    daytime_params["PAR+Wh+Tair"].loc[plot, "RMSE"] = rmse
    daytime_params["PAR+Wh+Tair"].loc[plot, "AIC"] = aic

    # ------ make gapfilled data ---
    pred_df.loc[:, "GPP_pred PAR+Wh+Tair"] = np.where(pred_day_mask, model(np.vstack([pred_df[PRED_PAR].values, pred_df[PRED_WH].values, pred_df[PRED_TAIR].values]), *popt), 0)

    # =====================================================================
    # STEP 4: PLOTS WITH PREDICTION
    # =====================================================================
    # ── Compute daily averages for plotting ──
    # fit_df_night / fit_df_day: MultiIndex with chamber + Timestamp → average over chambers then over hours
    ts_night = fit_df_night.index.get_level_values('Timestamp')
    daily_night = fit_df_night.copy()
    daily_night.index = ts_night
    daily_night = daily_night.groupby(pd.Grouper(freq='D')).mean()

    ts_day = fit_df_day.index.get_level_values('Timestamp')
    daily_day = fit_df_day.copy()
    daily_day.index = ts_day
    daily_day = daily_day.groupby(pd.Grouper(freq='D')).mean()

    # pred_df: split by PAR threshold before daily averaging
    daily_pred_index = pd.DatetimeIndex(pred_df.index.get_level_values('Timestamp'))
    pred_df_flat = pred_df.copy()
    pred_df_flat.index = daily_pred_index

    night_mask = pred_df_flat[PRED_PAR] < PAR_THRESHOLD
    day_mask   = pred_df_flat[PRED_PAR] > PAR_THRESHOLD

    reco_cols = [c for c in pred_df_flat.columns if c.startswith("Reco_pred")]
    gpp_cols  = [c for c in pred_df_flat.columns if c.startswith("GPP_pred")]
    other_cols = [c for c in pred_df_flat.columns if c not in reco_cols + gpp_cols]

    daily_pred_reco  = pred_df_flat.loc[night_mask, reco_cols].groupby(pd.Grouper(freq='D')).mean()
    daily_pred_gpp   = pred_df_flat.loc[day_mask,   gpp_cols ].groupby(pd.Grouper(freq='D')).mean()
    daily_pred_other = pred_df_flat[other_cols].groupby(pd.Grouper(freq='D')).mean()

    daily_pred = pd.concat([daily_pred_other, daily_pred_reco, daily_pred_gpp], axis=1)

    # ── Plot 1: Nighttime Reco ──
    fig1, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    colors = ['#1f77b4', '#ff7f0e']
    for i, (key, color) in enumerate(zip(["Q10 only", "Q10+Wh"], colors)):
        ax1 = axes[i]
        ax1.plot(daily_night.index, daily_night[COL_NEE], '.', color='grey',
                 markersize=5, alpha=0.5, label='Observed Reco (night, daily)')
        ax1.plot(daily_pred.index, daily_pred[f"Reco_pred {key}"], '-', color=colors[i],
                 alpha=0.5, label=f'{key} fitted (daily)')
        ax1.set_ylabel("Reco (µmol m⁻² s⁻¹)")
        ax1.legend(loc='upper right', markerscale=4)
        rmse = nighttime_params[key].loc[plot, "RMSE"]
        aic = nighttime_params[key].loc[plot, "AIC"]
        ax1.text(0.02, 0.95,
                 f"RMSE = {rmse:.3f}\nAIC = {aic:.1f}",
                 transform=ax1.transAxes, va='top', fontsize=9,
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    fig1.autofmt_xdate()
    fig1.tight_layout()
    fig1.savefig(os.path.join(os.environ['PROJDIR'], 'ELM_Nutrients', 'output', 'extract',
                              'extract_obs_small_collar', f"reco_fit_{warming}_{co2}.png"), dpi=150)
    print("\nSaved: reco_fit.png")

    # ── Plot 2: Daytime GPP — four models ──
    fig2, axes = plt.subplots(4, 1, figsize=(14, 16), sharex=True)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

    for i, (key, color) in enumerate(zip(["PAR only", "PAR+Wh", "PAR+Tair", "PAR+Wh+Tair"], colors)):
        ax = axes[i]
        ax.plot(daily_day.index, daily_day["GPP_obs"], '.', color='grey', markersize=5,
                alpha=0.4, label='Observed GPP (daily)')
        ax.plot(daily_pred.index, daily_pred[f"GPP_pred {key}"], '-', color=color, alpha=0.4, label=f'{key} fit (daily)')
        ax.set_ylabel("GPP (µmol m⁻² s⁻¹)")
        ax.legend(loc='upper right', markerscale=4)
        rmse = daytime_params[key].loc[plot, "RMSE"]
        aic = daytime_params[key].loc[plot, "AIC"]
        ax.text(0.02, 0.95, f"RMSE = {rmse:.3f}\nAIC = {aic:.1f}",
                transform=ax.transAxes, va='top', fontsize=9,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    fig2.autofmt_xdate()
    fig2.tight_layout()
    fig2.savefig(os.path.join(os.environ['PROJDIR'], 'ELM_Nutrients', 'output', 'extract',
                              'extract_obs_small_collar', f"gpp_fit_{warming}_{co2}.png"), dpi=150)
    print("Saved: gpp_fit_comparison.png")

    # ── Fitted and table ──
    pred_df['Plot'] = plot
    fit_df['Plot'] = plot # add plot labels
    save_pred.append(pred_df)
    save_fit.append(fit_df)

save_pred = pd.concat(save_pred, axis = 0)
save_fit = pd.concat(save_fit, axis = 0)

# Correct GPP values by setting non-growing seasons to 0
md = save_pred.index.get_level_values('Timestamp').month * 100 + save_pred.index.get_level_values('Timestamp').day
save_pred.loc[(md < 515) | (md > 1015), 'GPP_pred PAR only'] = 0.
save_pred.loc[(md < 515) | (md > 1015), 'GPP_pred PAR+Wh'] = 0.
save_pred.loc[(md < 515) | (md > 1015), 'GPP_pred PAR+Tair'] = 0.
save_pred.loc[(md < 515) | (md > 1015), 'GPP_pred PAR+Wh+Tair'] = 0.


# Save all results and parameters
save_fit.to_csv(os.path.join(os.environ['PROJDIR'], 'ELM_Nutrients', 'output', 'extract',
                             'extract_obs_small_collar', f"nee_fit.csv"))
save_pred.to_csv(os.path.join(os.environ['PROJDIR'], 'ELM_Nutrients', 'output', 'extract',
                              'extract_obs_small_collar', f"nee_pred.csv"))
for key in nighttime_params.keys():
    nighttime_params[key].to_csv(os.path.join(os.environ['PROJDIR'], 'ELM_Nutrients', 'output', 'extract',
                                              'extract_obs_small_collar', f"nighttime_params_{key}.csv"))
for key in daytime_params.keys():
    daytime_params[key].to_csv(os.path.join(os.environ['PROJDIR'], 'ELM_Nutrients', 'output', 'extract',
                                            'extract_obs_small_collar', f"daytime_params_{key}.csv"))