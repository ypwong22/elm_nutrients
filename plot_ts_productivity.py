""" Sensitivity of annual mean quantities to temperature, scatter plots """
import matplotlib.pyplot as plt
import os
import pandas as pd
from scipy.stats import linregress
import numpy as np
from utils.analysis import get_spruce_carbonfluxes


def fit_line(x, y):
    filt = ~np.isnan(x) & ~np.isnan(y)
    x = x[filt]
    y = y[filt]
    res = linregress(x, y)
    xnew = np.linspace(x.min(), x.max(), 3)
    ynew = res.slope * xnew + res.intercept
    r2 = res.rvalue**2  # coefficient of determination
    return xnew, ynew, res.slope, res.intercept, r2


prefix = "20231113_2"
#prefix = "UQ_20240315"
outdir = os.path.join(os.environ['PROJDIR'], 'ELM_Phenology', 'output', 'extract', prefix)

# T: 0, 2.25, 4.5, 6.75, 9; first ambient, then elevated CO2
chamber_list = {
    "amb": ["P06", "P20", "P13", "P08", "P17"],
    "elev": ["P19", "P11", "P04", "P16", "P10"],
}
sim_data = pd.read_csv(os.path.join(outdir, 'extract_ts_productivity.csv'), 
                       index_col=[0, 1, 2])
sim_data = sim_data.loc['average', :]
sim_data = sim_data.loc[sim_data.index.get_level_values(0) != 2020, :]
sim_varname = ["TOTVEGC_ABG_pima", "TOTVEGC_ABG_lala", "TOTVEGC_ABG_shrub", 
               "AGNPP_pima", "AGNPP_lala", "AGNPP_shrub", "AGNPP_tree", "AGNPP_tree_shrub", 
               "BGNPP_tree_shrub", "BG_to_AG", # ratio of AGNPP to BGNPP of tree+shrub
               "NPP_moss", "HR", "NEE"]
sim_var_unobs = ["MR_pima", 'MR_lala', 'MR_shrub', "SMINN_30", "SOLUTIONP_30",
                 "FPG_pima", "FPG_lala", "FPG_shrub", "FPG_P_pima", "FPG_P_lala",
                 "FPG_P_shrub", 'BG_to_AG_pima', 'BG_to_AG_lala', 'BG_to_AG_shrub',
                 "FPI","FPI_P", "TOTSOMC"]
sim_varname += sim_var_unobs

obs_data = get_spruce_carbonfluxes()
obs_varname = ["ABGbiomass evergreen conifer", "ABGbiomass deciduous conifer",
               "ABGbiomass shrub", "ABGnpp evergreen conifer",
               "ABGnpp deciduous conifer", "ANPP Shrub (~50%C)",
               "ANPP Tree (~48%C)", "AGNPP_tree_shrub", 
               'BNPP Tree & Shrub', "BG_to_AG", "NPP Sphag.", "RHCO2", "NCE"]
obs_varname += sim_var_unobs
title_list = ["Spruce biomass", "Tamarack biomass", "Shrub biomass",
              "Spruce ANPP", "Tamarack ANPP", "Shrub ANPP (Hanson | Salmon)", 
              "Tree ANPP (Hanson | Salmon)", "Tree+Shrub ANPP (Hanson | Salmon)",
              "Tree+Shrub BNPP", "BNPP:ANPP tree+shrub", "Sphagnum NPP", "HR", "NEE"]
title_list += sim_var_unobs


for varname, ov, title in zip(sim_varname, obs_varname, title_list):
    fig, ax = plt.subplots(figsize=(6, 4))

    if varname in ['AGNPP_tree', 'AGNPP_tree_shrub', 'AGNPP_shrub']:
        h = [None] * 6
        clist = ["r", "b", "g", "r", "b", "g"]
        mlist = ["r", "b", "g", "none", "none", "none"]
        llist = ["-", "-", "-", "--", "--", "--"]
    else:
        h = [None] * 4
        clist = ["r", "b", "r", "b"]
        mlist = ["r", "b", "none", "none"]
        llist = ["-", "-", "--", "--"]

    count = 0

    for co2, CO2 in zip(["amb", "elev"], ["ACO2", "ECO2"]):
        if varname == "BG_to_AG":
            sim_temp = sim_data.loc[chamber_list[co2], "BGNPP_tree_shrub"] / \
                (sim_data.loc[chamber_list[co2], "AGNPP_tree_shrub"])
        elif "BG_to_AG_" in varname:
            pft = varname.split('_')[-1]
            sim_temp = sim_data.loc[chamber_list[co2], f'BGNPP_{pft}'] / \
                       sim_data.loc[chamber_list[co2], f'AGNPP_{pft}']
        #elif 'delta' in varname:
        #    varname2 = varname.replace('delta_', '')
        #    # subtract the chamber 7 baseline
        #    sim_temp = sim_data.loc[chamber_list[co2], varname2] - \
        #        sim_data.loc[("P07", 2015), varname2]
        else:
            sim_temp = sim_data.loc[chamber_list[co2], varname]
            if varname in ["NEE", "HR"]:
                sim_temp = -1 * sim_temp

        sim_T = sim_data.loc[chamber_list[co2], "TBOT"]

        (h[count],) = ax.plot(
            sim_T, sim_temp, "o", color=clist[count], markerfacecolor=mlist[count]
        )
        xnew, ynew, slope, intercept, r2 = fit_line(sim_T, sim_temp)
        ax.plot(xnew, ynew, ls=llist[count], color=clist[count], markerfacecolor=mlist[count])
        ax.text(0.4, 0.65 + count * 0.08, "y$_{MOD\_" + CO2 + "}$="
            + f"{slope:.2f}x+{intercept:.2f}  "+ "R$^2$="+ f"{r2:.3f}", 
            color=clist[count], transform=ax.transAxes)
        count = count + 1

        if varname in sim_var_unobs:
            continue

        filt = obs_data["Plot"].isin(chamber_list[co2])

        if varname == "BG_to_AG":
            obs_temp = obs_data.loc[filt, "BNPP Tree & Shrub"] / \
                (obs_data.loc[filt, 'ANPP Tree (~48%C)'] + \
                 obs_data.loc[filt, 'ANPP Shrub (~50%C)'])
        elif varname == "AGNPP_tree_shrub":
            obs_temp = (obs_data.loc[filt, 'ANPP Tree (~48%C)'] + \
                 obs_data.loc[filt, 'ANPP Shrub (~50%C)'])
        else:
            obs_temp = obs_data.loc[filt, ov]
        obs_T = obs_data.loc[filt, "Mean Annual Temp. at 2 m"]

        (h[count],) = ax.plot(
            obs_T, obs_temp, "o", color=clist[count], markerfacecolor=mlist[count]
        )
        xnew, ynew, slope, intercept, r2 = fit_line(obs_T, obs_temp)
        ax.plot(xnew, ynew, ls=llist[count], color=clist[count], markerfacecolor=mlist[count])
        ax.text(0.4, 0.65 + count * 0.08,
            "y$_{OBS\_"+ CO2 + "}$=" + f"{slope:.2f}x+{intercept:.2f}  " + "R$^2$=" + f"{r2:.3f}",
            color=clist[count], transform=ax.transAxes)
        count = count + 1

        # add Verity's data apart from Paul's in case of dual observations
        if varname in ['AGNPP_tree', 'AGNPP_tree_shrub', 'AGNPP_shrub']:
            if varname == 'AGNPP_tree':
                obs_temp = (obs_data.loc[filt, 'ABGnpp evergreen conifer'] + \
                        obs_data.loc[filt, 'ABGnpp deciduous conifer'])
            elif varname == 'AGNPP_tree_shrub':
                obs_temp = (obs_data.loc[filt, 'ABGnpp evergreen conifer'] + \
                        obs_data.loc[filt, 'ABGnpp deciduous conifer'] + \
                        obs_data.loc[filt, 'ABGnpp shrub'])
            elif varname == 'AGNPP_shrub':
                obs_temp = obs_data.loc[filt, 'ABGnpp shrub']
            (h[count],) = ax.plot(
                obs_T, obs_temp, "o", color='green', markerfacecolor=mlist[count]
            )
            xnew, ynew, slope, intercept, r2 = fit_line(obs_T, obs_temp)
            ax.plot(xnew, ynew, ls=llist[count], color=clist[count], markerfacecolor=mlist[count])
            ax.text(0.4, 0.65 + count * 0.08,
                "y$_{OBS\_"+ CO2 + "}$=" + f"{slope:.2f}x+{intercept:.2f}  " + "R$^2$=" + f"{r2:.3f}",
                color=clist[count], transform=ax.transAxes)
            count = count + 1

    ax.axhline(0.0, ls=":", color="k", lw=0.5)

    if varname in sim_var_unobs:
        ax.legend(
            h, ["MOD_ACO2", "MOD_ECO2"], ncol=4, loc=[-0.1, -0.25]
        )
    else:
        if len(h) == 4:
            ax.legend(
                h, ["MOD_ACO2", "OBS_ACO2", "MOD_ECO2", "OBS_ECO2"], ncol=4, loc=[-0.1, -0.25]
            )
        else:
            ax.legend(h, ["MOD_ACO2", "OBS_ACO2", "OBS_2_ACO2", 
                          "MOD_ECO2", "OBS_ECO2", "OBS_2_ECO2"], 
                      ncol=4, loc=[-0.1, -0.25])
    ax.set_title(title)

    if "biomass" in title or "mineral" in title:
        ylabel = "g m-2"
    elif "N limit" in title or "P limit" in title or 'BG_to_AG' in title:
        ylabel = ""
    else:
        ylabel = title + " gC m-2 yr-1"
    ax.set_ylabel(ylabel)
    ax.set_xlabel("Mean Annual Temperature ($^o$C)")
    fig.savefig(
        os.path.join(outdir, f"result_{prefix}_plot_{varname}.png"), dpi=600.0, bbox_inches="tight"
    )
    plt.close(fig)
