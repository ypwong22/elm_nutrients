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


#prefix = "20231112"
#prefix = "20240315"
#prefix = "20240317"
prefix = "20240317_5"
#prefix = "UQ_20240315"
outdir = os.path.join(os.environ['PROJDIR'], 'ELM_Phenology', 'output', 
                      'extract', prefix)

# T: 0, 2.25, 4.5, 6.75, 9; first ambient, then elevated CO2
chamber_list = {
    "amb": ["P06", "P20", "P13", "P08", "P17"],
    "elev": ["P19", "P11", "P04", "P16", "P10"],
}
sim_data = pd.read_csv(os.path.join(outdir, 'extract_ts_productivity.csv'), 
                       index_col=[0, 1, 2])
sim_data = sim_data.loc['average', :]
sim_data = sim_data.loc[sim_data.index.get_level_values(0) != 2020, :]
sim_varname = ["NEE", "delta_TOTVEGC_ABG_pima", "delta_TOTVEGC_ABG_lala", "delta_TOTVEGC_ABG_shrub",
               "delta_AGNPP_pima", "delta_AGNPP_lala", "delta_AGNPP_shrub", 
               "TOTVEGC_ABG_pima", "TOTVEGC_ABG_lala", "TOTVEGC_ABG_shrub", 
               "AGNPP_pima", "AGNPP_lala", "AGNPP_shrub", "AGNPP_tree", "AGNPP_shrub_2", 
               "BGNPP_tree_shrub", "BG_to_AG", # ratio of AGNPP to BGNPP of tree+shrub
               "NPP_moss", "HR"]
sim_var_unobs = ["MR_pima", 'MR_lala', 'MR_shrub',"SMINN","SMINP",
                 "FPG_pima", "FPG_lala", "FPG_shrub", "FPG_P_pima", "FPG_P_lala",
                 "FPG_P_shrub", # "BGNPP_pima", "BGNPP_lala", "BGNPP_shrub",
                 "FPI","FPI_P"]
sim_varname += sim_var_unobs

obs_data = get_spruce_carbonfluxes()
obs_varname = ["NCE", "Delta ABGbiomass evergreen conifer", 
               "Delta ABGbiomass deciduous conifer",
               "Delta ABGbiomass shrub", "Delta ABGnpp evergreen conifer",
               "Delta ABGnpp deciduous conifer", "Delta ABGnpp shrub",
               "ABGbiomass evergreen conifer", "ABGbiomass deciduous conifer",
               "ABGbiomass shrub", "ABGnpp evergreen conifer",
               "ABGnpp deciduous conifer", "ABGnpp shrub",
               "ANPP Tree (~48%C)", "ANPP Shrub (~50%C)",
               "BNPP Tree & Shrub", "BG_to_AG",
               "NPP Sphag.", "RHCO2"]
obs_varname += sim_var_unobs
title_list = [
    "NEE", "Delta Spruce biomass", "Delta Tamarack biomass", 
    "Delta Shrub biomass", "Delta Spruce ANPP",
    "Delta Tamarack ANPP", "Delta Shrub ANPP", 
    "Spruce biomass", "Tamarack biomass", "Shrub biomass",
    "Spruce ANPP", "Tamarack ANPP", "Shrub ANPP Salmon", 
    "Tree ANPP", "Shrub ANPP Hansen", "Tree & Shrub BNPP",
    "BNPP:ANPP tree+shrub", "Sphagnum NPP", "HR",
]
title_unobs = ["Spruce MR", "Tamarack MR", "Shrub MR", "Soil mineral N", "Soil mineral P", 
               "Spruce N limit","Tamarack N limit","Shrub N limit",
               "Spruce P limit","Tamarack P limit","Shrub P limit", # "Spruce BNPP","Tamarack BNPP","Shrub BNPP",
               "Microbe N limit", "Microbe P limit"]
title_list += title_unobs
# clist = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a']
clist = ["r", "b", "r", "b"]
mlist = ["r", "b", "none", "none"]
llist = ["-", "-", "--", "--"]


for varname, ov, title in zip(sim_varname, obs_varname, title_list):
    fig, ax = plt.subplots(figsize=(6, 4))

    count = 0
    h = [None] * 4
    for co2, CO2 in zip(["amb", "elev"], ["ACO2", "ECO2"]):
        if varname == "BG_to_AG":
            sim_temp = sim_data.loc[chamber_list[co2], "BGNPP_tree_shrub"] / \
                (sim_data.loc[chamber_list[co2], "AGNPP_pima"] + \
                 sim_data.loc[chamber_list[co2], "AGNPP_lala"] + \
                 sim_data.loc[chamber_list[co2], "AGNPP_shrub"])
        elif 'delta' in varname:
            varname2 = varname.replace('delta_', '')
            # subtract the chamber 7 baseline
            sim_temp = sim_data.loc[chamber_list[co2], varname2] - \
                sim_data.loc[("P07", 2015), varname2]
        elif varname == 'AGNPP_shrub_2':
            sim_temp = sim_data.loc[chamber_list[co2], 'AGNPP_shrub']            
        else:
            sim_temp = sim_data.loc[chamber_list[co2], varname]
            if varname in ["NEE", "HR"]:
                sim_temp = -1 * sim_temp

        sim_T = sim_data.loc[chamber_list[co2], "TBOT"]

        (h[count],) = ax.plot(
            sim_T, sim_temp, "o", color=clist[count], markerfacecolor=mlist[count]
        )
        xnew, ynew, slope, intercept, r2 = fit_line(sim_T, sim_temp)
        ax.plot(
            xnew,
            ynew,
            ls=llist[count],
            color=clist[count],
            markerfacecolor=mlist[count],
        )
        ax.text(
            0.4,
            0.65 + count * 0.08,
            "y$_{MOD\_"
            + CO2
            + "}$="
            + f"{slope:.2f}x+{intercept:.2f}  "
            + "R$^2$="
            + f"{r2:.3f}",
            color=clist[count],
            transform=ax.transAxes,
        )
        count = count + 1

        if varname in sim_var_unobs:
            continue

        if varname == "BG_to_AG":
            filt = obs_data["Plot"].isin(chamber_list[co2])
            obs_temp = obs_data.loc[filt, "BNPP Tree & Shrub"] / \
                (obs_data.loc[filt, 'ANPP Tree (~48%C)'] + \
                 obs_data.loc[filt, 'ANPP Shrub (~50%C)'])
        else:
            obs_temp = obs_data.loc[obs_data["Plot"].isin(chamber_list[co2]), ov]
        obs_T = obs_data.loc[
            obs_data["Plot"].isin(chamber_list[co2]), "Mean Annual Temp. at 2 m"
        ]

        (h[count],) = ax.plot(
            obs_T, obs_temp, "o", color=clist[count], markerfacecolor=mlist[count]
        )
        xnew, ynew, slope, intercept, r2 = fit_line(obs_T, obs_temp)
        ax.plot(
            xnew,
            ynew,
            ls=llist[count],
            color=clist[count],
            markerfacecolor=mlist[count],
        )
        ax.text(
            0.4,
            0.65 + count * 0.08,
            "y$_{OBS\_"
            + CO2
            + "}$="
            + f"{slope:.2f}x+{intercept:.2f}  "
            + "R$^2$="
            + f"{r2:.3f}",
            color=clist[count],
            transform=ax.transAxes,
        )
        count = count + 1

    ax.axhline(0.0, ls=":", color="k", lw=0.5)

    if varname in sim_var_unobs:
        ax.legend(
            h, ["MOD_ACO2", "MOD_ECO2"], ncol=4, loc=[-0.1, -0.25]
        )
    else:
        ax.legend(
            h, ["MOD_ACO2", "OBS_ACO2", "MOD_ECO2", "OBS_ECO2"], ncol=4, loc=[-0.1, -0.25]
        )
    ax.set_title(title)

    if "biomass" in title or "mineral" in title:
        ylabel = "g m-2"
    elif "N limit" in title or "P limit" in title:
        ylabel = ""
    else:
        ylabel = title.split(" ")[-1] + " gC m-2 yr-1"
    ax.set_ylabel(ylabel)
    ax.set_xlabel("Mean Annual Temperature ($^o$C)")
    fig.savefig(
        os.path.join(outdir, f"result_{prefix}_plot_{varname}.png"), dpi=600.0, bbox_inches="tight"
    )
    plt.close(fig)
