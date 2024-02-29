import matplotlib.pyplot as plt
import os
import pandas as pd
from scipy.stats import linregress
import numpy as np


def fit_line(x, y):
    res = linregress(x, y)
    xnew = np.linspace(x.min(), x.max(), 3)
    ynew = res.slope * xnew + res.intercept
    r2 = res.rvalue**2  # coefficient of determination
    return xnew, ynew, res.slope, res.intercept, r2


prefix = "20231112"
#prefix = "20240227"
outdir = os.path.join(os.environ['PROJDIR'], 'ELM_Phenology', 'output', 
                      'extract', prefix)

# T: 0, 2.25, 4.5, 6.75, 9; first ambient, then elevated CO2
chamber_list = {
    "amb": ["P06", "P20", "P13", "P08", "P17"],
    "elev": ["P19", "P11", "P04", "P16", "P10"],
}
sim_data = pd.read_csv(os.path.join(outdir, 'extract_ts_productivity.csv'), 
                       index_col=[0, 1, 2])
sim_data = sim_data.loc['hummock', :] # TEMPORARY debug
sim_data = sim_data.loc[sim_data.index.get_level_values(0) != 2020, :]
sim_varname = [
    "NEE",
    "AGNPP_tree",
    "AGNPP_shrub",
    "BGNPP_tree_shrub",
    "NPP_moss",
    "HR",
]

obs_data = pd.read_excel(
    "SPRUCE C Budget Summary 28Apr2022EXP.xlsx",
    sheet_name="DataForPythonRead",
    skiprows=1,
    engine="openpyxl",
)
obs_data = obs_data.loc[obs_data["Year"] != 2020, :]
obs_varname = [
    "NCE",
    "ANPP Tree (~48%C)",
    "ANPP Shrub (~50%C)",
    "BNPP Tree & Shrub",
    "NPP Sphag.",
    "RHCO2",
]
title_list = [
    "NEE",
    "Tree ANPP",
    "Shrub ANPP",
    "Tree & Shrub BNPP",
    "Sphagnum NPP",
    "HR",
]
# clist = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a']
clist = ["b", "r", "b", "r"]
mlist = ["b", "r", "none", "none"]
llist = ["-", "-", "--", "--"]


for varname, ov, title in zip(sim_varname, obs_varname, title_list):
    fig, ax = plt.subplots(figsize=(6, 4))

    count = 0
    h = [None] * 4
    for co2, CO2 in zip(["amb", "elev"], ["ACO2", "ECO2"]):
        sim_temp = sim_data.loc[chamber_list[co2], varname]
        if varname in ["NEE", "HR"]:
            sim_temp = -1 * sim_temp

        sim_T = sim_data.loc[chamber_list[co2], "TBOT"]

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

    ax.axhline(0.0, ls=":", color="k", lw=0.5)

    ax.legend(
        h, ["OBS_ACO2", "MOD_ACO2", "OBS_ECO2", "MOD_ECO2"], ncol=4, loc=[-0.1, -0.25]
    )
    ax.set_title(title)

    ylabel = title.split(" ")[-1] + " gC m-2 yr-1"
    ax.set_ylabel(ylabel)
    ax.set_xlabel("Mean Annual Temperature ($^o$C)")
    fig.savefig(
        os.path.join(outdir, f"result_{prefix}_plot_{varname}.png"), dpi=600.0, bbox_inches="tight"
    )
    plt.close(fig)
