
""" Sensitivity of annual mean quantities to temperature, scatter plots """
import matplotlib.pyplot as plt
import os
import pandas as pd
from scipy.stats import linregress
import numpy as np
from utils.analysis import get_sim_carbonfluxes


def fit_line(x, y):
    filt = ~np.isnan(x) & ~np.isnan(y)
    x = x[filt]
    y = y[filt]

    if sum(filt) == 0:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    res = linregress(x, y)
    xnew = np.linspace(x.min(), x.max(), 3)
    ynew = res.slope * xnew + res.intercept
    r2 = res.rvalue**2  # coefficient of determination
    return xnew, ynew, res.slope, res.intercept, r2


#prefix = "20240311_3_2"
prefix = "UQ_20240311_3"
outdir = os.path.join(os.environ['PROJDIR'], 'ELM_Phenology', 'output', 'extract', prefix)

# T: 0, 2.25, 4.5, 6.75, 9; first ambient, then elevated CO2
chamber_list = {
    "amb": ["P06", "P20", "P13", "P08", "P17"],
    "elev": ["P19", "P11", "P04", "P16", "P10"],
}


sim_data = pd.read_csv(os.path.join(outdir, 'extract_ts_productivity.csv'), 
                       index_col=[0, 1, 2])
sim_data = sim_data.loc['average', :]
# sim_data = sim_data.loc[sim_data.index.get_level_values(0) != 2020, :]
sim_varname = ['AGBiomass_Spruce', 'AGBiomass_Tamarack', 'AGBiomass_Shrub',
               'AGNPPtoBiomass_Spruce', 'AGNPPtoBiomass_Tamarack', 'AGNPPtoBiomass_Shrub',
               'AGNPP_Spruce', 'AGNPP_Tamarack', 'AGNPP_Shrub', 'NPP_moss',
               'BGNPP_TreeShrub', 'BGtoAG_TreeShrub', 'HR', 'NEE']
units_list = ['gC m-2', 'gC m-2', 'gC m-2', 'yr-1', 'yr-1', 'yr-1',
              'gC m-2 yr-1', 'gC m-2 yr-1', 'gC m-2 yr-1', 'gC m-2 yr-1',
              'gC m-2 yr-1', '', 'gC m-2 yr-1', 'gC m-2 yr-1']
sim_var_unobs = ["SMINN_30", "SOLUTIONP_30",
                 "FPG_Spruce", "FPG_Tamarack", "FPG_Shrub", "FPG_P_Spruce", "FPG_P_Tamarack",
                 "FPG_P_Shrub", # 'BG_to_AG_pima', 'BG_to_AG_lala', 'BG_to_AG_shrub',
                 "FPI","FPI_P", "TOTSOMC",
                 "LEAFC_ALLOC_TO_TOTVEGC_ABG_Spruce",
                 "LEAFC_ALLOC_TO_TOTVEGC_ABG_Tamarack",
                 "LEAFC_ALLOC_TO_TOTVEGC_ABG_Shrub"]
units_list_unobs = ['gN m-3', 'gP m-3', '','','','','','','','','','gC m-2',"yr-1","yr-1","yr-1"]
sim_varname += sim_var_unobs
units_list = units_list + units_list_unobs


obs_data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ELM_Phenology', 'output', 'extract',
                                    'extract_obs_productivity.csv'))
obs_varname = ['AGBiomass_Spruce', 'AGBiomass_Tamarack', 'AGBiomass_Shrub',
               'AGNPPtoBiomass_Spruce', 'AGNPPtoBiomass_Tamarack', 'AGNPPtoBiomass_Shrub',
               'AGNPP_Spruce', 'AGNPP_Tamarack', 'AGNPP_Shrub', 'NPP_moss',
               'BGNPP_TreeShrub', 'BGtoAG_TreeShrub', 'HR', 'NEE']
obs_varname += sim_var_unobs
title_list = ["Spruce Biomass", "Tamarack Biomass", "Shrub Biomass",
              "Spruce Aboveground NPP:Biomass", "Tamarack Aboveground NPP:Biomass", 
              "Shrub Aboveground NPP:Biomass", 
              "Spruce Aboveground NPP", "Tamarack Aboveground NPP", "Shrub Aboveground NPP", 
              'Sphagnum NPP','Trees+Shrub Fine Root NPP', 'Trees+Shrub Root:Aboveground',
              'HR', 'NEE']
title_list += sim_var_unobs


for varname, ov, title, unit in zip(sim_varname, obs_varname, title_list, units_list):
    print(varname)
    fig, ax = plt.subplots(figsize=(6, 4))

    #if varname in ['AGNPP_tree', 'AGNPP_tree_shrub', 'AGNPP_shrub']:
    #    h = [None] * 6
    #    clist = ["r", "b", "g", "r", "b", "g"]
    #    mlist = ["r", "b", "g", "none", "none", "none"]
    #    llist = ["-", "-", "-", "--", "--", "--"]
    #else:
    h = [None] * 4
    clist = ["r", "b", "r", "b"]
    mlist = ["r", "b", "none", "none"]
    llist = ["-", "-", "--", "--"]

    count = 0

    for co2, CO2 in zip(["amb", "elev"], ["ACO2", "ECO2"]):
        #if "BG_to_AG_" in varname:
        #    pft = varname.split('_')[-1]
        #    sim_temp = sim_data.loc[chamber_list[co2], f'BGNPP_{pft}'] / \
        #               np.maximum(sim_data.loc[chamber_list[co2], f'AGNPP_{pft}'], 1e-10)
        sim_temp = sim_data.loc[chamber_list[co2], varname]
        if varname in ["NEE", "HR"]:
            sim_temp = -1 * sim_temp

        sim_T = sim_data.loc[chamber_list[co2], "Tair"]

        (h[count],) = ax.plot(
            sim_T, sim_temp, "o", color=clist[count], markerfacecolor=mlist[count]
        )
        xnew, ynew, slope, intercept, r2 = fit_line(sim_T, sim_temp)
        ax.plot(xnew, ynew, ls=llist[count], color=clist[count], markerfacecolor=mlist[count])
        ax.text(0.4, 0.65 + count * 0.08, "y$_{MOD\_" + CO2 + "}$="
            + f"{slope:.2f}x+{intercept:.2f}  "+ "R$^2$="+ f"{r2:.3f}", 
            color=clist[count], transform=ax.transAxes)
        count = count + 1

        dummy()

        if varname in sim_var_unobs:
            continue

        filt = obs_data["Plot"].isin(chamber_list[co2])
        obs_temp = obs_data.loc[filt, ov]
        obs_T = obs_data.loc[filt, "Tair"]

        (h[count],) = ax.plot(
            obs_T, obs_temp, "o", color=clist[count], markerfacecolor=mlist[count]
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
    ax.set_ylabel(unit)
    ax.set_xlabel("Mean Annual Temperature ($^o$C)")
    fig.savefig(
        os.path.join(outdir, f"result_{prefix}_plot_{varname}.png"), dpi=600.0, bbox_inches="tight"
    )
    plt.close(fig)
