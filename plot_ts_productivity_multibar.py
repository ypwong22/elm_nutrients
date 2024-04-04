""" Sensitivity of annual mean quantities to temperature, bar plots 
    of slope & intercept to facilitate comparison across multiple runs """
import matplotlib.pyplot as plt
import os
import pandas as pd
import statsmodels.api as stats
import numpy as np
import seaborn as sns
from utils.analysis import get_spruce_carbonfluxes


def fit_line(x,y):
    filt = ~np.isnan(x) & ~np.isnan(y)
    x = x[filt]
    y = y[filt]
    res = stats.OLS(y, stats.add_constant(x)).fit()
    intercept = res.params[0]
    intercept_ci = res.conf_int()[0, :]
    slope = res.params[1]
    slope_ci = res.conf_int()[1, :]
    return slope, slope_ci, intercept, intercept_ci


def barplot_data(prefix_list, xlabel_list, prefix_out):
    # T: 0, 2.25, 4.5, 6.75, 9; first ambient, then elevated CO2
    chamber_list = {'amb': ['P06', 'P20', 'P13', 'P08', 'P17'], 
                    'elev': ['P19', 'P11', 'P04', 'P16', 'P10']}
    sim_varname = ["NEE", "TOTVEGC_ABG_pima", "TOTVEGC_ABG_lala", "TOTVEGC_ABG_shrub", 
                   "AGNPP_pima", "AGNPP_lala", "AGNPP_shrub", # "AGNPP_tree",
                   "BGNPP_tree_shrub", "BG_to_AG", "NPP_moss", "HR"]

    obs_data = get_spruce_carbonfluxes()
    obs_varname = ["NCE", "ABGbiomass evergreen conifer", "ABGbiomass deciduous conifer",
                   "ABGbiomass shrub", "ABGnpp evergreen conifer",
                   "ABGnpp deciduous conifer", "ABGnpp shrub", # "ANPP Shrub (~50%C)" (merged)
                   "BNPP Tree & Shrub", # "ANPP Tree (~48%C)", 
                   "BG_to_AG", "NPP Sphag.", "RHCO2"]
    title_list = ["NEE (gC m-2 yr-1)", "Spruce biomass (gC m-2)",
                  "Tamarack biomass (gC m-2)", "Shrub biomas (gC m-2)",
                  "Spruce ANPP (gC m-2 yr-1)", "Tamarack ANPP (gC m-2 yr-1)", 
                  "Shrub ANPP (gC m-2 yr-1)",  # "Tree ANPP", 
                  "Trees + shrub BNPP (gC m-2 yr-1)",
                  "BNPP:ANPP tree+shrub", "Sphagnum NPP (gC m-2 yr-1)", "HR (gC m-2 yr-1)"]

    for varname, ov, title in zip(sim_varname, obs_varname, title_list):
        slope_collection = pd.DataFrame(np.nan,
                                        index = pd.MultiIndex.from_product([['Observation'] \
                + prefix_list, ['amb', 'elev']]),
                                        columns = ['slope', 'slope_low', 'slope_up'])
        slope_collection.index.names = ['Source', 'CO2']

        mean_collection = slope_collection.copy()

        for prefix in prefix_list:

            outdir = os.path.join(os.environ['PROJDIR'], 'ELM_Phenology', 'output', 
                                'extract', prefix)
            sim_data = pd.read_csv(os.path.join(outdir, 'extract_ts_productivity.csv'), 
                                index_col=[0, 1, 2])
            sim_data = sim_data.loc['average', :]
            sim_data = sim_data.loc[sim_data.index.get_level_values(0) != 2020, :] # 2020 bad year in obs

            for co2 in ['amb','elev']:
                if varname == "BG_to_AG":
                    sim_temp = (sim_data.loc[chamber_list[co2], "BGNPP_tree_shrub"] / \
                        (sim_data.loc[chamber_list[co2], "AGNPP_pima"] + \
                        sim_data.loc[chamber_list[co2], "AGNPP_lala"] + \
                        sim_data.loc[chamber_list[co2], "AGNPP_shrub"])).values
                else:
                    sim_temp = sim_data.loc[chamber_list[co2], varname].values
                    if varname in ["NEE", "HR"]:
                        sim_temp = -1 * sim_temp
                sim_T    = sim_data.loc[chamber_list[co2], 'TBOT'].values

                filt = obs_data["Plot"].isin(chamber_list[co2])
                if varname == "BG_to_AG":
                    obs_temp = (obs_data.loc[filt, "BNPP Tree & Shrub"] / \
                        (obs_data.loc[filt, 'ANPP Tree (~48%C)'] + \
                        obs_data.loc[filt, 'ANPP Shrub (~50%C)'])).values
                elif varname == 'AGNPP_shrub':
                    obs_temp = obs_data.loc[filt, ["ANPP Shrub (~50%C)", 
                        "ABGnpp shrub"]].values.reshape(-1)
                else:
                    obs_temp = obs_data.loc[filt, ov].values
                obs_T = obs_data.loc[filt, "Mean Annual Temp. at 2 m"].values
                if varname == 'AGNPP_shrub':
                    obs_T = np.vstack([obs_T, obs_T]).T.reshape(-1)

                slope, ci, intercept, intercept_ci = fit_line(obs_T, obs_temp)
                slope_collection.loc[('Observation', co2), 'slope'] = slope
                slope_collection.loc[('Observation', co2), 'slope_low'] = ci[0]
                slope_collection.loc[('Observation', co2), 'slope_up'] = ci[1]

                # better use mean instead of intercept which is biased by the slope
                mean_collection.loc[('Observation', co2), 'mean'] = np.nanmean(obs_temp)
                mean_collection.loc[('Observation', co2), 'mean_low'] = np.nanpercentile(obs_temp, 5)
                mean_collection.loc[('Observation', co2), 'mean_up'] = np.nanpercentile(obs_temp, 95)

                slope, ci, intercept, intercept_ci = fit_line(sim_T, sim_temp)
                slope_collection.loc[(prefix, co2), 'slope'] = slope
                slope_collection.loc[(prefix, co2), 'slope_low'] = ci[0]
                slope_collection.loc[(prefix, co2), 'slope_up'] = ci[1]
                mean_collection.loc[(prefix, co2), 'mean'] = np.mean(sim_temp)
                mean_collection.loc[(prefix, co2), 'mean_low'] = np.percentile(sim_temp, 5)
                mean_collection.loc[(prefix, co2), 'mean_up'] = np.percentile(sim_temp, 95)

        slope_collection = slope_collection.reset_index()
        mean_collection = mean_collection.reset_index()

        for metric, metric_name in zip([slope_collection, mean_collection], 
                                       ['slope', 'mean']):
            fig, ax = plt.subplots(figsize = (6,4))
            sns.barplot(data = metric, x = 'Source', y = metric_name, hue = 'CO2', ax = ax,
                        hue_order = ['amb', 'elev']) # , palette = ['#b3cde3','#fbb4ae'], )

            x = [ii + jj*0.5 - 0.25 for ii in range(len(prefix_list)+1) for jj in range(2)]
            ax.errorbar(x = x, y = metric[metric_name].values, 
                        yerr = np.abs(metric[[f'{metric_name}_low', f'{metric_name}_up']
                                             ].values - metric[metric_name].values.reshape(-1,1)).T,
                        fmt = 'none', elinewidth = 1, ecolor = 'k')
            ax.set_xticklabels(["Observation"] + xlabel_list, rotation = 45)
            if metric_name == 'slope':
                ax.set_ylabel(f'{title.replace(")","")} $^o$C$^{-1}$)')
            else:
                ax.set_ylabel(title)
            ax.set_xlabel('')
            fig.savefig(os.path.join(outdir, "..", f"{prefix_out}_{metric_name}_{varname}.png"), 
                        dpi = 600., bbox_inches = 'tight')
            plt.close(fig)


if __name__ == '__main__':
    barplot_data(["20231113_3", "20240311", "20240316"], 
                 ["Default", "+ Roots uptake", "+ Roots uptake\n+ phenology"],
                 "20240401")