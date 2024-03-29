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
    intercept = res.params.values[0]
    intercept_ci = res.conf_int().values[0, :]
    slope = res.params.values[1]
    slope_ci = res.conf_int().values[1, :]
    return slope, slope_ci, intercept, intercept_ci


def barplot_data(prefix_list, prefix_out):
    # T: 0, 2.25, 4.5, 6.75, 9; first ambient, then elevated CO2
    chamber_list = {'amb': ['P06', 'P20', 'P13', 'P08', 'P17'], 'elev': ['P19', 'P11', 'P04', 'P16', 'P10']}
    sim_varname = ["NEE", "TOTVEGC_ABG_pima", "TOTVEGC_ABG_lala", "TOTVEGC_ABG_shrub", 
                   "AGNPP_pima", "AGNPP_lala", "AGNPP_shrub", "AGNPP_tree", "AGNPP_shrub_2", 
                   "BGNPP_tree_shrub", "BG_to_AG", "NPP_moss", "HR"]

    obs_data = get_spruce_carbonfluxes()
    obs_varname = ["NCE", "ABGbiomass evergreen conifer", "ABGbiomass deciduous conifer",
                   "ABGbiomass shrub", "ABGnpp evergreen conifer",
                   "ABGnpp deciduous conifer", "ABGnpp shrub",
                   "ANPP Tree (~48%C)", "ANPP Shrub (~50%C)", "BNPP Tree & Shrub", 
                   "BG_to_AG", "NPP Sphag.", "RHCO2"]
    title_list = ["NEE", "Spruce biomass", "Tamarack biomass", "Shrub biomass",
                  "Spruce ANPP", "Tamarack ANPP", "Shrub ANPP Salmon", 
                  "Tree ANPP", "Shrub ANPP Hansen", "Tree & Shrub BNPP",
                  "BNPP:ANPP tree+shrub", "Sphagnum NPP", "HR"]

    for varname, ov, title in zip(sim_varname, obs_varname, title_list):
        slope_collection = pd.DataFrame(np.nan,
                                        index = pd.MultiIndex.from_product([['Observation'] \
                + prefix_list, ['amb', 'elev']]),
                                        columns = ['slope', 'slope_low', 'slope_up'])
        slope_collection.index.names = ['Source', 'CO2']

        intercept_collection = slope_collection.copy()

        for prefix in prefix_list:

            outdir = os.path.join(os.environ['PROJDIR'], 'ELM_Phenology', 'output', 
                                'extract', prefix)
            sim_data = pd.read_csv(os.path.join(outdir, 'extract_ts_productivity.csv'), 
                                index_col=[0, 1, 2])
            sim_data = sim_data.loc['average', :]
            sim_data = sim_data.loc[sim_data.index.get_level_values(0) != 2020, :] # 2020 bad year in obs

            count = 0
            for co2 in ['amb','elev']:
                if varname == "BG_to_AG":
                    sim_temp = sim_data.loc[chamber_list[co2], "BGNPP_tree_shrub"] / \
                        (sim_data.loc[chamber_list[co2], "AGNPP_pima"] + \
                        sim_data.loc[chamber_list[co2], "AGNPP_lala"] + \
                        sim_data.loc[chamber_list[co2], "AGNPP_shrub"])
                elif varname == 'AGNPP_shrub_2':
                    sim_temp = sim_data.loc[chamber_list[co2], 'AGNPP_shrub']            
                else:
                    sim_temp = sim_data.loc[chamber_list[co2], varname]
                    if varname in ["NEE", "HR"]:
                        sim_temp = -1 * sim_temp

                sim_T    = sim_data.loc[chamber_list[co2], 'TBOT']

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

                slope, ci, intercept, intercept_ci = fit_line(obs_T, obs_temp)
                slope_collection.loc[('Observation', co2), 'slope'] = slope
                slope_collection.loc[('Observation', co2), 'slope_low'] = ci[0]
                slope_collection.loc[('Observation', co2), 'slope_up'] = ci[1]
                intercept_collection.loc[('Observation', co2), 'intercept'] = intercept
                intercept_collection.loc[('Observation', co2), 'intercept_low'] = intercept_ci[0]
                intercept_collection.loc[('Observation', co2), 'intercept_up'] = intercept_ci[1]

                slope, ci, intercept, intercept_ci = fit_line(sim_T, sim_temp)
                slope_collection.loc[(prefix, co2), 'slope'] = slope
                slope_collection.loc[(prefix, co2), 'slope_low'] = ci[0]
                slope_collection.loc[(prefix, co2), 'slope_up'] = ci[1]
                intercept_collection.loc[(prefix, co2), 'intercept'] = intercept
                intercept_collection.loc[(prefix, co2), 'intercept_low'] = intercept_ci[0]
                intercept_collection.loc[(prefix, co2), 'intercept_up'] = intercept_ci[1]

        slope_collection = slope_collection.reset_index()
        intercept_collection = intercept_collection.reset_index()

        for metric, metric_name in zip([slope_collection, intercept_collection], 
                                       ['slope', 'intercept']):
            fig, ax = plt.subplots(figsize = (6,4))
            sns.barplot(data = metric, x = 'Source', y = metric_name, hue = 'CO2', ax = ax,
                        hue_order = ['amb', 'elev'], palette = ['#66c2a5', '#fc8d62', '#8da0cb'])

            x = [ii + jj*0.5 - 0.25 for ii in range(len(prefix_list)+1) for jj in range(2)]
            ax.errorbar(x = x, y = metric[metric_name].values, 
                        yerr = np.abs(metric[[f'{metric_name}_low', f'{metric_name}_up']
                                             ].values - metric[metric_name].values.reshape(-1,1)).T,
                        fmt = 'none', elinewidth = 1, ecolor = 'k')

            ax.set_xticklabels(ax.get_xticklabels(), rotation = 45)

            ylabel = title.split(' ')[-1] + ' gC m-2 yr-1 $\degree$C$^{-1}$'
            ax.set_ylabel(ylabel)
            ax.set_xlabel('Mean Annual Temperature ($^o$C)')
            fig.savefig(os.path.join(outdir, "..", f"{prefix_out}_{metric_name}_{varname}.png"), 
                        dpi = 600., bbox_inches = 'tight')
            plt.close(fig)


if __name__ == '__main__':

    barplot_data(["20231112", "20240317_5", "20240316"], "20240327")