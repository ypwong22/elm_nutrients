""" Sensitivity of annual mean quantities to temperature, bar plots 
    of slope & intercept to facilitate comparison across multiple runs """
import matplotlib.pyplot as plt
import os
import pandas as pd
import statsmodels.api as stats
import numpy as np
import seaborn as sns
from matplotlib import rcParams


rcParams['font.size'] = 16
rcParams['axes.titlesize'] = 16


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
    sim_varname = ['AGBiomass_Spruce', 'AGBiomass_Tamarack', 'AGBiomass_Shrub',
                   'AGNPPtoBiomass_Spruce', 'AGNPPtoBiomass_Tamarack', 'AGNPPtoBiomass_Shrub',
                   'AGNPP_Spruce', 'AGNPP_Tamarack', 'AGNPP_Shrub', 'NPP_moss',
                   'BGNPP_TreeShrub', 'BGtoAG_TreeShrub', 'HR', 'NEE', 'AGNPP_Tree']

    obs_data = pd.read_csv(os.path.join(os.environ['PROJDIR'], 'ELM_Nutrients', 'output', 'extract',
                                        'extract_obs_productivity.csv'))
    obs_varname = ['AGBiomass_Spruce', 'AGBiomass_Tamarack', 'AGBiomass_Shrub',
                   'AGNPPtoBiomass_Spruce', 'AGNPPtoBiomass_Tamarack', 'AGNPPtoBiomass_Shrub',
                   'AGNPP_Spruce', 'AGNPP_Tamarack', 'AGNPP_Shrub', 'NPP_moss',
                   'BGNPP_TreeShrub', 'BGtoAG_TreeShrub', 'HR', 'NEE', 'AGNPP_Tree']
    title_list = ["Spruce Biomass", "Tamarack Biomass", "Shrub Biomass",
                  "Spruce ANPP:Biomass", "Tamarack ANPP:Biomass", 
                  "Shrub ANPP:Biomass", 
                  "Spruce ANPP", "Tamarack ANPP", "Shrub ANPP", 
                  'Sphagnum NPP','Trees+Shrub Fine Root NPP', 'Trees+Shrub Root:Aboveground',
                  'HR', 'NEE', 'Tree ANPP']
    units_list = ['gC m-2', 'gC m-2', 'gC m-2', 'yr-1', 'yr-1', 'yr-1',
                  'gC m-2 yr-1', 'gC m-2 yr-1', 'gC m-2 yr-1', 'gC m-2 yr-1',
                  'gC m-2 yr-1', '', 'gC m-2 yr-1', 'gC m-2 yr-1', 'gC m-2 yr-1']

    for varname, ov, title, unit in zip(sim_varname, obs_varname, title_list, units_list):
        slope_collection = pd.DataFrame(np.nan,
                                        index = pd.MultiIndex.from_product([['Observation'] \
                + prefix_list, ['Ambient', 'Elevated']]),
                                        columns = ['slope', 'slope_low', 'slope_up'])
        slope_collection.index.names = ['Source', 'CO2']

        mean_collection = slope_collection.copy()

        for prefix in prefix_list:

            outdir = os.path.join(os.environ['PROJDIR'], 'ELM_Nutrients', 'output', 
                                  'extract', prefix)
            sim_data = pd.read_csv(os.path.join(outdir, 'extract_ts_productivity.csv'), 
                                   index_col=[0, 1, 2])
            sim_data = sim_data.loc['average', :]
            sim_data = sim_data.loc[sim_data.index.get_level_values(0) != 2020, :] # 2020 bad year in obs

            for co2,co2_ in zip(['amb','elev'], ['Ambient','Elevated']):
                if varname == 'AGNPP_Tree':
                    sim_temp = sim_data.loc[chamber_list[co2], 'AGNPP_Spruce'].values + \
                        sim_data.loc[chamber_list[co2], 'AGNPP_Tamarack'].values
                else:
                    sim_temp = sim_data.loc[chamber_list[co2], varname].values
                if varname in ["NEE", "HR"]:
                    sim_temp = -1 * sim_temp

                sim_T = sim_data.loc[chamber_list[co2], "Tair"].values

                filt = obs_data["Plot"].isin(chamber_list[co2])
                if ov == 'AGNPP_Tree':
                    obs_temp = obs_data.loc[filt, 'AGNPP_Spruce'].values + \
                        obs_data.loc[filt, 'AGNPP_Tamarack'].values
                else:
                    obs_temp = obs_data.loc[filt, ov].values
                obs_T = obs_data.loc[filt, "Tair"].values

                slope, ci, intercept, intercept_ci = fit_line(obs_T, obs_temp)
                slope_collection.loc[('Observation', co2_), 'slope'] = slope
                slope_collection.loc[('Observation', co2_), 'slope_low'] = ci[0]
                slope_collection.loc[('Observation', co2_), 'slope_up'] = ci[1]

                # better use mean instead of intercept which is biased by the slope
                mean_collection.loc[('Observation', co2_), 'mean'] = np.nanmean(obs_temp)
                mean_collection.loc[('Observation', co2_), 'mean_low'] = np.nanpercentile(obs_temp, 5)
                mean_collection.loc[('Observation', co2_), 'mean_up'] = np.nanpercentile(obs_temp, 95)

                print(prefix, varname)

                slope, ci, intercept, intercept_ci = fit_line(sim_T, sim_temp)
                slope_collection.loc[(prefix, co2_), 'slope'] = slope
                slope_collection.loc[(prefix, co2_), 'slope_low'] = ci[0]
                slope_collection.loc[(prefix, co2_), 'slope_up'] = ci[1]
                mean_collection.loc[(prefix, co2_), 'mean'] = np.mean(sim_temp)
                mean_collection.loc[(prefix, co2_), 'mean_low'] = np.percentile(sim_temp, 5)
                mean_collection.loc[(prefix, co2_), 'mean_up'] = np.percentile(sim_temp, 95)

        slope_collection = slope_collection.reset_index()
        mean_collection = mean_collection.reset_index()

        for metric, metric_name in zip([slope_collection, mean_collection], 
                                       ['slope', 'mean']):
            fig, ax = plt.subplots(figsize = (6,6))
            sns.barplot(data = metric, x = 'Source', y = metric_name, hue = 'CO2', ax = ax,
                        hue_order = ['Ambient', 'Elevated'], palette = ['#377eb8','#ff7f00'])

            x = [ii + jj*0.5 - 0.25 for ii in range(len(prefix_list)+1) for jj in range(2)]

            ax.errorbar(x = x, y = metric[metric_name].values, 
                        yerr = np.abs(metric[[f'{metric_name}_low', f'{metric_name}_up']
                                             ].values - metric[metric_name].values.reshape(-1,1)).T,
                        fmt = 'none', elinewidth = 1, ecolor = 'k')
            ax.set_xticklabels(["Observation"] + xlabel_list) # , rotation = 45)
            ax.set_title(f'{title}')
            if metric_name == 'slope':
                ax.set_ylabel(f'Temperature sensitivity ({unit} $^o$C$^{-1}$)')
            else:
                ax.set_ylabel(f'Mean ({unit})')
            ax.set_xlabel('')
            fig.savefig(os.path.join(outdir, "..", f"{prefix_out}_{metric_name}_{varname}.png"), 
                        dpi = 600., bbox_inches = 'tight')
            plt.close(fig)


if __name__ == '__main__':
    barplot_data(["20231113", "UQ_20231113", "20240313", "20240314", "20240315"], 
                 ["Default", "Default-opt", "v1", "v2", "v3"],
                 "20241014")
    #barplot_data(["20231113_4", "20240316_1"], 
    #             ["Default Model", "New Model"], "20240429")
