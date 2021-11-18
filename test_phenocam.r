library(phenor)
source("phenor-master/R/pr_fm_phenocam.R") # with debug changes
source("utils/utils.r")

path_data = "/lustre/haven/proj/UTK0134/DATA/Vegetation/PhenoCam_V2_1674/"
path_out  = "/lustre/haven/proj/UTK0134/Phenology_ELM/output_phenor/"

# MODIFY
direction = "falling"

####################################################################################################
# Load formated phenocam data
####################################################################################################
read_anew = T

if (read_anew){
    path_transition = paste0(path_data, "data/data_record_5/3day")
    phenocam_data = pr_fm_phenocam(path = path_transition, direction = direction, 
                                   gcc_value = "gcc_90", threshold = 25, offset = 264, 
                                   spread = 30, internal = T)
    save.image(paste0(path_out, "phenocam_data_", direction, ".RData"))
} else {
    load(paste0(path_out, "phenocam_data_", direction, ".RData"))
}

####################################################################################################
# Extract the phenology for a fit 
####################################################################################################
plot_gcc_series <- function(site, result_model){
    path0                  = paste0(path_data, "data/data_record_4/")
    tsfiles                = list.files(path0, glob2rx( paste0(site, "_EN_*_3day.csv") ))
    csv_list               = lapply(paste0(path0, tsfiles), readLines)
    for (i in seq(1, length(csv_list))){
        if (i == 1){
            csv_list[[1]]          = csv_list[[1]][25:length(csv_list[[1]])]
        } else {
            csv_list[[i]]     = csv_list[[i]][26:length(csv_list[[i]])]
        }
    }
    csv_list               = unlist(csv_list)
    writeLines(text        = csv_list, con = "test_phenocam_temp.csv")
    series                 = read.csv("test_phenocam_temp.csv")

    pdf(paste0(path_out, paste0(direction, "_", site, "_", result_model$model, ".pdf")))
    par(mfrow = c(2,1))
    plot(as.Date(series[["date"]]), series[["smooth_gcc_90"]])

    plot(series[["doy"]], series[["smooth_gcc_90"]])
    abline(v = result_model$predicted, col = "red")
    abline(v = result_model$measured, col = "green")

    legend("bottomright", c("predicted transition dates",
                            "measured transition dates"), pch=1, lty=1,
           inset=c(0,1), xpd=T, horiz=T, bty="n", 
           col = c("red", "green"))
    dev.off()
}

EN_sites = get_EN_sites()
print(paste0("Number of EN sites: ", length(EN_sites)))
model_list = c("LIN","TT","TTs","PTT","PTTs",
               "M1","M1s","AT","SQ","SQb","SM1",
               "SM1b","PA","PAb","PM1",
               "PM1b","UN","UM1","SGSI","AGSI")

for (i in seq(2, length(EN_sites))){
    site = EN_sites[[i]]

    for (j in seq(1, length(model_list))){
        model = model_list[j]

        if (is.null(phenocam_data[[site]])){ next }

        result_model = pr_fit(model, phenocam_data[[site]], control = list(max.call = 40000))

        # Plot the goodness-of-fit of the extracted phenological dates
        plot_gcc_series(site, result_model)
    }
}

####################################################################################################
# Plot the goodness-of-fit of the extracted phenological dates
####################################################################################################
results = pr_fit_comparison(data = phenocam_data[site], model = model_list,  method = "GenSA", 
                            control = list(max.call = 40000))
pr_plot_arrows(results, models = model_list)
dev.off()
pr_plot_comparison(results)
dev.off()
