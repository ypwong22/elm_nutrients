library(daymetr)
source("utils.r")

setwd("/lustre/haven/proj/UTK0134/Phenology_ELM/Daymet")


################################################################################
# Obtain the evergreen needleleaf sites' information
################################################################################
EN_sites = get_EN_sites()

metadata = jsonlite::fromJSON("https://phenocam.sr.unh.edu/webcam/network/siteinfo/")
sitelist = metadata$site
metadata = data.frame( cbind( metadata$lat, metadata$lon) )
colnames(metadata) = c("lat", "lon")
rownames(metadata) = sitelist
metadata = metadata[EN_sites, ]
metadata = metadata[!is.na(metadata["lat"]), ]
metadata$site = row.names(metadata)
row.names(metadata) = seq(1, nrow(metadata))

################################################################################
# Download the Daymet data for evergreen sites
################################################################################
download <- function(row){
    print(row)
    dir.create(file.path(".", row["site"]), showWarnings = F)

    options(show.error.messages = F)
    try( download_daymet_tiles( location = c(as.numeric(row["lat"]), as.numeric(row["lon"])),
                                start = 1980, end = 2020, 
                                path = paste0(getwd(), "/", row["site"]) ), T )
    msg = .Last.value[1]
    options(show.error.messages = T)

    if (grepl("Error", msg, fixed = T)){ 
        print(paste0( "Unsuccessful download for site ", row["site"] ))
        print(msg)
    }
}

apply( metadata, 1, download )

################################################################################
# Create the tmean variable required by phenor
################################################################################
proc_tmean <- function(row){
    path_site = paste0(getwd(), "/", row["site"])
    tile = as.numeric(strsplit(strsplit(list.files(path_site)[1], ".", fixed = T)[[1]][1], "_")[[1]][3])
    for (year in seq(1980, 2020)){
        daymet_grid_tmean(path = path_site, product = tile, year = year, internal = F)
    }
}

proc_tmean(metadata[1, ])

setwd("~/Git/phenology_elm")