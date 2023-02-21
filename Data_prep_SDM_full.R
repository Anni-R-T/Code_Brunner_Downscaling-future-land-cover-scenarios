if(Sys.info()[["nodename"]]=="vmdomisch04") DIR="/mnt/brunner/brunner/R" 
setwd(DIR)
R_temp_delete <- paste0(DIR, "/R_temp_delete")
rasterOptions(tmpdir=R_temp_delete)

rasterOptions()

if (!require("raster")) { install.packages("raster", dependencies = TRUE) ; library(raster)}
if (!require("rgdal")) { install.packages("rgdal", dependencies = TRUE) ; library(rgdal)}
if (!require("maptools")) { install.packages("maptools", dependencies = TRUE) ; library(maptools)}
if (!require("rgeos")) { install.packages("rgeos", dependencies = TRUE) ; library(rgeos)}
if (!require("foreach")) { install.packages("foreach", dependencies = TRUE) ; library(foreach)}
if (!require("doParallel")) { install.packages("doParallel", dependencies = TRUE) ; library(doParallel)}
if (!require("plyr")) { install.packages("plyr", dependencies = TRUE) ; library(plyr)}
if (!require("dplyr")) { install.packages("dplyr", dependencies = TRUE) ; library(dplyr)}
if (!require("reshape")) { install.packages("reshape", dependencies = TRUE) ; library(reshape)}
if (!require("data.table")) { install.packages("data.table", dependencies = TRUE) ; library(data.table)}
if (!require("biomod2")) { BiocManager::install("biomod2") ; library(biomod2)}
if (!require("viridis")) { install.packages("viridis", dependencies = TRUE) ; library(viridis)}
if (!require("RColorBrewer")) { install.packages("RColorBrewer", dependencies = TRUE) ; library(RColorBrewer)}
if (!require("ggplot2")) { install.packages("ggplot2", dependencies = TRUE) ; library(ggplot2)}
if (!require("ggthemes")) { install.packages("ggthemes", dependencies = TRUE) ; library(ggthemes)} # theme_map()
if (!require("caret")) { install.packages("caret", dependencies = TRUE) ; library(caret)} # theme_map()
if (!require("scales")) { install.packages("scales", dependencies = TRUE) ; library(scales)} # scaling

#=======================================load study_area shape and basins===============================#
study_area <- shapefile("/mnt/brunner/brunner/stream_network/basins_study_area_outline.shp")

#load raster of basin
basin <- raster("/mnt/brunner/brunner/stream_network/basins_study_area.tif")
#basin <- crop(basin, extent(study_area))
#basin <- mask(basin, study_area)
#x11();plot(basin)




#========================1: present climate, topographical and 2010 Landcover===================================#
#========================load maps
#TOPOGRAPHICAL DATA
topo_global_list <- list.files(path = "/mnt/brunner/brunner/SDM_data/topo/", 
                               pattern = "*.tif$", 
                               full.names = TRUE)
topo_global <- stack(topo_global_list)

topography <- crop(topo_global, extent(basin))
topography<- mask(topography, basin)


#BIOCLIM PRESENT DATA
bioclim_global_list <- list.files(path = "/mnt/brunner/brunner/SDM_data/Bioclim_present/", 
                                  pattern = "*.tif$", 
                                  full.names = TRUE)
bioclim_global <- stack(bioclim_global_list)

bioclim_present <- crop(bioclim_global, extent(basin))
bioclim_present <- mask(bioclim_present, basin)
#NAvalue(bioclim_present)<-(-9999)

#========================calculate statistics per basin

#MEAN
bioclim_present_mean <- zonal(bioclim_present, basin, fun = mean, digits = 0, na.rm = TRUE)
bioclim_present_mean <-as.data.frame(bioclim_present_mean)
colnames(bioclim_present_mean)[1] <- "basin_id" #rename zone to basin_id to fit ensemble model
colnames(bioclim_present_mean)[2:10] <- paste("mean", colnames(bioclim_present_mean[,c(2:10)]), sep = "_")
write.csv(bioclim_present_mean, "/mnt/brunner/brunner/SDM_data/tables/bioclim_present_mean.csv", row.names = F)

topography_mean <- zonal(topography, basin, fun = mean, digits = 0, na.rm = TRUE)
topography_mean <-as.data.frame(topography_mean)
colnames(topography_mean)[1] <- "basin_id"
colnames(topography_mean)[2:6] <- paste("mean", colnames(topography_mean[,c(2:6)]), sep = "_")
write.csv(topography_mean, "/mnt/brunner/brunner/SDM_data/tables/topography_mean.csv", row.names = F)

#MIN#
bioclim_present_min <- zonal(bioclim_present, basin, fun = min, digits = 0, na.rm = TRUE)
bioclim_present_min <-as.data.frame(bioclim_present_min)
colnames(bioclim_present_min)[1] <- "basin_id"
colnames(bioclim_present_min)[2:10] <- paste("min", colnames(bioclim_present_min[,c(2:10)]), sep = "_")

topography_min <- zonal(topography, basin, fun = min, digits = 0, na.rm = TRUE)
topography_min <-as.data.frame(topography_min)
colnames(topography_min)[1] <- "basin_id"
colnames(topography_min)[2:6] <- paste("min", colnames(topography_min[,c(2:6)]), sep = "_")


#MAX#
bioclim_present_max <- zonal(bioclim_present, basin, fun = max, digits = 0, na.rm = TRUE)
bioclim_present_max  <-as.data.frame(bioclim_present_max)
colnames(bioclim_present_max)[1] <- "basin_id"
colnames(bioclim_present_max)[2:10] <- paste("max", colnames(bioclim_present_max[,c(2:10)]), sep = "_")

topography_max <- zonal(topography, basin, fun = max, digits = 0, na.rm = TRUE)
topography_max  <-as.data.frame(topography_max )
colnames(topography_max)[1] <- "basin_id"
colnames(topography_max)[2:6] <- paste("max", colnames(topography_max[,c(2:6)]), sep = "_")

#RANGE
bioclim_present_range <- bioclim_present_max %>%  dplyr::mutate(range_ann_mean_temp = bioclim_present_max[, 2] - bioclim_present_min[, 2], 
                                                                range_ann_prec = bioclim_present_max[, 3] - bioclim_present_min[, 3],
                                                                range_mean_temp_coldest_quarter = bioclim_present_max[, 4] - bioclim_present_min[, 4],
                                                                range_mean_temp_warmest_quarter = bioclim_present_max[, 5] - bioclim_present_min[, 5],
                                                                range_prec_driest_quarter = bioclim_present_max[, 6] - bioclim_present_min[, 6], 
                                                                range_prec_seasonality = bioclim_present_max[, 7] - bioclim_present_min[, 7],
                                                                range_prec_wettest_quarter = bioclim_present_max[, 8] - bioclim_present_min[, 8],
                                                                #temp_range = bioclim_present_max[, 9] - bioclim_present_min[, 9],
                                                                range_temp_seasonality = bioclim_present_max[, 10] - bioclim_present_min[, 10])
bioclim_present_range <- bioclim_present_range[-c(2:10)]
write.csv(bioclim_present_range, "/mnt/brunner/brunner/SDM_data/tables/bioclim_present_range.csv", row.names = F)


topography_range <- topography_max %>%  dplyr::mutate(range_eastness  = topography_max[, 2] - topography_min[, 2], 
                                                      range_elevation = topography_max[, 3] - topography_min[, 3],
                                                      range_northness = topography_max[, 4] - topography_min[, 4],
                                                      range_roughness = topography_max[, 5] - topography_min[, 5],
                                                      range_slope = topography_max[, 6] - topography_min[, 6]) 

topography_range <- topography_range[-c(2:6)]
write.csv(topography_range, "/mnt/brunner/brunner/SDM_data/tables/topography_range.csv", row.names = F)

#add tables together
topography_bioclim_present <- Reduce(merge, list(bioclim_present_mean, bioclim_present_range, topography_mean, topography_range))
write.csv(topography_bioclim_present, "/mnt/brunner/brunner/SDM_data/tables/topography_bioclim_presen.csv", row.names = F)


#============LC percent per basin

#LANDCOVER ESA 2010
LC_ESA_2010_LAC <- raster("/mnt/brunner/brunner/maps/1km/LCLA_2010_1km.tif")
LC_ESA_2010 <- crop(LC_ESA_2010_LAC, extent(study_area))
LC_ESA_2010 <- mask(LC_ESA_2010, study_area)
extent(LC_ESA_2010) <- alignExtent(LC_ESA_2010, basin)
NAvalue(LC_ESA_2010)<--9999
#x11(); plot(LC_ESA_2010)

# my zonal function to get number of pixels per LC category per basin
myZonal <- function (x, z, digits = 0, na.rm = TRUE) { 
  vals <- raster::getValues(x) 
  zones <- round(raster::getValues(z), digits = digits) 
  rDT <- data.table::data.table(vals, z = zones) 
  plyr::count(rDT) }

# pixelcounts for full basin
zstat_basin <- myZonal(basin, basin)
results_basin <- subset(zstat_basin, select =c("vals", "freq"))
colnames(results_basin)[1] <- "basin_id"
results_basin[is.na(results_basin)] <- 0
write.csv(results_basin,"/mnt/brunner/brunner/SDM_data/tables/pixels_basin.csv", row.names = F )

##pixelcount for the ESA map categories per basin
ESA_2010_zonal <-myZonal(LC_ESA_2010, basin)
results_10 <- data.table::dcast(ESA_2010_zonal, z ~ vals)
colnames(results_10)[1] <- "basin_id"
results_10[is.na(results_10)] <- 0

#bind basin count and LC 10 count
#all_results <- cbind(results_basin, results_10 )
basin_2010_LC <- merge(results_basin, results_10, by = "basin_id", all= T )
basin_2010_LC <- basin_2010_LC[-1,]



LC_percentage_basin_2010 <- (basin_2010_LC[,3:8]/basin_2010_LC[,2])
LC_percentage_basin_2010 <- round(LC_percentage_basin_2010, digits = 5)
LC_percentage_basin_2010 <- cbind(basin_2010_LC[1], LC_percentage_basin_2010)

LC_percentage_basin_2010 <- as.data.frame(LC_percentage_basin_2010)
setnames(LC_percentage_basin_2010, old = c('1', '2', '3', '4', '5', '6'), new = c('crops', 'forest','pasture','urban', 'other', 'water'))
write.csv(LC_percentage_basin_2010,"/mnt/brunner/brunner/SDM_data/tables/LC_percentage_2010.csv", row.names = F )

topography_bioclim_landcover_present <- Reduce(merge, list(bioclim_present_mean, bioclim_present_range, topography_mean, topography_range, LC_percentage_basin_2010))
write.csv(topography_bioclim_present, "/mnt/brunner/brunner/SDM_data/tables/topography_bioclim_landcover_present.csv", row.names = F)

#====================stream_topo
#load the data from the GRASS stream vector
stream_vector_data <- read.csv("/mnt/brunner/brunner/SDM_data/tables/stream_vector_data_study_area.csv")

#make list with variables needed for the model 8from ensemble example)
get_these <- c("basin_id",  # rename to "basin_id"
               "strahler", # Strahler's stream order
               "length", # stream length
               "stright", # length of stream as stright line
               "sinosoid", # fractal dimension: stream length/stright stream length
               "cum_length", # length of stream from source
               "flow_accum", # upstream contributing area
               "out_dist",  #  distance of current stream init from outlet
               "source_elev", # elevation of stream init
               "outlet_elev", #elevation of stream outlet
               "elev_drop", # difference between source_elev and outlet_elev + drop outlet
               "out_drop", # drop at the outlet of the stream
               "gradient") # drop/length

#reduce the stream vector data to those variables needed
stream_topo <- stream_vector_data[get_these]
#add togeter topo bio lc and the stream vector data
#topography_bioclim_landcover_present <- read.csv("/mnt/brunner/brunner/SDM_data/tables/topo_bio_lc_present.csv")
topo_bio_LC_stream_present <- merge (topography_bioclim_landcover_present, stream_topo, by = "basin_id")
write.csv(topo_bio_LC_stream_present, "/mnt/brunner/brunner/SDM_data/tables/topo_bio_lc_stream_present.csv", row.names = F)




#========================check_correlation


corr_topo_bioclim_LC <- cor(topo_bio_LC_stream_present[2:46])
corr_topo_bioclim_LC<- as.data.frame(corr_topo_bioclim_LC)
write.csv(corr_topo_bioclim_LC, "/mnt/brunner/brunner/SDM_data/tables/cor_topography_bioclim_landcover_present.csv", )

#install.packages("Hmisc")
library("Hmisc")
rcorr_topo_bioclim_LC <- rcorr(as.matrix(topo_bio_LC_stream_present[2:46]))

rcorr_topo_bioclim_LC_coeff <- rcorr_topo_bioclim_LC$r
rcorr_topo_bioclim_LC_coeff <- as.data.frame(rcorr_topo_bioclim_LC_coeff)
write.csv(rcorr_topo_bioclim_LC_coeff, "/mnt/brunner/brunner/SDM_data/tables/rcorr_topo_bio_LC_stream_present.csv" )

rcorr_topo_bioclim_LC_p = rcorr_topo_bioclim_LC$P
rcorr_topo_bioclim_LC_p <- as.data.frame(rcorr_topo_bioclim_LC_p)
write.csv(rcorr_topo_bioclim_LC_p, "/mnt/brunner/brunner/SDM_data/tables/pvalue_topo_bio_LC_stream_present.csv")

#install.packages("corrplot")
library(corrplot)
x11();corrplot(corr_topo_bioclim_LC)

#========================2: climate, topographical and Landcover 2050 ===================================#

#RCP 6

#=============================== BIOCLIM CCSM rcp6 2041-2060
CCSM_60_2050_list <- list.files(path = "/mnt/brunner/brunner/SDM_data/Bioclim_2050/New_NA/CCSM_global_60/", 
                                pattern = "*.tif$", 
                                full.names = TRUE)
CCSM_60_2050_global<- stack(CCSM_60_2050_list)
CCSM_60_2050 <- crop(CCSM_60_2050_global, extent(basin))
CCSM_60_2050 <- mask(CCSM_60_2050, basin)
x11(); plot(CCSM_60_2050)

CCSM_60_2050_mean <- zonal(CCSM_60_2050, basin, fun = mean, digits = 0, na.rm = TRUE)
CCSM_60_2050_mean <-as.data.frame(CCSM_60_2050_mean)
colnames(CCSM_60_2050_mean)[1] <- "basin_id" #rename zone to basin_id to fit ensemble model
colnames(CCSM_60_2050_mean)[2:10] <- paste("mean", colnames(CCSM_60_2050_mean[,c(2:10)]), sep = "_")
write.csv(CCSM_60_2050_mean, "/mnt/brunner/brunner/SDM_data/tables/CCSM_60_2050_mean.csv", row.names = F)

CCSM_60_2050_min <- zonal(CCSM_60_2050, basin, fun = min, digits = 0, na.rm = TRUE)
CCSM_60_2050_min <-as.data.frame(CCSM_60_2050_min)
colnames(CCSM_60_2050_min)[1] <- "basin_id"
colnames(CCSM_60_2050_min)[2:10] <- paste("min", colnames(CCSM_60_2050_min[,c(2:10)]), sep = "_")

CCSM_60_2050_max <- zonal(CCSM_60_2050, basin, fun = max, digits = 0, na.rm = TRUE)
CCSM_60_2050_max  <-as.data.frame(CCSM_60_2050_max)
colnames(CCSM_60_2050_max)[1] <- "basin_id"
colnames(CCSM_60_2050_max)[2:10] <- paste("max", colnames(CCSM_60_2050_max[,c(2:10)]), sep = "_")

CCSM_60_2050_range <- CCSM_60_2050_max %>%  dplyr::mutate(range_CCSM_60_2050_ann_mean_temp = CCSM_60_2050_max[, 2] - CCSM_60_2050_min[, 2], 
                                                                range_CCSM_60_2050_ann_prec = CCSM_60_2050_max[, 3] - CCSM_60_2050_min[, 3],
                                                                range_CCSM_60_2050_mean_temp_coldest_quarter = CCSM_60_2050_max[, 4] - CCSM_60_2050_min[, 4],
                                                                range_CCSM_60_2050_mean_temp_warmest_quarter = CCSM_60_2050_max[, 5] - CCSM_60_2050_min[, 5],
                                                                range_CCSM_60_2050_prec_driest_quarter = CCSM_60_2050_max[, 6] - CCSM_60_2050_min[, 6], 
                                                                range_CCSM_60_2050_prec_seasonality = CCSM_60_2050_max[, 7] - CCSM_60_2050_min[, 7],
                                                                range_CCSM_60_2050_prec_wettest_quarter = CCSM_60_2050_max[, 8] - CCSM_60_2050_min[, 8],
                                                                #temp_range = CCSM_60_2050_max[, 9] - CCSM_60_2050_min[, 9],
                                                                range_CCSM_60_2050_temp_seasonality = CCSM_60_2050_max[, 10] - CCSM_60_2050_min[, 10])
CCSM_60_2050_range <- CCSM_60_2050_range[-c(2:10)]
write.csv(CCSM_60_2050_range, "/mnt/brunner/brunner/SDM_data/tables/CCSM_60_2050_range.csv", row.names = F)

#=============================== BIOCLIM IPSL rcp6 2041-2060
IPSL_60_2050_list <- list.files(path = "/mnt/brunner/brunner/SDM_data/Bioclim_2050/New_NA/IPSL_global_60/", 
                                pattern = "*.tif$", 
                                full.names = TRUE)
IPSL_60_2050_global<- stack(IPSL_60_2050_list)
IPSL_60_2050 <- crop(IPSL_60_2050_global, extent(basin))
IPSL_60_2050 <- mask(IPSL_60_2050, basin)

IPSL_60_2050_mean <- zonal(IPSL_60_2050, basin, fun = mean, digits = 0, na.rm = TRUE)
IPSL_60_2050_mean <-as.data.frame(IPSL_60_2050_mean)
colnames(IPSL_60_2050_mean)[1] <- "basin_id" #rename zone to basin_id to fit ensemble model
colnames(IPSL_60_2050_mean)[2:10] <- paste("mean", colnames(IPSL_60_2050_mean[,c(2:10)]), sep = "_")
write.csv(IPSL_60_2050_mean, "/mnt/brunner/brunner/SDM_data/tables/IPSL_60_2050_mean.csv", row.names = F)

IPSL_60_2050_min <- zonal(IPSL_60_2050, basin, fun = min, digits = 0, na.rm = TRUE)
IPSL_60_2050_min <-as.data.frame(IPSL_60_2050_min)
colnames(IPSL_60_2050_min)[1] <- "basin_id"
colnames(IPSL_60_2050_min)[2:10] <- paste("min", colnames(IPSL_60_2050_min[,c(2:10)]), sep = "_")

IPSL_60_2050_max <- zonal(IPSL_60_2050, basin, fun = max, digits = 0, na.rm = TRUE)
IPSL_60_2050_max  <-as.data.frame(IPSL_60_2050_max)
colnames(IPSL_60_2050_max)[1] <- "basin_id"




colnames(IPSL_60_2050_max)[2:10] <- paste("max", colnames(IPSL_60_2050_max[,c(2:10)]), sep = "_")

IPSL_60_2050_range <- IPSL_60_2050_max %>%  dplyr::mutate(range_IPSL_60_2050_ann_mean_temp = IPSL_60_2050_max[, 2] - IPSL_60_2050_min[, 2], 
                                                          range_IPSL_60_2050_ann_prec = IPSL_60_2050_max[, 3] - IPSL_60_2050_min[, 3],
                                                          range_IPSL_60_2050_mean_temp_coldest_quarter = IPSL_60_2050_max[, 4] - IPSL_60_2050_min[, 4],
                                                          range_IPSL_60_2050_mean_temp_warmest_quarter = IPSL_60_2050_max[, 5] - IPSL_60_2050_min[, 5],
                                                          range_IPSL_60_2050_prec_driest_quarter = IPSL_60_2050_max[, 6] - IPSL_60_2050_min[, 6], 
                                                          range_IPSL_60_2050_prec_seasonality = IPSL_60_2050_max[, 7] - IPSL_60_2050_min[, 7],
                                                          range_IPSL_60_2050_prec_wettest_quarter = IPSL_60_2050_max[, 8] - IPSL_60_2050_min[, 8],
                                                          #temp_range = IPSL_60_2050_max[, 9] - IPSL_60_2050_min[, 9],
                                                          range_IPSL_60_2050_temp_seasonality = IPSL_60_2050_max[, 10] - IPSL_60_2050_min[, 10])
IPSL_60_2050_range <- IPSL_60_2050_range[-c(2:10)]
write.csv(IPSL_60_2050_range, "/mnt/brunner/brunner/SDM_data/tables/IPSL_60_2050_range.csv", row.names = F)

#=============================== BIOCLIM MIROC rcp6 2041-2060
MIROC_60_2050_list <- list.files(path = "/mnt/brunner/brunner/SDM_data/Bioclim_2050/New_NA/MIROC_global_60/", 
                                 pattern = "*.tif$", 
                                 full.names = TRUE)
MIROC_60_2050_global<- stack(MIROC_60_2050_list)
MIROC_60_2050 <- crop(MIROC_60_2050_global, extent(basin))
MIROC_60_2050 <- mask(MIROC_60_2050, basin)

MIROC_60_2050_mean <- zonal(MIROC_60_2050, basin, fun = mean, digits = 0, na.rm = TRUE)
MIROC_60_2050_mean <-as.data.frame(MIROC_60_2050_mean)
colnames(MIROC_60_2050_mean)[1] <- "basin_id" #rename zone to basin_id to fit ensemble model
colnames(MIROC_60_2050_mean)[2:10] <- paste("mean", colnames(MIROC_60_2050_mean[,c(2:10)]), sep = "_")
write.csv(MIROC_60_2050_mean, "/mnt/brunner/brunner/SDM_data/tables/MIROC_60_2050_mean.csv", row.names = F)

MIROC_60_2050_min <- zonal(MIROC_60_2050, basin, fun = min, digits = 0, na.rm = TRUE)
MIROC_60_2050_min <-as.data.frame(MIROC_60_2050_min)
colnames(MIROC_60_2050_min)[1] <- "basin_id"
colnames(MIROC_60_2050_min)[2:10] <- paste("min", colnames(MIROC_60_2050_min[,c(2:10)]), sep = "_")

MIROC_60_2050_max <- zonal(MIROC_60_2050, basin, fun = max, digits = 0, na.rm = TRUE)
MIROC_60_2050_max  <-as.data.frame(MIROC_60_2050_max)
colnames(MIROC_60_2050_max)[1] <- "basin_id"
colnames(MIROC_60_2050_max)[2:10] <- paste("max", colnames(MIROC_60_2050_max[,c(2:10)]), sep = "_")

MIROC_60_2050_range <- MIROC_60_2050_max %>%  dplyr::mutate(range_MIROC_60_2050_ann_mean_temp = MIROC_60_2050_max[, 2] - MIROC_60_2050_min[, 2], 
                                                          range_MIROC_60_2050_ann_prec = MIROC_60_2050_max[, 3] - MIROC_60_2050_min[, 3],
                                                          range_MIROC_60_2050_mean_temp_coldest_quarter = MIROC_60_2050_max[, 4] - MIROC_60_2050_min[, 4],
                                                          range_MIROC_60_2050_mean_temp_warmest_quarter = MIROC_60_2050_max[, 5] - MIROC_60_2050_min[, 5],
                                                          range_MIROC_60_2050_prec_driest_quarter = MIROC_60_2050_max[, 6] - MIROC_60_2050_min[, 6], 
                                                          range_MIROC_60_2050_prec_seasonality = MIROC_60_2050_max[, 7] - MIROC_60_2050_min[, 7],
                                                          range_MIROC_60_2050_prec_wettest_quarter = MIROC_60_2050_max[, 8] - MIROC_60_2050_min[, 8],
                                                          #temp_range = MIROC_60_2050_max[, 9] - MIROC_60_2050_min[, 9],
                                                          range_MIROC_60_2050_temp_seasonality = MIROC_60_2050_max[, 10] - MIROC_60_2050_min[, 10])
MIROC_60_2050_range <- MIROC_60_2050_range[-c(2:10)]
write.csv(MIROC_60_2050_range, "/mnt/brunner/brunner/SDM_data/tables/MIROC_60_2050_range.csv", row.names = F)

#RCP 8.5
#=============================== BIOCLIM CCSM rcp8.5 2041-2060
CCSM_85_2050_list <- list.files(path = "/mnt/brunner/brunner/SDM_data/Bioclim_2050/New_NA/CCSM_global_85/", 
                                pattern = "*.tif$", 
                                full.names = TRUE)
CCSM_85_2050_global<- stack(CCSM_85_2050_list)
CCSM_85_2050 <- crop(CCSM_85_2050_global, extent(basin))
CCSM_85_2050 <- mask(CCSM_85_2050, basin)

CCSM_85_2050_mean <- zonal(CCSM_85_2050, basin, fun = mean, digits = 0, na.rm = TRUE)
CCSM_85_2050_mean <-as.data.frame(CCSM_85_2050_mean)
colnames(CCSM_85_2050_mean)[1] <- "basin_id" #rename zone to basin_id to fit ensemble model
colnames(CCSM_85_2050_mean)[2:10] <- paste("mean", colnames(CCSM_85_2050_mean[,c(2:10)]), sep = "_")
write.csv(CCSM_85_2050_mean, "/mnt/brunner/brunner/SDM_data/tables/CCSM_85_2050_mean.csv", row.names = F)

CCSM_85_2050_min <- zonal(CCSM_85_2050, basin, fun = min, digits = 0, na.rm = TRUE)
CCSM_85_2050_min <-as.data.frame(CCSM_85_2050_min)
colnames(CCSM_85_2050_min)[1] <- "basin_id"
colnames(CCSM_85_2050_min)[2:10] <- paste("min", colnames(CCSM_85_2050_min[,c(2:10)]), sep = "_")

CCSM_85_2050_max <- zonal(CCSM_85_2050, basin, fun = max, digits = 0, na.rm = TRUE)
CCSM_85_2050_max  <-as.data.frame(CCSM_85_2050_max)
colnames(CCSM_85_2050_max)[1] <- "basin_id"
colnames(CCSM_85_2050_max)[2:10] <- paste("max", colnames(CCSM_85_2050_max[,c(2:10)]), sep = "_")

CCSM_85_2050_range <- CCSM_85_2050_max %>%  dplyr::mutate(range_CCSM_85_2050_ann_mean_temp = CCSM_85_2050_max[, 2] - CCSM_85_2050_min[, 2], 
                                                          range_CCSM_85_2050_ann_prec = CCSM_85_2050_max[, 3] - CCSM_85_2050_min[, 3],
                                                          range_CCSM_85_2050_mean_temp_coldest_quarter = CCSM_85_2050_max[, 4] - CCSM_85_2050_min[, 4],
                                                          range_CCSM_85_2050_mean_temp_warmest_quarter = CCSM_85_2050_max[, 5] - CCSM_85_2050_min[, 5],
                                                          range_CCSM_85_2050_prec_driest_quarter = CCSM_85_2050_max[, 6] - CCSM_85_2050_min[, 6], 
                                                          range_CCSM_85_2050_prec_seasonality = CCSM_85_2050_max[, 7] - CCSM_85_2050_min[, 7],
                                                          range_CCSM_85_2050_prec_wettest_quarter = CCSM_85_2050_max[, 8] - CCSM_85_2050_min[, 8],
                                                          #temp_range = CCSM_85_2050_max[, 9] - CCSM_85_2050_min[, 9],
                                                          range_CCSM_85_2050_temp_seasonality = CCSM_85_2050_max[, 10] - CCSM_85_2050_min[, 10])
CCSM_85_2050_range <- CCSM_85_2050_range[-c(2:10)]
write.csv(CCSM_85_2050_range, "/mnt/brunner/brunner/SDM_data/tables/CCSM_85_2050_range.csv", row.names = F)





#=============================== BIOCLIM IPSL rcp8.5 2041-2060
IPSL_85_2050_list <- list.files(path = "/mnt/brunner/brunner/SDM_data/Bioclim_2050/New_NA/IPSL_global_85/", 
                                pattern = "*.tif$", 
                                full.names = TRUE)
IPSL_85_2050_global<- stack(IPSL_85_2050_list)
IPSL_85_2050 <- crop(IPSL_85_2050_global, extent(basin))
IPSL_85_2050 <- mask(IPSL_85_2050, basin)

IPSL_85_2050_mean <- zonal(IPSL_85_2050, basin, fun = mean, digits = 0, na.rm = TRUE)
IPSL_85_2050_mean <-as.data.frame(IPSL_85_2050_mean)
colnames(IPSL_85_2050_mean)[1] <- "basin_id" #rename zone to basin_id to fit ensemble model
colnames(IPSL_85_2050_mean)[2:10] <- paste("mean", colnames(IPSL_85_2050_mean[,c(2:10)]), sep = "_")
write.csv(IPSL_85_2050_mean, "/mnt/brunner/brunner/SDM_data/tables/IPSL_85_2050_mean.csv", row.names = F)

IPSL_85_2050_min <- zonal(IPSL_85_2050, basin, fun = min, digits = 0, na.rm = TRUE)
IPSL_85_2050_min <-as.data.frame(IPSL_85_2050_min)
colnames(IPSL_85_2050_min)[1] <- "basin_id"
colnames(IPSL_85_2050_min)[2:10] <- paste("min", colnames(IPSL_85_2050_min[,c(2:10)]), sep = "_")

IPSL_85_2050_max <- zonal(IPSL_85_2050, basin, fun = max, digits = 0, na.rm = TRUE)
IPSL_85_2050_max  <-as.data.frame(IPSL_85_2050_max)
colnames(IPSL_85_2050_max)[1] <- "basin_id"
colnames(IPSL_85_2050_max)[2:10] <- paste("max", colnames(IPSL_85_2050_max[,c(2:10)]), sep = "_")

IPSL_85_2050_range <- IPSL_85_2050_max %>%  dplyr::mutate(range_IPSL_85_2050_ann_mean_temp = IPSL_85_2050_max[, 2] - IPSL_85_2050_min[, 2], 
                                                          range_IPSL_85_2050_ann_prec = IPSL_85_2050_max[, 3] - IPSL_85_2050_min[, 3],
                                                          range_IPSL_85_2050_mean_temp_coldest_quarter = IPSL_85_2050_max[, 4] - IPSL_85_2050_min[, 4],
                                                          range_IPSL_85_2050_mean_temp_warmest_quarter = IPSL_85_2050_max[, 5] - IPSL_85_2050_min[, 5],
                                                          range_IPSL_85_2050_prec_driest_quarter = IPSL_85_2050_max[, 6] - IPSL_85_2050_min[, 6], 
                                                          range_IPSL_85_2050_prec_seasonality = IPSL_85_2050_max[, 7] - IPSL_85_2050_min[, 7],
                                                          range_IPSL_85_2050_prec_wettest_quarter = IPSL_85_2050_max[, 8] - IPSL_85_2050_min[, 8],
                                                          #temp_range = IPSL_85_2050_max[, 9] - IPSL_85_2050_min[, 9],
                                                          range_IPSL_85_2050_temp_seasonality = IPSL_85_2050_max[, 10] - IPSL_85_2050_min[, 10])
IPSL_85_2050_range <- IPSL_85_2050_range[-c(2:10)]
write.csv(IPSL_85_2050_range, "/mnt/brunner/brunner/SDM_data/tables/IPSL_85_2050_range.csv", row.names = F)


#=============================== BIOCLIM MIROC rcp8.5 2041-2060
MIROC_85_2050_list <- list.files(path = "/mnt/brunner/brunner/SDM_data/Bioclim_2050/New_NA/MIROC_global_85/", 
                                 pattern = "*.tif$", 
                                 full.names = TRUE)
MIROC_85_2050_global<- stack(MIROC_85_2050_list)
MIROC_85_2050 <- crop(MIROC_85_2050_global, extent(basin))
MIROC_85_2050 <- mask(MIROC_85_2050, basin)

MIROC_85_2050_mean <- zonal(MIROC_85_2050, basin, fun = mean, digits = 0, na.rm = TRUE)
MIROC_85_2050_mean <-as.data.frame(MIROC_85_2050_mean)
colnames(MIROC_85_2050_mean)[1] <- "basin_id" #rename zone to basin_id to fit ensemble model
colnames(MIROC_85_2050_mean)[2:10] <- paste("mean", colnames(MIROC_85_2050_mean[,c(2:10)]), sep = "_")
write.csv(MIROC_85_2050_mean, "/mnt/brunner/brunner/SDM_data/tables/MIROC_85_2050_mean.csv", row.names = F)
  
  MIROC_85_2050_min <- zonal(MIROC_85_2050, basin, fun = min, digits = 0, na.rm = TRUE)
  MIROC_85_2050_min <-as.data.frame(MIROC_85_2050_min)
  colnames(MIROC_85_2050_min)[1] <- "basin_id"
  colnames(MIROC_85_2050_min)[2:10] <- paste("min", colnames(MIROC_85_2050_min[,c(2:10)]), sep = "_")
  
  MIROC_85_2050_max <- zonal(MIROC_85_2050, basin, fun = max, digits = 0, na.rm = TRUE)
  MIROC_85_2050_max  <-as.data.frame(MIROC_85_2050_max)
  colnames(MIROC_85_2050_max)[1] <- "basin_id"
  colnames(MIROC_85_2050_max)[2:10] <- paste("max", colnames(MIROC_85_2050_max[,c(2:10)]), sep = "_")
  
  MIROC_85_2050_range <- MIROC_85_2050_max %>%  dplyr::mutate(range_MIROC_85_2050_ann_mean_temp = MIROC_85_2050_max[, 2] - MIROC_85_2050_min[, 2], 
                                                              range_MIROC_85_2050_ann_prec = MIROC_85_2050_max[, 3] - MIROC_85_2050_min[, 3],
                                                              range_MIROC_85_2050_mean_temp_coldest_quarter = MIROC_85_2050_max[, 4] - MIROC_85_2050_min[, 4],
                                                              range_MIROC_85_2050_mean_temp_warmest_quarter = MIROC_85_2050_max[, 5] - MIROC_85_2050_min[, 5],
                                                              range_MIROC_85_2050_prec_driest_quarter = MIROC_85_2050_max[, 6] - MIROC_85_2050_min[, 6], 
                                                              range_MIROC_85_2050_prec_seasonality = MIROC_85_2050_max[, 7] - MIROC_85_2050_min[, 7],
                                                              range_MIROC_85_2050_prec_wettest_quarter = MIROC_85_2050_max[, 8] - MIROC_85_2050_min[, 8],
                                                              #temp_range = MIROC_85_2050_max[, 9] - MIROC_85_2050_min[, 9],
                                                              range_MIROC_85_2050_temp_seasonality = MIROC_85_2050_max[, 10] - MIROC_85_2050_min[, 10])
  MIROC_85_2050_range <- MIROC_85_2050_range[-c(2:10)]
  write.csv(MIROC_85_2050_range, "/mnt/brunner/brunner/SDM_data/tables/MIROC_85_2050_range.csv", row.names = F)

  #=============================== LC percent per basin 2050
  #=========SSP1
  #LANDCOVER SSP1 2050
  LC_SSP1_2050_LAC <- raster("/mnt/brunner/brunner/dinamica/final_model_output/SSP1/SSP1_05.tif")
  LC_SSP1_2050 <- crop(LC_SSP1_2050_LAC, extent(study_area))
  LC_SSP1_2050 <- mask(LC_SSP1_2050, study_area)
  extent(LC_SSP1_2050) <- alignExtent(LC_SSP1_2050, basin)
  #NAvalue(LC_SSP1_2050)<-(-9999)
  #x11(); plot(LC_SSP1_2050)

  
  ##pixelcount for the SSP1 LC 2050 map categories per basin
  SSP1_2050_zonal <-myZonal(LC_SSP1_2050, basin)
  results_SSP1_2050 <- data.table::dcast(SSP1_2050_zonal, z ~ vals)
  colnames(results_SSP1_2050)[1] <- "basin_id"
  results_SSP1_2050[is.na(results_SSP1_2050)] <- 0
  
  #bind basin count and LC 10 count
  #all_results <- cbind(results_basin, results_SSP1_2050 )
  SSP1_2050_LC <- merge(results_basin, results_SSP1_2050, by = "basin_id", all= T )
  SSP1_2050_LC <- SSP1_2050_LC[-1,]
  
  
  LC_percentage_SSP1_2050 <- (SSP1_2050_LC[,3:8]/SSP1_2050_LC[,2])
  LC_percentage_SSP1_2050 <- round(LC_percentage_SSP1_2050, digits = 5)
  LC_percentage_SSP1_2050 <- cbind(SSP1_2050_LC[1], LC_percentage_SSP1_2050)

  
  LC_percentage_SSP1_2050 <- as.data.frame(LC_percentage_SSP1_2050)
  setnames(LC_percentage_SSP1_2050, old = c('1', '2', '3', '4', '5', '6'), new = c('crops', 'forest','pasture','urban', 'other', 'water'))
  write.csv(LC_percentage_SSP1_2050,"/mnt/brunner/brunner/SDM_data/tables/LC_percentage_SSP1_2050.csv", row.names = F )
  
  #=========SSP2
  #LANDCOVER SSP2 2050
  LC_SSP2_2050_LAC <- raster("/mnt/brunner/brunner/dinamica/final_model_output/SSP2/SSP2_05.tif")
  LC_SSP2_2050 <- crop(LC_SSP2_2050_LAC, extent(study_area))
  LC_SSP2_2050 <- mask(LC_SSP2_2050, study_area)
  extent(LC_SSP2_2050) <- alignExtent(LC_SSP2_2050, basin)
  #NAvalue(LC_SSP2_2050)<-(-9999)
 # x11(); plot(LC_SSP2_2050)
  
  
  ##pixelcount for the SSP2 LC 2050 map categories per basin
  SSP2_2050_zonal <-myZonal(LC_SSP2_2050, basin)
  results_SSP2_2050 <- data.table::dcast(SSP2_2050_zonal, z ~ vals)
  colnames(results_SSP2_2050)[1] <- "basin_id"
  results_SSP2_2050[is.na(results_SSP2_2050)] <- 0
  
  #bind basin count and LC 10 count
  #all_results <- cbind(results_basin, results_SSP2_2050 )
  SSP2_2050_LC <- merge(results_basin, results_SSP2_2050, by = "basin_id", all= T )
  SSP2_2050_LC <- SSP2_2050_LC[-1,]
  
  
  LC_percentage_SSP2_2050 <- (SSP2_2050_LC[,3:8]/SSP2_2050_LC[,2])
  LC_percentage_SSP2_2050 <- round(LC_percentage_SSP2_2050, digits = 5)
  LC_percentage_SSP2_2050 <- cbind(SSP2_2050_LC[1], LC_percentage_SSP2_2050)
  
  
  LC_percentage_SSP2_2050 <- as.data.frame(LC_percentage_SSP2_2050)
  setnames(LC_percentage_SSP2_2050, old = c('1', '2', '3', '4', '5', '6'), new = c('crops', 'forest','pasture','urban', 'other', 'water'))
  write.csv(LC_percentage_SSP2_2050,"/mnt/brunner/brunner/SDM_data/tables/LC_percentage_SSP2_2050.csv", row.names = F )
  
  
  #=========SSP3
  #LANDCOVER SSP3 2050
  LC_SSP3_2050_LAC <- raster("/mnt/brunner/brunner/dinamica/final_model_output/SSP3/SSP3_05.tif")
  LC_SSP3_2050 <- crop(LC_SSP3_2050_LAC, extent(study_area))
  LC_SSP3_2050 <- mask(LC_SSP3_2050, study_area)
  extent(LC_SSP3_2050) <- alignExtent(LC_SSP3_2050, basin)
  #NAvalue(LC_SSP3_2050)<-(-9999)
  #x11(); plot(LC_SSP3_2050)
  
  
  ##pixelcount for the SSP3 LC 2050 map categories per basin
  SSP3_2050_zonal <-myZonal(LC_SSP3_2050, basin)
  results_SSP3_2050 <- data.table::dcast(SSP3_2050_zonal, z ~ vals)
  colnames(results_SSP3_2050)[1] <- "basin_id"
  results_SSP3_2050[is.na(results_SSP3_2050)] <- 0
  
  #bind basin count and LC 10 count
  #all_results <- cbind(results_basin, results_SSP3_2050 )
  SSP3_2050_LC <- merge(results_basin, results_SSP3_2050, by = "basin_id", all= T )
  SSP3_2050_LC <- SSP3_2050_LC[-1,]
  
  
  LC_percentage_SSP3_2050 <- (SSP3_2050_LC[,3:8]/SSP3_2050_LC[,2])
  LC_percentage_SSP3_2050 <- round(LC_percentage_SSP3_2050, digits = 3)
  LC_percentage_SSP3_2050 <- cbind(SSP3_2050_LC[1], LC_percentage_SSP3_2050)
  
  
  LC_percentage_SSP3_2050 <- as.data.frame(LC_percentage_SSP3_2050)
  setnames(LC_percentage_SSP3_2050, old = c('1', '2', '3', '4', '5', '6'), new = c('crops', 'forest','pasture','urban', 'other', 'water'))
  write.csv(LC_percentage_SSP3_2050,"/mnt/brunner/brunner/SDM_data/tables/LC_percentage_SSP3_2050.csv", row.names = F )
  
  #=========SSP4
  #LANDCOVER SSP4 2050
  LC_SSP4_2050_LAC <- raster("/mnt/brunner/brunner/dinamica/final_model_output/SSP4/SSP4_05.tif")
  LC_SSP4_2050 <- crop(LC_SSP4_2050_LAC, extent(study_area))
  LC_SSP4_2050 <- mask(LC_SSP4_2050, study_area)
  extent(LC_SSP4_2050) <- alignExtent(LC_SSP4_2050, basin)
  #NAvalue(LC_SSP4_2050)<-(-9999)
 # x11(); plot(LC_SSP4_2050)
  
  
  ##pixelcount for the SSP4 LC 2050 map categories per basin
  SSP4_2050_zonal <-myZonal(LC_SSP4_2050, basin)
  results_SSP4_2050 <- data.table::dcast(SSP4_2050_zonal, z ~ vals)
  colnames(results_SSP4_2050)[1] <- "basin_id"
  results_SSP4_2050[is.na(results_SSP4_2050)] <- 0
  
  #bind basin count and LC 10 count
  #all_results <- cbind(results_basin, results_SSP4_2050 )
  SSP4_2050_LC <- merge(results_basin, results_SSP4_2050, by = "basin_id", all= T )
  SSP4_2050_LC <- SSP4_2050_LC[-1,]
  
  
  LC_percentage_SSP4_2050 <- (SSP4_2050_LC[,3:8]/SSP4_2050_LC[,2])
  LC_percentage_SSP4_2050 <- round(LC_percentage_SSP4_2050, digits = 5)
  LC_percentage_SSP4_2050 <- cbind(SSP4_2050_LC[1], LC_percentage_SSP4_2050)
  
  
  LC_percentage_SSP4_2050 <- as.data.frame(LC_percentage_SSP4_2050)
  setnames(LC_percentage_SSP4_2050, old = c('1', '2', '3', '4', '5', '6'), new = c('crops', 'forest','pasture','urban', 'other', 'water'))
  write.csv(LC_percentage_SSP4_2050,"/mnt/brunner/brunner/SDM_data/tables/LC_percentage_SSP4_2050.csv", row.names = F )
  
  #=========SSP5
  #LANDCOVER SSP5 2050
  LC_SSP5_2050_LAC <- raster("/mnt/brunner/brunner/dinamica/final_model_output/SSP5/SSP5_05.tif")
  LC_SSP5_2050 <- crop(LC_SSP5_2050_LAC, extent(study_area))
  LC_SSP5_2050 <- mask(LC_SSP5_2050, study_area)
  extent(LC_SSP5_2050) <- alignExtent(LC_SSP5_2050, basin)
  #NAvalue(LC_SSP5_2050)<-(-9999)
 # x11(); plot(LC_SSP5_2050)
  
  
  ##pixelcount for the SSP5 LC 2050 map categories per basin
  SSP5_2050_zonal <-myZonal(LC_SSP5_2050, basin)
  results_SSP5_2050 <- data.table::dcast(SSP5_2050_zonal, z ~ vals)
  colnames(results_SSP5_2050)[1] <- "basin_id"
  results_SSP5_2050[is.na(results_SSP5_2050)] <- 0
  
  #bind basin count and LC 10 count
  #all_results <- cbind(results_basin, results_SSP5_2050 )
  SSP5_2050_LC <- merge(results_basin, results_SSP5_2050, by = "basin_id", all= T )
  SSP5_2050_LC <- SSP5_2050_LC[-1,]
  
  
  LC_percentage_SSP5_2050 <- (SSP5_2050_LC[,3:8]/SSP5_2050_LC[,2])
  LC_percentage_SSP5_2050 <- round(LC_percentage_SSP5_2050, digits = 5)
  LC_percentage_SSP5_2050 <- cbind(SSP5_2050_LC[1], LC_percentage_SSP5_2050)
  
  
  LC_percentage_SSP5_2050 <- as.data.frame(LC_percentage_SSP5_2050)
  setnames(LC_percentage_SSP5_2050, old = c('1', '2', '3', '4', '5', '6'), new = c('crops', 'forest','pasture','urban', 'other', 'water'))
  write.csv(LC_percentage_SSP5_2050,"/mnt/brunner/brunner/SDM_data/tables/LC_percentage_SSP5_2050.csv", row.names = F )
  
  
  
  
  #========================3: climate, topographical and Landcover 2070 ===================================#
  
  #RCP 6
  
  #=============================== BIOCLIM CCSM rcp6 2061-2080
  CCSM_60_2070_list <- list.files(path = "/mnt/brunner/brunner/SDM_data/Bioclim_2070/New_NA_2070/CCSM_global_60/", 
                                  pattern = "*.tif$", 
                                  full.names = TRUE)
  CCSM_60_2070_global<- stack(CCSM_60_2070_list)
  CCSM_60_2070 <- crop(CCSM_60_2070_global, extent(basin))
  CCSM_60_2070 <- mask(CCSM_60_2070, basin)
  #x11(); plot(CCSM_60_2070)
  
  CCSM_60_2070_mean <- zonal(CCSM_60_2070, basin, fun = mean, digits = 0, na.rm = TRUE)
  CCSM_60_2070_mean <-as.data.frame(CCSM_60_2070_mean)
  colnames(CCSM_60_2070_mean)[1] <- "basin_id" #rename zone to basin_id to fit ensemble model
  colnames(CCSM_60_2070_mean)[2:10] <- paste("mean", colnames(CCSM_60_2070_mean[,c(2:10)]), sep = "_")
  write.csv(CCSM_60_2070_mean, "/mnt/brunner/brunner/SDM_data/tables/CCSM_60_2070_mean.csv", row.names = F)
  
  CCSM_60_2070_min <- zonal(CCSM_60_2070, basin, fun = min, digits = 0, na.rm = TRUE)
  CCSM_60_2070_min <-as.data.frame(CCSM_60_2070_min)
  colnames(CCSM_60_2070_min)[1] <- "basin_id"
  colnames(CCSM_60_2070_min)[2:10] <- paste("min", colnames(CCSM_60_2070_min[,c(2:10)]), sep = "_")
  
  CCSM_60_2070_max <- zonal(CCSM_60_2070, basin, fun = max, digits = 0, na.rm = TRUE)
  CCSM_60_2070_max  <-as.data.frame(CCSM_60_2070_max)
  colnames(CCSM_60_2070_max)[1] <- "basin_id"
  colnames(CCSM_60_2070_max)[2:10] <- paste("max", colnames(CCSM_60_2070_max[,c(2:10)]), sep = "_")
  
  CCSM_60_2070_range <- CCSM_60_2070_max %>%  dplyr::mutate(range_CCSM_60_2070_ann_mean_temp = CCSM_60_2070_max[, 2] - CCSM_60_2070_min[, 2], 
                                                            range_CCSM_60_2070_ann_prec = CCSM_60_2070_max[, 3] - CCSM_60_2070_min[, 3],
                                                            range_CCSM_60_2070_mean_temp_coldest_quarter = CCSM_60_2070_max[, 4] - CCSM_60_2070_min[, 4],
                                                            range_CCSM_60_2070_mean_temp_warmest_quarter = CCSM_60_2070_max[, 5] - CCSM_60_2070_min[, 5],
                                                            range_CCSM_60_2070_prec_driest_quarter = CCSM_60_2070_max[, 6] - CCSM_60_2070_min[, 6], 
                                                            range_CCSM_60_2070_prec_seasonality = CCSM_60_2070_max[, 7] - CCSM_60_2070_min[, 7],
                                                            range_CCSM_60_2070_prec_wettest_quarter = CCSM_60_2070_max[, 8] - CCSM_60_2070_min[, 8],
                                                            #temp_range = CCSM_60_2070_max[, 9] - CCSM_60_2070_min[, 9],
                                                            range_CCSM_60_2070_temp_seasonality = CCSM_60_2070_max[, 10] - CCSM_60_2070_min[, 10])
  CCSM_60_2070_range <- CCSM_60_2070_range[-c(2:10)]
  write.csv(CCSM_60_2070_range, "/mnt/brunner/brunner/SDM_data/tables/CCSM_60_2070_range.csv", row.names = F)
  
  #=============================== BIOCLIM IPSL rcp6 2061-2080
  IPSL_60_2070_list <- list.files(path = "/mnt/brunner/brunner/SDM_data/Bioclim_2070/New_NA_2070/IPSL_global_60/", 
                                  pattern = "*.tif$", 
                                  full.names = TRUE)
  IPSL_60_2070_global<- stack(IPSL_60_2070_list)
  IPSL_60_2070 <- crop(IPSL_60_2070_global, extent(basin))
  IPSL_60_2070 <- mask(IPSL_60_2070, basin)
  
  IPSL_60_2070_mean <- zonal(IPSL_60_2070, basin, fun = mean, digits = 0, na.rm = TRUE)
  IPSL_60_2070_mean <-as.data.frame(IPSL_60_2070_mean)
  colnames(IPSL_60_2070_mean)[1] <- "basin_id" #rename zone to basin_id to fit ensemble model
  colnames(IPSL_60_2070_mean)[2:10] <- paste("mean", colnames(IPSL_60_2070_mean[,c(2:10)]), sep = "_")
  write.csv(IPSL_60_2070_mean, "/mnt/brunner/brunner/SDM_data/tables/IPSL_60_2070_mean.csv", row.names = F)
  
  IPSL_60_2070_min <- zonal(IPSL_60_2070, basin, fun = min, digits = 0, na.rm = TRUE)
  IPSL_60_2070_min <-as.data.frame(IPSL_60_2070_min)
  colnames(IPSL_60_2070_min)[1] <- "basin_id"
  colnames(IPSL_60_2070_min)[2:10] <- paste("min", colnames(IPSL_60_2070_min[,c(2:10)]), sep = "_")
  
  IPSL_60_2070_max <- zonal(IPSL_60_2070, basin, fun = max, digits = 0, na.rm = TRUE)
  IPSL_60_2070_max  <-as.data.frame(IPSL_60_2070_max)
  colnames(IPSL_60_2070_max)[1] <- "basin_id"
  
  
  
  
  colnames(IPSL_60_2070_max)[2:10] <- paste("max", colnames(IPSL_60_2070_max[,c(2:10)]), sep = "_")
  
  IPSL_60_2070_range <- IPSL_60_2070_max %>%  dplyr::mutate(range_IPSL_60_2070_ann_mean_temp = IPSL_60_2070_max[, 2] - IPSL_60_2070_min[, 2], 
                                                            range_IPSL_60_2070_ann_prec = IPSL_60_2070_max[, 3] - IPSL_60_2070_min[, 3],
                                                            range_IPSL_60_2070_mean_temp_coldest_quarter = IPSL_60_2070_max[, 4] - IPSL_60_2070_min[, 4],
                                                            range_IPSL_60_2070_mean_temp_warmest_quarter = IPSL_60_2070_max[, 5] - IPSL_60_2070_min[, 5],
                                                            range_IPSL_60_2070_prec_driest_quarter = IPSL_60_2070_max[, 6] - IPSL_60_2070_min[, 6], 
                                                            range_IPSL_60_2070_prec_seasonality = IPSL_60_2070_max[, 7] - IPSL_60_2070_min[, 7],
                                                            range_IPSL_60_2070_prec_wettest_quarter = IPSL_60_2070_max[, 8] - IPSL_60_2070_min[, 8],
                                                            #temp_range = IPSL_60_2070_max[, 9] - IPSL_60_2070_min[, 9],
                                                            range_IPSL_60_2070_temp_seasonality = IPSL_60_2070_max[, 10] - IPSL_60_2070_min[, 10])
  IPSL_60_2070_range <- IPSL_60_2070_range[-c(2:10)]
  write.csv(IPSL_60_2070_range, "/mnt/brunner/brunner/SDM_data/tables/IPSL_60_2070_range.csv", row.names = F)
  
  #=============================== BIOCLIM MIROC rcp6 2061-2080
  MIROC_60_2070_list <- list.files(path = "/mnt/brunner/brunner/SDM_data/Bioclim_2070/New_NA_2070/MIROC_global_60/", 
                                   pattern = "*.tif$", 
                                   full.names = TRUE)
  MIROC_60_2070_global<- stack(MIROC_60_2070_list)
  MIROC_60_2070 <- crop(MIROC_60_2070_global, extent(basin))
  MIROC_60_2070 <- mask(MIROC_60_2070, basin)
  
  MIROC_60_2070_mean <- zonal(MIROC_60_2070, basin, fun = mean, digits = 0, na.rm = TRUE)
  MIROC_60_2070_mean <-as.data.frame(MIROC_60_2070_mean)
  colnames(MIROC_60_2070_mean)[1] <- "basin_id" #rename zone to basin_id to fit ensemble model
  colnames(MIROC_60_2070_mean)[2:10] <- paste("mean", colnames(MIROC_60_2070_mean[,c(2:10)]), sep = "_")
  write.csv(MIROC_60_2070_mean, "/mnt/brunner/brunner/SDM_data/tables/MIROC_60_2070_mean.csv", row.names = F)
  
  MIROC_60_2070_min <- zonal(MIROC_60_2070, basin, fun = min, digits = 0, na.rm = TRUE)
  MIROC_60_2070_min <-as.data.frame(MIROC_60_2070_min)
  colnames(MIROC_60_2070_min)[1] <- "basin_id"
  colnames(MIROC_60_2070_min)[2:10] <- paste("min", colnames(MIROC_60_2070_min[,c(2:10)]), sep = "_")
  
  MIROC_60_2070_max <- zonal(MIROC_60_2070, basin, fun = max, digits = 0, na.rm = TRUE)
  MIROC_60_2070_max  <-as.data.frame(MIROC_60_2070_max)
  colnames(MIROC_60_2070_max)[1] <- "basin_id"
  colnames(MIROC_60_2070_max)[2:10] <- paste("max", colnames(MIROC_60_2070_max[,c(2:10)]), sep = "_")
  
  MIROC_60_2070_range <- MIROC_60_2070_max %>%  dplyr::mutate(range_MIROC_60_2070_ann_mean_temp = MIROC_60_2070_max[, 2] - MIROC_60_2070_min[, 2], 
                                                              range_MIROC_60_2070_ann_prec = MIROC_60_2070_max[, 3] - MIROC_60_2070_min[, 3],
                                                              range_MIROC_60_2070_mean_temp_coldest_quarter = MIROC_60_2070_max[, 4] - MIROC_60_2070_min[, 4],
                                                              range_MIROC_60_2070_mean_temp_warmest_quarter = MIROC_60_2070_max[, 5] - MIROC_60_2070_min[, 5],
                                                              range_MIROC_60_2070_prec_driest_quarter = MIROC_60_2070_max[, 6] - MIROC_60_2070_min[, 6], 
                                                              range_MIROC_60_2070_prec_seasonality = MIROC_60_2070_max[, 7] - MIROC_60_2070_min[, 7],
                                                              range_MIROC_60_2070_prec_wettest_quarter = MIROC_60_2070_max[, 8] - MIROC_60_2070_min[, 8],
                                                              #temp_range = MIROC_60_2070_max[, 9] - MIROC_60_2070_min[, 9],
                                                              range_MIROC_60_2070_temp_seasonality = MIROC_60_2070_max[, 10] - MIROC_60_2070_min[, 10])
  MIROC_60_2070_range <- MIROC_60_2070_range[-c(2:10)]
  write.csv(MIROC_60_2070_range, "/mnt/brunner/brunner/SDM_data/tables/MIROC_60_2070_range.csv", row.names = F)
  
  #RCP 8.5
  #=============================== BIOCLIM CCSM rcp8.5 2061-2080
  CCSM_85_2070_list <- list.files(path = "/mnt/brunner/brunner/SDM_data/Bioclim_2070/New_NA_2070/CCSM_global_85/", 
                                  pattern = "*.tif$", 
                                  full.names = TRUE)
  CCSM_85_2070_global<- stack(CCSM_85_2070_list)
  CCSM_85_2070 <- crop(CCSM_85_2070_global, extent(basin))
  CCSM_85_2070 <- mask(CCSM_85_2070, basin)
  
  CCSM_85_2070_mean <- zonal(CCSM_85_2070, basin, fun = mean, digits = 0, na.rm = TRUE)
  CCSM_85_2070_mean <-as.data.frame(CCSM_85_2070_mean)
  colnames(CCSM_85_2070_mean)[1] <- "basin_id" #rename zone to basin_id to fit ensemble model
  colnames(CCSM_85_2070_mean)[2:10] <- paste("mean", colnames(CCSM_85_2070_mean[,c(2:10)]), sep = "_")
  write.csv(CCSM_85_2070_mean, "/mnt/brunner/brunner/SDM_data/tables/CCSM_85_2070_mean.csv", row.names = F)
  
  CCSM_85_2070_min <- zonal(CCSM_85_2070, basin, fun = min, digits = 0, na.rm = TRUE)
  CCSM_85_2070_min <-as.data.frame(CCSM_85_2070_min)
  colnames(CCSM_85_2070_min)[1] <- "basin_id"
  colnames(CCSM_85_2070_min)[2:10] <- paste("min", colnames(CCSM_85_2070_min[,c(2:10)]), sep = "_")
  
  CCSM_85_2070_max <- zonal(CCSM_85_2070, basin, fun = max, digits = 0, na.rm = TRUE)
  CCSM_85_2070_max  <-as.data.frame(CCSM_85_2070_max)
  colnames(CCSM_85_2070_max)[1] <- "basin_id"
  colnames(CCSM_85_2070_max)[2:10] <- paste("max", colnames(CCSM_85_2070_max[,c(2:10)]), sep = "_")
  
  CCSM_85_2070_range <- CCSM_85_2070_max %>%  dplyr::mutate(range_CCSM_85_2070_ann_mean_temp = CCSM_85_2070_max[, 2] - CCSM_85_2070_min[, 2], 
                                                            range_CCSM_85_2070_ann_prec = CCSM_85_2070_max[, 3] - CCSM_85_2070_min[, 3],
                                                            range_CCSM_85_2070_mean_temp_coldest_quarter = CCSM_85_2070_max[, 4] - CCSM_85_2070_min[, 4],
                                                            range_CCSM_85_2070_mean_temp_warmest_quarter = CCSM_85_2070_max[, 5] - CCSM_85_2070_min[, 5],
                                                            range_CCSM_85_2070_prec_driest_quarter = CCSM_85_2070_max[, 6] - CCSM_85_2070_min[, 6], 
                                                            range_CCSM_85_2070_prec_seasonality = CCSM_85_2070_max[, 7] - CCSM_85_2070_min[, 7],
                                                            range_CCSM_85_2070_prec_wettest_quarter = CCSM_85_2070_max[, 8] - CCSM_85_2070_min[, 8],
                                                            #temp_range = CCSM_85_2070_max[, 9] - CCSM_85_2070_min[, 9],
                                                            range_CCSM_85_2070_temp_seasonality = CCSM_85_2070_max[, 10] - CCSM_85_2070_min[, 10])
  CCSM_85_2070_range <- CCSM_85_2070_range[-c(2:10)]
  write.csv(CCSM_85_2070_range, "/mnt/brunner/brunner/SDM_data/tables/CCSM_85_2070_range.csv", row.names = F)
  
  
  
  
  
  #=============================== BIOCLIM IPSL rcp8.5 2061-2080
  IPSL_85_2070_list <- list.files(path = "/mnt/brunner/brunner/SDM_data/Bioclim_2070/New_NA_2070/IPSL_global_85/", 
                                  pattern = "*.tif$", 
                                  full.names = TRUE)
  IPSL_85_2070_global<- stack(IPSL_85_2070_list)
  IPSL_85_2070 <- crop(IPSL_85_2070_global, extent(basin))
  IPSL_85_2070 <- mask(IPSL_85_2070, basin)
  
  IPSL_85_2070_mean <- zonal(IPSL_85_2070, basin, fun = mean, digits = 0, na.rm = TRUE)
  IPSL_85_2070_mean <-as.data.frame(IPSL_85_2070_mean)
  colnames(IPSL_85_2070_mean)[1] <- "basin_id" #rename zone to basin_id to fit ensemble model
  colnames(IPSL_85_2070_mean)[2:10] <- paste("mean", colnames(IPSL_85_2070_mean[,c(2:10)]), sep = "_")
  write.csv(IPSL_85_2070_mean, "/mnt/brunner/brunner/SDM_data/tables/IPSL_85_2070_mean.csv", row.names = F)
  
  IPSL_85_2070_min <- zonal(IPSL_85_2070, basin, fun = min, digits = 0, na.rm = TRUE)
  IPSL_85_2070_min <-as.data.frame(IPSL_85_2070_min)
  colnames(IPSL_85_2070_min)[1] <- "basin_id"
  colnames(IPSL_85_2070_min)[2:10] <- paste("min", colnames(IPSL_85_2070_min[,c(2:10)]), sep = "_")
  
  IPSL_85_2070_max <- zonal(IPSL_85_2070, basin, fun = max, digits = 0, na.rm = TRUE)
  IPSL_85_2070_max  <-as.data.frame(IPSL_85_2070_max)
  colnames(IPSL_85_2070_max)[1] <- "basin_id"
  colnames(IPSL_85_2070_max)[2:10] <- paste("max", colnames(IPSL_85_2070_max[,c(2:10)]), sep = "_")
  
  IPSL_85_2070_range <- IPSL_85_2070_max %>%  dplyr::mutate(range_IPSL_85_2070_ann_mean_temp = IPSL_85_2070_max[, 2] - IPSL_85_2070_min[, 2], 
                                                            range_IPSL_85_2070_ann_prec = IPSL_85_2070_max[, 3] - IPSL_85_2070_min[, 3],
                                                            range_IPSL_85_2070_mean_temp_coldest_quarter = IPSL_85_2070_max[, 4] - IPSL_85_2070_min[, 4],
                                                            range_IPSL_85_2070_mean_temp_warmest_quarter = IPSL_85_2070_max[, 5] - IPSL_85_2070_min[, 5],
                                                            range_IPSL_85_2070_prec_driest_quarter = IPSL_85_2070_max[, 6] - IPSL_85_2070_min[, 6], 
                                                            range_IPSL_85_2070_prec_seasonality = IPSL_85_2070_max[, 7] - IPSL_85_2070_min[, 7],
                                                            range_IPSL_85_2070_prec_wettest_quarter = IPSL_85_2070_max[, 8] - IPSL_85_2070_min[, 8],
                                                            #temp_range = IPSL_85_2070_max[, 9] - IPSL_85_2070_min[, 9],
                                                            range_IPSL_85_2070_temp_seasonality = IPSL_85_2070_max[, 10] - IPSL_85_2070_min[, 10])
  IPSL_85_2070_range <- IPSL_85_2070_range[-c(2:10)]
  write.csv(IPSL_85_2070_range, "/mnt/brunner/brunner/SDM_data/tables/IPSL_85_2070_range.csv", row.names = F)
  
  
  #=============================== BIOCLIM MIROC rcp8.5 2061-2080
  MIROC_85_2070_list <- list.files(path = "/mnt/brunner/brunner/SDM_data/Bioclim_2070/New_NA_2070/MIROC_global_85/", 
                                   pattern = "*.tif$", 
                                   full.names = TRUE)
  MIROC_85_2070_global<- stack(MIROC_85_2070_list)
  MIROC_85_2070 <- crop(MIROC_85_2070_global, extent(basin))
  MIROC_85_2070 <- mask(MIROC_85_2070, basin)
  
  MIROC_85_2070_mean <- zonal(MIROC_85_2070, basin, fun = mean, digits = 0, na.rm = TRUE)
  MIROC_85_2070_mean <-as.data.frame(MIROC_85_2070_mean)
  colnames(MIROC_85_2070_mean)[1] <- "basin_id" #rename zone to basin_id to fit ensemble model
  colnames(MIROC_85_2070_mean)[2:10] <- paste("mean", colnames(MIROC_85_2070_mean[,c(2:10)]), sep = "_")
  write.csv(MIROC_85_2070_mean, "/mnt/brunner/brunner/SDM_data/tables/MIROC_85_2070_mean.csv", row.names = F)
  
  MIROC_85_2070_min <- zonal(MIROC_85_2070, basin, fun = min, digits = 0, na.rm = TRUE)
  MIROC_85_2070_min <-as.data.frame(MIROC_85_2070_min)
  colnames(MIROC_85_2070_min)[1] <- "basin_id"
  colnames(MIROC_85_2070_min)[2:10] <- paste("min", colnames(MIROC_85_2070_min[,c(2:10)]), sep = "_")
  
  MIROC_85_2070_max <- zonal(MIROC_85_2070, basin, fun = max, digits = 0, na.rm = TRUE)
  MIROC_85_2070_max  <-as.data.frame(MIROC_85_2070_max)
  colnames(MIROC_85_2070_max)[1] <- "basin_id"
  colnames(MIROC_85_2070_max)[2:10] <- paste("max", colnames(MIROC_85_2070_max[,c(2:10)]), sep = "_")
  
  MIROC_85_2070_range <- MIROC_85_2070_max %>%  dplyr::mutate(range_MIROC_85_2070_ann_mean_temp = MIROC_85_2070_max[, 2] - MIROC_85_2070_min[, 2], 
                                                              range_MIROC_85_2070_ann_prec = MIROC_85_2070_max[, 3] - MIROC_85_2070_min[, 3],
                                                              range_MIROC_85_2070_mean_temp_coldest_quarter = MIROC_85_2070_max[, 4] - MIROC_85_2070_min[, 4],
                                                              range_MIROC_85_2070_mean_temp_warmest_quarter = MIROC_85_2070_max[, 5] - MIROC_85_2070_min[, 5],
                                                              range_MIROC_85_2070_prec_driest_quarter = MIROC_85_2070_max[, 6] - MIROC_85_2070_min[, 6], 
                                                              range_MIROC_85_2070_prec_seasonality = MIROC_85_2070_max[, 7] - MIROC_85_2070_min[, 7],
                                                              range_MIROC_85_2070_prec_wettest_quarter = MIROC_85_2070_max[, 8] - MIROC_85_2070_min[, 8],
                                                              #temp_range = MIROC_85_2070_max[, 9] - MIROC_85_2070_min[, 9],
                                                              range_MIROC_85_2070_temp_seasonality = MIROC_85_2070_max[, 10] - MIROC_85_2070_min[, 10])
  MIROC_85_2070_range <- MIROC_85_2070_range[-c(2:10)]
  write.csv(MIROC_85_2070_range, "/mnt/brunner/brunner/SDM_data/tables/MIROC_85_2070_range.csv", row.names = F)
  
  
  #=============================== LC percent per basin 2070
  #=========SSP1
  #LANDCOVER SSP1 2070
  LC_SSP1_2070_LAC <- raster("/mnt/brunner/brunner/dinamica/final_model_output/SSP1/SSP1_07.tif")
  LC_SSP1_2070 <- crop(LC_SSP1_2070_LAC, extent(study_area))
  LC_SSP1_2070 <- mask(LC_SSP1_2070, study_area)
  extent(LC_SSP1_2070) <- alignExtent(LC_SSP1_2070, basin)
  #NAvalue(LC_SSP1_2070)<-(-9999)
  #x11(); plot(LC_SSP1_2070)
  
  
  ##pixelcount for the SSP1 LC 2070 map categories per basin
  SSP1_2070_zonal <-myZonal(LC_SSP1_2070, basin)
  results_SSP1_2070 <- data.table::dcast(SSP1_2070_zonal, z ~ vals)
  colnames(results_SSP1_2070)[1] <- "basin_id"
  results_SSP1_2070[is.na(results_SSP1_2070)] <- 0
  
  #bind basin count and LC 10 count
  #all_results <- cbind(results_basin, results_SSP1_2070 )
  SSP1_2070_LC <- merge(results_basin, results_SSP1_2070, by = "basin_id", all= T )
  SSP1_2070_LC <- SSP1_2070_LC[-1,]
  
  
  LC_percentage_SSP1_2070 <- (SSP1_2070_LC[,3:8]/SSP1_2070_LC[,2])
  LC_percentage_SSP1_2070 <- round(LC_percentage_SSP1_2070, digits = 5)
  LC_percentage_SSP1_2070 <- cbind(SSP1_2070_LC[1], LC_percentage_SSP1_2070)
  
  
  LC_percentage_SSP1_2070 <- as.data.frame(LC_percentage_SSP1_2070)
  setnames(LC_percentage_SSP1_2070, old = c('1', '2', '3', '4', '5', '6'), new = c('crops', 'forest','pasture','urban', 'other', 'water'))
  write.csv(LC_percentage_SSP1_2070,"/mnt/brunner/brunner/SDM_data/tables/LC_percentage_SSP1_2070.csv", row.names = F )
  
  #=========SSP2
  #LANDCOVER SSP2 2070
  LC_SSP2_2070_LAC <- raster("/mnt/brunner/brunner/dinamica/final_model_output/SSP2/SSP2_07.tif")
  LC_SSP2_2070 <- crop(LC_SSP2_2070_LAC, extent(study_area))
  LC_SSP2_2070 <- mask(LC_SSP2_2070, study_area)
  extent(LC_SSP2_2070) <- alignExtent(LC_SSP2_2070, basin)
  #NAvalue(LC_SSP2_2070)<-(-9999)
  # x11(); plot(LC_SSP2_2070)
  
  
  ##pixelcount for the SSP2 LC 2070 map categories per basin
  SSP2_2070_zonal <-myZonal(LC_SSP2_2070, basin)
  results_SSP2_2070 <- data.table::dcast(SSP2_2070_zonal, z ~ vals)
  colnames(results_SSP2_2070)[1] <- "basin_id"
  results_SSP2_2070[is.na(results_SSP2_2070)] <- 0
  
  #bind basin count and LC 10 count
  #all_results <- cbind(results_basin, results_SSP2_2070 )
  SSP2_2070_LC <- merge(results_basin, results_SSP2_2070, by = "basin_id", all= T )
  SSP2_2070_LC <- SSP2_2070_LC[-1,]
  
  
  LC_percentage_SSP2_2070 <- (SSP2_2070_LC[,3:8]/SSP2_2070_LC[,2])
  LC_percentage_SSP2_2070 <- round(LC_percentage_SSP2_2070, digits = 5)
  LC_percentage_SSP2_2070 <- cbind(SSP2_2070_LC[1], LC_percentage_SSP2_2070)
  
  
  LC_percentage_SSP2_2070 <- as.data.frame(LC_percentage_SSP2_2070)
  setnames(LC_percentage_SSP2_2070, old = c('1', '2', '3', '4', '5', '6'), new = c('crops', 'forest','pasture','urban', 'other', 'water'))
  write.csv(LC_percentage_SSP2_2070,"/mnt/brunner/brunner/SDM_data/tables/LC_percentage_SSP2_2070.csv", row.names = F )
  
  
  #=========SSP3
  #LANDCOVER SSP3 2070
  LC_SSP3_2070_LAC <- raster("/mnt/brunner/brunner/dinamica/final_model_output/SSP3/SSP3_07.tif")
  LC_SSP3_2070 <- crop(LC_SSP3_2070_LAC, extent(study_area))
  LC_SSP3_2070 <- mask(LC_SSP3_2070, study_area)
  extent(LC_SSP3_2070) <- alignExtent(LC_SSP3_2070, basin)
  #NAvalue(LC_SSP3_2070)<-(-9999)
  #x11(); plot(LC_SSP3_2070)
  
  
  ##pixelcount for the SSP3 LC 2070 map categories per basin
  SSP3_2070_zonal <-myZonal(LC_SSP3_2070, basin)
  results_SSP3_2070 <- data.table::dcast(SSP3_2070_zonal, z ~ vals)
  colnames(results_SSP3_2070)[1] <- "basin_id"
  results_SSP3_2070[is.na(results_SSP3_2070)] <- 0
  
  #bind basin count and LC 10 count
  #all_results <- cbind(results_basin, results_SSP3_2070 )
  SSP3_2070_LC <- merge(results_basin, results_SSP3_2070, by = "basin_id", all= T )
  SSP3_2070_LC <- SSP3_2070_LC[-1,]
  
  
  LC_percentage_SSP3_2070 <- (SSP3_2070_LC[,3:8]/SSP3_2070_LC[,2])
  LC_percentage_SSP3_2070 <- round(LC_percentage_SSP3_2070, digits = 5)
  LC_percentage_SSP3_2070 <- cbind(SSP3_2070_LC[1], LC_percentage_SSP3_2070)
  
  
  LC_percentage_SSP3_2070 <- as.data.frame(LC_percentage_SSP3_2070)
  setnames(LC_percentage_SSP3_2070, old = c('1', '2', '3', '4', '5', '6'), new = c('crops', 'forest','pasture','urban', 'other', 'water'))
  write.csv(LC_percentage_SSP3_2070,"/mnt/brunner/brunner/SDM_data/tables/LC_percentage_SSP3_2070.csv", row.names = F )
  
  #=========SSP4
  #LANDCOVER SSP4 2070
  LC_SSP4_2070_LAC <- raster("/mnt/brunner/brunner/dinamica/final_model_output/SSP4/SSP4_07.tif")
  LC_SSP4_2070 <- crop(LC_SSP4_2070_LAC, extent(study_area))
  LC_SSP4_2070 <- mask(LC_SSP4_2070, study_area)
  extent(LC_SSP4_2070) <- alignExtent(LC_SSP4_2070, basin)
  #NAvalue(LC_SSP4_2070)<-(-9999)
  # x11(); plot(LC_SSP4_2070)
  
  
  ##pixelcount for the SSP4 LC 2070 map categories per basin
  SSP4_2070_zonal <-myZonal(LC_SSP4_2070, basin)
  results_SSP4_2070 <- data.table::dcast(SSP4_2070_zonal, z ~ vals)
  colnames(results_SSP4_2070)[1] <- "basin_id"
  results_SSP4_2070[is.na(results_SSP4_2070)] <- 0
  
  #bind basin count and LC 10 count
  #all_results <- cbind(results_basin, results_SSP4_2070 )
  SSP4_2070_LC <- merge(results_basin, results_SSP4_2070, by = "basin_id", all= T )
  SSP4_2070_LC <- SSP4_2070_LC[-1,]
  
  
  LC_percentage_SSP4_2070 <- (SSP4_2070_LC[,3:8]/SSP4_2070_LC[,2])
  LC_percentage_SSP4_2070 <- round(LC_percentage_SSP4_2070, digits = 5)
  LC_percentage_SSP4_2070 <- cbind(SSP4_2070_LC[1], LC_percentage_SSP4_2070)
  
  
  LC_percentage_SSP4_2070 <- as.data.frame(LC_percentage_SSP4_2070)
  setnames(LC_percentage_SSP4_2070, old = c('1', '2', '3', '4', '5', '6'), new = c('crops', 'forest','pasture','urban', 'other', 'water'))
  write.csv(LC_percentage_SSP4_2070,"/mnt/brunner/brunner/SDM_data/tables/LC_percentage_SSP4_2070.csv", row.names = F )
  
  #=========SSP5
  #LANDCOVER SSP5 2070
  LC_SSP5_2070_LAC <- raster("/mnt/brunner/brunner/dinamica/final_model_output/SSP5/SSP5_07.tif")
  LC_SSP5_2070 <- crop(LC_SSP5_2070_LAC, extent(study_area))
  LC_SSP5_2070 <- mask(LC_SSP5_2070, study_area)
  extent(LC_SSP5_2070) <- alignExtent(LC_SSP5_2070, basin)
  #NAvalue(LC_SSP5_2070)<-(-9999)
  # x11(); plot(LC_SSP5_2070)
  
  
  ##pixelcount for the SSP5 LC 2070 map categories per basin
  SSP5_2070_zonal <-myZonal(LC_SSP5_2070, basin)
  results_SSP5_2070 <- data.table::dcast(SSP5_2070_zonal, z ~ vals)
  colnames(results_SSP5_2070)[1] <- "basin_id"
  results_SSP5_2070[is.na(results_SSP5_2070)] <- 0
  
  #bind basin count and LC 10 count
  #all_results <- cbind(results_basin, results_SSP5_2070 )
  SSP5_2070_LC <- merge(results_basin, results_SSP5_2070, by = "basin_id", all= T )
  SSP5_2070_LC <- SSP5_2070_LC[-1,]
  
  
  LC_percentage_SSP5_2070 <- (SSP5_2070_LC[,3:8]/SSP5_2070_LC[,2])
  LC_percentage_SSP5_2070 <- round(LC_percentage_SSP5_2070, digits = 5)
  LC_percentage_SSP5_2070 <- cbind(SSP5_2070_LC[1], LC_percentage_SSP5_2070)
  
  
  LC_percentage_SSP5_2070 <- as.data.frame(LC_percentage_SSP5_2070)
  setnames(LC_percentage_SSP5_2070, old = c('1', '2', '3', '4', '5', '6'), new = c('crops', 'forest','pasture','urban', 'other', 'water'))
  write.csv(LC_percentage_SSP5_2070,"/mnt/brunner/brunner/SDM_data/tables/LC_percentage_SSP5_2070.csv", row.names = F )
