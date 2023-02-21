if(Sys.info()[["nodename"]]=="vmdomisch04") 
  DIR="/mnt/brunner/brunner/R" 
setwd(DIR)
R_temp_delete <- paste0(DIR, "/R_temp_delete")
rasterOptions(tmpdir=R_temp_delete)

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
if (!require("viridis")) { install.packages("viridis", dependencies = TRUE) ; library(viridis)}
if (!require("RColorBrewer")) { install.packages("RColorBrewer", dependencies = TRUE) ; library(RColorBrewer)}
if (!require("ggplot2")) { install.packages("ggplot2", dependencies = TRUE) ; library(ggplot2)}
if (!require("ggthemes")) { install.packages("ggthemes", dependencies = TRUE) ; library(ggthemes)} # theme_map()
if (!require("caret")) { install.packages("caret", dependencies = TRUE) ; library(caret)} # theme_map()
if (!require("scales")) { install.packages("scales", dependencies = TRUE) ; library(scales)} # scaling
if (!require("tidyr")) { install.packages("tidyr", dependencies = TRUE) ; library(tidyr)} # scaling

#as data table raster function
as.data.table.raster <- function(x, row.names = NULL, optional = FALSE, xy=FALSE, inmem = canProcessInMemory(x, 2), ...) {
  stopifnot(require("data.table"))
  if(inmem) {
    v <- as.data.table(raster::as.data.frame(x, row.names=row.names, optional=optional, xy=xy, ...))
  } else {
    tr <- blockSize(x, n=2)
    l <- lapply(1:tr$n, function(i)
      as.data.table(as.data.frame(getValues(x,
                                            row=tr$row[i],
                                            nrows=tr$nrows[i]),
                                  row.names=row.names, optional=optional, xy=xy, ...)))
    v <- rbindlist(l)
  }
  coln <- names(x)
  if(xy) coln <- c("x", "y", coln)
  setnames(v, coln)
  v
}
if (!isGeneric("as.data.table")) {
  setGeneric("as.data.table", function(x, ...)
    standardGeneric("as.data.table"))
}
setMethod('as.data.table', signature(x='data.frame'), data.table::as.data.table)
setMethod('as.data.table', signature(x='Raster'), as.data.table.raster)
#================================= Data

#=====load environmental data
env_present <- read.csv("/mnt/brunner/brunner/SDM_data/tables/bio_lc_stream_present.csv")

env_SSP1_2050 <- read.csv("/mnt/brunner/brunner/SDM_data/tables/SSP1_env_mean_2050.csv")
env_SSP2_2050 <- read.csv("/mnt/brunner/brunner/SDM_data/tables/SSP2_env_mean_2050.csv")
env_SSP3_2050 <- read.csv("/mnt/brunner/brunner/SDM_data/tables/SSP3_env_mean_2050.csv")
env_SSP4_2050 <- read.csv("/mnt/brunner/brunner/SDM_data/tables/SSP4_env_mean_2050.csv")
env_SSP5_2050 <- read.csv("/mnt/brunner/brunner/SDM_data/tables/SSP5_env_mean_2050.csv")

env_SSP1_2070 <- read.csv("/mnt/brunner/brunner/SDM_data/tables/SSP1_env_mean_2070.csv")
env_SSP2_2070 <- read.csv("/mnt/brunner/brunner/SDM_data/tables/SSP2_env_mean_2070.csv")
env_SSP3_2070 <- read.csv("/mnt/brunner/brunner/SDM_data/tables/SSP3_env_mean_2070.csv")
env_SSP4_2070 <- read.csv("/mnt/brunner/brunner/SDM_data/tables/SSP4_env_mean_2070.csv")
env_SSP5_2070 <- read.csv("/mnt/brunner/brunner/SDM_data/tables/SSP5_env_mean_2070.csv")

#=====load richness data

present_richness_table  <-read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/present_richness_table_bin_dist_CO.csv")


SSP1_2050_richness_table <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/SSP1_2050_richness_table_bin_dist_CO.csv")
SSP2_2050_richness_table <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/SSP2_2050_richness_table_bin_dist_CO.csv")
SSP3_2050_richness_table <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/SSP3_2050_richness_table_bin_dist_CO.csv")
SSP4_2050_richness_table <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/SSP4_2050_richness_table_bin_dist_CO.csv")
SSP5_2050_richness_table <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/SSP5_2050_richness_table_bin_dist_CO.csv")

SSP1_2070_richness_table <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/SSP1_2070_richness_table_bin_dist_CO.csv")
SSP2_2070_richness_table <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/SSP2_2070_richness_table_bin_dist_CO.csv")
SSP3_2070_richness_table <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/SSP3_2070_richness_table_bin_dist_CO.csv")
SSP4_2070_richness_table <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/SSP4_2070_richness_table_bin_dist_CO.csv")
SSP5_2070_richness_table <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/SSP5_2070_richness_table_bin_dist_CO.csv")

#===merge richness & env data
present_richness_env <-merge(present_richness_table, env_present, by = "basin_id")

SSP1_2050_richness_env <-merge(SSP1_2050_richness_table, env_SSP1_2050, by = "basin_id")
SSP2_2050_richness_env <-merge(SSP2_2050_richness_table, env_SSP2_2050, by = "basin_id")
SSP3_2050_richness_env <-merge(SSP3_2050_richness_table, env_SSP3_2050, by = "basin_id")
SSP4_2050_richness_env <-merge(SSP4_2050_richness_table, env_SSP4_2050, by = "basin_id")
SSP5_2050_richness_env <-merge(SSP5_2050_richness_table, env_SSP5_2050, by = "basin_id")

SSP1_2070_richness_env <-merge(SSP1_2070_richness_table, env_SSP1_2070, by = "basin_id")
SSP2_2070_richness_env <-merge(SSP2_2070_richness_table, env_SSP2_2070, by = "basin_id")
SSP3_2070_richness_env <-merge(SSP3_2070_richness_table, env_SSP3_2070, by = "basin_id")
SSP4_2070_richness_env <-merge(SSP4_2070_richness_table, env_SSP4_2070, by = "basin_id")
SSP5_2070_richness_env <-merge(SSP5_2070_richness_table, env_SSP5_2070, by = "basin_id")

#====load binary data for single species
present_binary_dist_all <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/present_richness_all_bin_dist_CO.csv")
species_count <- data.frame(species=names(present_binary_dist_all), basin_sum= colSums(present_binary_dist_all))
write.csv(species_count,"/mnt/brunner/brunner/Species_results/binary_dist_CO/env/count_species_present.csv" )

SSP1_2050_binary_dist_all <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/SSP1_2050_all_means_bin_dist_CO.csv")
SSP2_2050_binary_dist_all <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/SSP2_2050_all_means_bin_dist_CO.csv")
SSP3_2050_binary_dist_all <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/SSP3_2050_all_means_bin_dist_CO.csv")
SSP4_2050_binary_dist_all <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/SSP4_2050_all_means_bin_dist_CO.csv")
SSP5_2050_binary_dist_all <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/SSP5_2050_all_means_bin_dist_CO.csv")

SSP1_2070_binary_dist_all <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/SSP1_2070_all_means_bin_dist_CO.csv")
SSP2_2070_binary_dist_all <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/SSP2_2070_all_means_bin_dist_CO.csv")
SSP3_2070_binary_dist_all <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/SSP3_2070_all_means_bin_dist_CO.csv")
SSP4_2070_binary_dist_all <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/SSP4_2070_all_means_bin_dist_CO.csv")
SSP5_2070_binary_dist_all <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/SSP5_2070_all_means_bin_dist_CO.csv")

#===load means per species
temperature <- env_present[, c("basin_id", "mean_ann_mean_temp")]
temperature$temp_C <- temperature$mean_ann_mean_temp/10

mean_temp <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/mean_sd_temp_species_present.csv")
mean_temp$mean_C <- mean_temp$mean/10
mean_temp$sd_C <- mean_temp$sd/10
mean_temp_SSP2_50 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/mean_sd_temp_species_SSP2_50.csv")
mean_temp_SSP2_50$mean_C <- mean_temp_SSP2_50$mean/10
mean_temp_SSP2_50$sd_C <- mean_temp_SSP2_50$sd/10
mean_temp_SSP5_50 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/mean_sd_temp_species_SSP5_50.csv")
mean_temp_SSP5_50$mean_C <- mean_temp_SSP5_50$mean/10
mean_temp_SSP5_50$sd_C <- mean_temp_SSP5_50$sd/10
mean_temp_SSP2_70 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/mean_sd_temp_species_SSP2_70.csv")
mean_temp_SSP2_70$mean_C <- mean_temp_SSP2_70$mean/10
mean_temp_SSP2_70$sd_C <- mean_temp_SSP2_70$sd/10
mean_temp_SSP5_70 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/mean_sd_temp_species_SSP5_70.csv")
mean_temp_SSP5_70$mean_C <- mean_temp_SSP5_70$mean/10
mean_temp_SSP5_70$sd_C <- mean_temp_SSP5_70$sd/10

#mean temperature

env_species_plot1 <-ggplot(mean_temp,aes(x=reorder(species,mean_C),mean_C))+geom_errorbar(aes(ymin=mean_C-sd_C,ymax=mean_C+sd_C),width=0.2, color="#57bfff")+geom_point()+
  geom_hline(aes(yintercept = mean(mean_C)),color="blue")+
  scale_y_continuous(limits = c(9, 32), breaks = seq(0, 30, by = 5))+
  labs(y="C", x="Species")+
  x11();plot(env_species_plot1)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/temp_species_present_2.png", width = 20, height = 10, units = "cm")

# change basins according to temperature
basin_sum_present <-data.frame(species=names(present_binary_dist_all), basin_sum= colSums(present_binary_dist_all))

SSP2_50_sum_present <-data.frame(species=names(SSP2_2050_binary_dist_all), basin_sum_SSP2_50= colSums(SSP2_2050_binary_dist_all))
SSP5_50_sum_present <-data.frame(species=names(SSP5_2050_binary_dist_all), basin_sum_SSP5_50= colSums(SSP5_2050_binary_dist_all))

SSP2_70_sum_present <-data.frame(species=names(SSP2_2070_binary_dist_all), basin_sum_SSP2_70= colSums(SSP2_2070_binary_dist_all))
SSP5_70_sum_present <-data.frame(species=names(SSP5_2070_binary_dist_all), basin_sum_SSP5_70= colSums(SSP5_2070_binary_dist_all))


mean_temp_basins <- Reduce(merge, list(mean_temp, basin_sum_present, SSP2_50_sum_present, SSP5_50_sum_present, SSP2_70_sum_present, SSP5_70_sum_present))
mean_temp_basins <- mean_temp_basins[-c(2:3, 5)]
mean_temp_basins$change_50_SSP2 <- ((mean_temp_basins$basin_sum_SSP2_50/mean_temp_basins$basin_sum)-1)*100
mean_temp_basins$change_50_SSP5 <- ((mean_temp_basins$basin_sum_SSP5_50/mean_temp_basins$basin_sum)-1)*100

mean_temp_basins$change_70_SSP2 <- ((mean_temp_basins$basin_sum_SSP2_70/mean_temp_basins$basin_sum)-1)*100
mean_temp_basins$change_70_SSP5 <- ((mean_temp_basins$basin_sum_SSP5_70/mean_temp_basins$basin_sum)-1)*100

write.csv(mean_temp_basins,"/mnt/brunner/brunner/Species_results/binary_dist_CO/env/basins_temp_change.csv")


basins_temp_changes2 <- ggplot()+
  geom_point(mean_temp_basins, mapping=aes(x=reorder(species,mean_C),change_50_SSP2, color= "blue"), pch= 20)+
  geom_point(mean_temp_basins, mapping=aes(x=reorder(species,mean_C),change_50_SSP5, color = "#636363"), pch =20)+
  scale_color_identity(name="projection",
                       breaks = c( "blue", "#636363"),
                       labels = c("SSP2 2050", "SSP5 2050"),
                       guide = "legend")+
  labs(y="Change of suitable basin per species (%)", x="Species sorted by base year mean temperature (ascending) ")+
  scale_y_continuous(limits = c(-100, 1000), breaks = seq(-100, 1000, by = 100))+
  x11();plot(basins_temp_changes2)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/basins_temp_change_2050.png", width = 15, height = 10, units = "cm")

basins_temp_changes3 <- ggplot()+
  geom_point(mean_temp_basins, mapping=aes(x=reorder(species,mean_C),change_70_SSP2, color= "blue"), pch= 20)+
  geom_point(mean_temp_basins, mapping=aes(x=reorder(species,mean_C),change_70_SSP5, color = "#636363"), pch =20)+
  scale_color_identity(name="projection",
                       breaks = c( "blue", "#636363"),
                       labels = c("SSP2 2070", "SSP5 2070"),
                       guide = "legend")+
  labs(y="Change of suitable basin per species (%)", x="Species sorted by base year mean temperature (ascending) ")+
  scale_y_continuous(limits = c(-100, 1000), breaks = seq(-100, 1000, by = 100))+
  x11();plot(basins_temp_changes3)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/basins_temp_change_2070.png", width = 15, height = 10, units = "cm")


# precipitation 

mean_prec <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/mean_sd_prec_species_present.csv")
mean_prec_SSP2_50 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/mean_sd_prec_species_SSP2_50.csv")
mean_prec_SSP5_50 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/mean_sd_prec_species_SSP5_50.csv")
mean_prec_SSP2_70 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/mean_sd_prec_species_SSP2_70.csv")
mean_prec_SSP5_70 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/mean_sd_prec_species_SSP5_70.csv")

env_species_plot2 <-ggplot(mean_prec,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#2c6b92")+geom_point()+
  geom_hline(aes(yintercept = mean(mean)),color="blue")+
  scale_y_continuous(limits = c(500, 8000), breaks = seq(500, 8000, by = 500))+
  labs(y="Kg m^-2", x="Species")+
  x11();plot(env_species_plot2)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/prec_species_present_2.png", width = 20, height = 10, units = "cm")




basin_sum_present <-data.frame(species=names(present_binary_dist_all), basin_sum= colSums(present_binary_dist_all))

SSP2_50_sum_present <-data.frame(species=names(SSP2_2050_binary_dist_all), basin_sum_SSP2_50= colSums(SSP2_2050_binary_dist_all))
SSP5_50_sum_present <-data.frame(species=names(SSP5_2050_binary_dist_all), basin_sum_SSP5_50= colSums(SSP5_2050_binary_dist_all))

SSP2_70_sum_present <-data.frame(species=names(SSP2_2070_binary_dist_all), basin_sum_SSP2_70= colSums(SSP2_2070_binary_dist_all))
SSP5_70_sum_present <-data.frame(species=names(SSP5_2070_binary_dist_all), basin_sum_SSP5_70= colSums(SSP5_2070_binary_dist_all))


mean_prec_basins <- Reduce(merge, list(mean_prec, basin_sum_present, SSP2_50_sum_present, SSP5_50_sum_present, SSP2_70_sum_present, SSP5_70_sum_present))
mean_prec_basins <- mean_prec_basins[-c(3)]
mean_prec_basins$change_50_SSP2 <- ((mean_prec_basins$basin_sum_SSP2_50/mean_prec_basins$basin_sum)-1)*100
mean_prec_basins$change_50_SSP5 <- ((mean_prec_basins$basin_sum_SSP5_50/mean_prec_basins$basin_sum)-1)*100

mean_prec_basins$change_70_SSP2 <- ((mean_prec_basins$basin_sum_SSP2_70/mean_prec_basins$basin_sum)-1)*100
mean_prec_basins$change_70_SSP5 <- ((mean_prec_basins$basin_sum_SSP5_70/mean_prec_basins$basin_sum)-1)*100

basins_prec_changes2 <- ggplot()+
  geom_point(mean_prec_basins, mapping=aes(x=reorder(species,mean),change_50_SSP2, color= "blue"), pch= 20)+
  geom_point(mean_prec_basins, mapping=aes(x=reorder(species,mean),change_50_SSP5, color = "#636363"), pch =20)+
  scale_color_identity(name="projection",
                       breaks = c( "blue", "#636363"),
                       labels = c("SSP2 2050", "SSP5 2050"),
                       guide = "legend")+
  labs(y="Change of suitable basin per species (%)", x="Species sorted by base year mean precipitation (ascending) ")+
  scale_y_continuous(limits = c(-100, 2000), breaks = seq(-100, 2000, by = 100))+
  x11();plot(basins_prec_changes2)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/basins_prec_change_2050.png", width = 15, height = 10, units = "cm")

basins_prec_changes3 <- ggplot()+
  geom_point(mean_prec_basins, mapping=aes(x=reorder(species,mean),change_70_SSP2, color= "blue"), pch= 20)+
  geom_point(mean_prec_basins, mapping=aes(x=reorder(species,mean),change_70_SSP5, color = "#636363"), pch =20)+
  scale_color_identity(name="projection",
                       breaks = c( "blue", "#636363"),
                       labels = c("SSP2 2070", "SSP5 2070"),
                       guide = "legend")+
  labs(y="Change of suitable basin per species (%)", x="Species sorted by base year mean precipitation (ascending) ")+
  scale_y_continuous(limits = c(-100, 2000), breaks = seq(-100, 2000, by = 100))+
  x11();plot(basins_prec_changes3)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/basins_prec_change_2070.png", width = 15, height = 10, units = "cm")



#temp & prec
temp_prec <- cbind(mean_temp [,1:2], mean_prec=mean_prec[,2])

env_species_plot3 <-ggplot()+ 
  geom_point(temp_prec,mapping=aes(x=reorder(species,mean),mean_prec), color ="blue", pch= 20)+
  # scale_y_continuous(limits = c(500, 8000), breaks = seq(500, 8000, by = 500))+
  labs(y="Kg m^-2", x="Species sorted by base year mean temperature (ascending)")+
  x11();plot(env_species_plot3)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/prec_species_present_2.png", width = 20, height = 10, units = "cm")

#mean temperature and elevation change

env_present_not_used <- read.csv("/mnt/brunner/brunner/SDM_data/tables/all_var/topo_bio_lc_stream_present.csv")
#ELEVATION
present_richness_elev <-merge(present_richness_table, env_present_not_used[,c("basin_id","mean_elevation_global")], by = "basin_id")
SSP2_2050_richness_elev <-merge(SSP2_2050_richness_table, env_present_not_used[,c("basin_id","mean_elevation_global")], by = "basin_id")
SSP2_2070_richness_elev <-merge(SSP2_2070_richness_table, env_present_not_used[,c("basin_id","mean_elevation_global")], by = "basin_id")
SSP5_2050_richness_elev <-merge(SSP5_2050_richness_table, env_present_not_used[,c("basin_id","mean_elevation_global")], by = "basin_id")
SSP5_2070_richness_elev <-merge(SSP5_2070_richness_table, env_present_not_used[,c("basin_id","mean_elevation_global")], by = "basin_id")


mean_elev <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/mean_sd_elev_species_present.csv")
mean_elev_SSP2_50 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/mean_sd_elev_species_SSP2_50.csv")
mean_elev_SSP5_50 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/mean_sd_elev_species_SSP5_50.csv")
mean_elev_SSP2_70 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/mean_sd_elev_species_SSP2_70.csv")
mean_elev_SSP5_70 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/mean_sd_elev_species_SSP5_70.csv")



elev_change <- cbind(mean_elev[,1:2],mean_temp=mean_temp[,2], mean_prec=mean_prec[,2], mean_elev_SSP2_50=mean_elev_SSP2_50[,2], mean_elev_SSP5_50= mean_elev_SSP5_50[,2],
                     mean_elev_SSP2_70=mean_elev_SSP2_70[,2], mean_elev_SSP5_70= mean_elev_SSP5_70[,2])
elev_change$SSP2_50 <- elev_change$mean_elev_SSP2_50 - elev_change$mean
elev_change$SSP5_50 <- elev_change$mean_elev_SSP5_50 - elev_change$mean
elev_change$SSP2_70 <- elev_change$mean_elev_SSP2_70 - elev_change$mean
elev_change$SSP5_70 <- elev_change$mean_elev_SSP5_70 - elev_change$mean

write.csv(elev_change, "/mnt/brunner/brunner/Species_results/binary_dist_CO/env/elev_change.csv", row.names = F)


species_elev_change4 <- ggplot()+
  geom_point(elev_change, mapping=aes(x=reorder(species, mean_temp), SSP2_50, color="blue"), pch =20)+
  geom_point(elev_change, mapping=aes(x=reorder(species, mean_temp), SSP5_50, color= "#636363"), pch=20)+
  scale_color_identity(name="projection",
                       breaks = c( "blue", "#636363"),
                       labels = c("SSP2 2050", "SSP5 2050"),
                       guide = "legend")+
  scale_y_continuous(limits = c(-500, 1000), breaks = seq(-500, 1000, by = 100))+
  labs(y="Change in mean elevation(m)", x="Species sorted by base year mean temperature (ascending)")+
  x11();plot(species_elev_change4)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/elev_change_temp_2050.png", width = 15, height = 10, units = "cm")


species_elev_change5 <- ggplot()+
  geom_point(elev_change, mapping=aes(x=reorder(species, mean_temp), SSP2_70, color="blue"),pch =20)+
  geom_point(elev_change, mapping=aes(x=reorder(species, mean_temp), SSP5_70, color= "#636363"),pch =20)+
  scale_color_identity(name="projection",
                       breaks = c( "blue", "#636363"),
                       labels = c("SSP2 2070", "SSP5 2070"),
                       guide = "legend")+
  scale_y_continuous(limits = c(-500, 1000), breaks = seq(-500, 1000, by = 100))+
  labs(y="Change in mean elevation(m)", x="Species sorted by base year mean temperature (ascending)")+
  x11();plot(species_elev_change5)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/elev_change_temp_2070.png", width = 15, height = 10, units = "cm")


species_elev_change6 <- ggplot()+
  geom_point(elev_change, mapping=aes(x=reorder(species, mean_prec), SSP2_50, color="blue"),pch =20)+
  geom_point(elev_change, mapping=aes(x=reorder(species, mean_prec), SSP5_50, color= "#636363"),pch =20)+
  scale_color_identity(name="projection",
                       breaks = c( "blue", "#636363"),
                       labels = c("SSP2 2050", "SSP5 2050"),
                       guide = "legend")+
  scale_y_continuous(limits = c(-500, 1000), breaks = seq(-500, 1000, by = 100))+
  labs(y="Change in mean elevation(m)", x="Species sorted by base year mean precipitation (ascending)")+
  x11();plot(species_elev_change6)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/elev_change_prec_2050.png", width = 15, height = 10, units = "cm")


species_elev_change7 <- ggplot()+
  geom_point(elev_change, mapping=aes(x=reorder(species, mean_prec), SSP2_70, color="blue"),pch =20)+
  geom_point(elev_change, mapping=aes(x=reorder(species, mean_prec), SSP5_70, color= "#636363"),pch =20)+
  scale_color_identity(name="projection",
                       breaks = c( "blue", "#636363"),
                       labels = c("SSP2 2070", "SSP5 2070"),
                       guide = "legend")+
  scale_y_continuous(limits = c(-500, 1000), breaks = seq(-500, 1000, by = 100))+
  labs(y="Change in mean elevation(m)", x="Species sorted by base year mean precipitation (ascending)")+
  x11();plot(species_elev_change7)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/elev_change_prec_2070.png", width = 15, height = 10, units = "cm")


#mean temperature and flow accumulation
mean_flow_accum <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/mean_sd_flow_accum_species_present.csv")
mean_flow_accum_SSP2_50 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/mean_sd_flow_accum_species_SSP2_50.csv")
mean_flow_accum_SSP5_50 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/mean_sd_flow_accum_species_SSP5_50.csv")
mean_flow_accum_SSP2_70 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/mean_sd_flow_accum_species_SSP2_70.csv")
mean_flow_accum_SSP5_70 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/mean_sd_flow_accum_species_SSP5_70.csv")

flow_accum_change <- cbind(mean_flow_accum[,1:2],mean_temp=mean_temp[,2], mean_prec=mean_prec[,2], mean_flow_accum_SSP2_50=mean_flow_accum_SSP2_50[,2], mean_flow_accum_SSP5_50= mean_flow_accum_SSP5_50[,2],
                           mean_flow_accum_SSP2_70=mean_flow_accum_SSP2_70[,2], mean_flow_accum_SSP5_70= mean_flow_accum_SSP5_70[,2])
flow_accum_change$SSP2_50 <- flow_accum_change$mean_flow_accum_SSP2_50 - flow_accum_change$mean
flow_accum_change$SSP5_50 <- flow_accum_change$mean_flow_accum_SSP5_50 - flow_accum_change$mean
flow_accum_change$SSP2_70 <- flow_accum_change$mean_flow_accum_SSP2_70 - flow_accum_change$mean
flow_accum_change$SSP5_70 <- flow_accum_change$mean_flow_accum_SSP5_70 - flow_accum_change$mean





species_flow_accum_change4 <- ggplot()+
  geom_point(flow_accum_change, mapping=aes(x=reorder(species, mean_temp), SSP2_50, color="blue"),pch =20)+
  geom_point(flow_accum_change, mapping=aes(x=reorder(species, mean_temp), SSP5_50, color= "#636363"),pch =20)+
  scale_color_identity(name="projection",
                       breaks = c( "blue", "#636363"),
                       labels = c("SSP2 2050", "SSP5 2050"),
                       guide = "legend")+
  labs(y="Change in mean flow accumulation", x="Species sorted by base year mean temperature (ascending)")+
  
  x11(); plot(species_flow_accum_change4)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/flow_accum_change_temp_2050.png", width = 15, height = 10, units = "cm")

species_flow_accum_change5 <- ggplot()+
  geom_point(flow_accum_change, mapping=aes(x=reorder(species, mean_temp), SSP2_70, color="blue"),pch =20)+
  geom_point(flow_accum_change, mapping=aes(x=reorder(species, mean_temp), SSP5_70, color= "#636363"),pch =20)+
  scale_color_identity(name="projection",
                       breaks = c( "blue", "#636363"),
                       labels = c("SSP2 2070", "SSP5 2070"),
                       guide = "legend")+
  labs(y="Change in mean flow accumulation", x="Species sorted by base year mean temperature (ascending)")+
  
  x11(); plot(species_flow_accum_change5)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/flow_accum_change_temp_2070.png", width = 15, height = 10, units = "cm")



species_flow_accum_change4 <- ggplot()+
  geom_point(flow_accum_change, mapping=aes(x=reorder(species, mean_prec), SSP2_50, color="blue"),pch =20)+
  geom_point(flow_accum_change, mapping=aes(x=reorder(species, mean_prec), SSP5_50, color= "#636363"),pch =20)+
  scale_color_identity(name="projection",
                       breaks = c( "blue", "#636363"),
                       labels = c("SSP2 2050", "SSP5 2050"),
                       guide = "legend")+
  labs(y="Change in mean flow accumulation", x="Species sorted by base year mean precipitation (ascending)")+
  
  x11(); plot(species_flow_accum_change4)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/flow_accum_change_prec_2050.png", width = 15, height = 10, units = "cm")

species_flow_accum_change4 <- ggplot()+
  geom_point(flow_accum_change, mapping=aes(x=reorder(species, mean_prec), SSP2_70, color="blue"),pch =20)+
  geom_point(flow_accum_change, mapping=aes(x=reorder(species, mean_prec), SSP5_70, color= "#636363"),pch =20)+
  scale_color_identity(name="projection",
                       breaks = c( "blue", "#636363"),
                       labels = c("SSP2 2070", "SSP5 2070"),
                       guide = "legend")+
  labs(y="Change in mean flow accumulation", x="Species sorted by base year mean precipitation (ascending)")+
  
  x11(); plot(species_flow_accum_change4)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/flow_accum_change_prec_2070.png", width = 15, height = 10, units = "cm")


#forest
mean_forest <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/mean_sd_forest_species_present.csv")
mean_forest$mean_pcnt <- mean_forest$mean*100
mean_forest$sd_pcnt <-mean_forest$sd*100
mean_forest_SSP2_50 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/mean_sd_forest_species_SSP2_50.csv")
mean_forest_SSP2_50$mean_pcnt <- mean_forest_SSP2_50$mean*100
mean_forest_SSP2_50$sd_pcnt <-mean_forest_SSP2_50$sd*100
mean_forest_SSP5_50 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/mean_sd_forest_species_SSP5_50.csv")
mean_forest_SSP5_50$mean_pcnt <- mean_forest_SSP5_50$mean*100
mean_forest_SSP5_50$sd_pcnt <-mean_forest_SSP5_50$sd*100
mean_forest_SSP2_70 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/mean_sd_forest_species_SSP2_70.csv")
mean_forest_SSP2_70$mean_pcnt <- mean_forest_SSP2_70$mean*100
mean_forest_SSP2_70$sd_pcnt <-mean_forest_SSP2_70$sd*100
mean_forest_SSP5_70 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/mean_sd_forest_species_SSP5_70.csv")
mean_forest_SSP5_70$mean_pcnt <- mean_forest_SSP5_70$mean*100
mean_forest_SSP5_70$sd_pcnt <-mean_forest_SSP5_70$sd*100

forest_species_plot1 <-ggplot(mean_forest,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#099a1a")+geom_point()+
  geom_hline(aes(yintercept = mean(mean)),color="blue")+
  scale_y_continuous(limits = c(-0.15, 1.1), breaks = seq(0, 1, by = 0.1))+
  labs(y="forest cover (%)", x="Species")+
  ggtitle("Forest cover of species occurences- present")
x11();plot(forest_species_plot1)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/forest_species_present.png", width = 20, height = 10, units = "cm")

forest_species_plot2 <-ggplot(mean_forest,aes(x=reorder(species,mean_pcnt),mean_pcnt))+
  geom_errorbar(aes(ymin=mean_pcnt-sd_pcnt,ymax=mean_pcnt+sd_pcnt),width=0.2, color="#099a1a")+geom_point()+
  geom_hline(aes(yintercept = mean(mean_pcnt)),color="blue")+
  scale_y_continuous(limits = c(-15, 110), breaks = seq(-10, 100, by = 10))+
  labs(y="forest cover (%)", x="Species")+
  ggtitle("Forest cover of species occurences- present")
x11();plot(forest_species_plot2)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/forest_species_present_pcnt.png", width = 20, height = 10, units = "cm")




# change basins according to forest
basin_sum_present <-data.frame(species=names(present_binary_dist_all), basin_sum= colSums(present_binary_dist_all))

SSP2_50_sum_present <-data.frame(species=names(SSP2_2050_binary_dist_all), basin_sum_SSP2_50= colSums(SSP2_2050_binary_dist_all))
SSP5_50_sum_present <-data.frame(species=names(SSP5_2050_binary_dist_all), basin_sum_SSP5_50= colSums(SSP5_2050_binary_dist_all))

SSP2_70_sum_present <-data.frame(species=names(SSP2_2070_binary_dist_all), basin_sum_SSP2_70= colSums(SSP2_2070_binary_dist_all))
SSP5_70_sum_present <-data.frame(species=names(SSP5_2070_binary_dist_all), basin_sum_SSP5_70= colSums(SSP5_2070_binary_dist_all))


mean_forest_basins <- Reduce(merge, list(mean_forest, basin_sum_present, SSP2_50_sum_present, SSP5_50_sum_present, SSP2_70_sum_present, SSP5_70_sum_present))
mean_forest_basins <- mean_forest_basins[-c(2:3,5)]
mean_forest_basins$change_50_SSP2 <- ((mean_forest_basins$basin_sum_SSP2_50/mean_forest_basins$basin_sum)-1)*100
mean_forest_basins$change_50_SSP5 <- ((mean_forest_basins$basin_sum_SSP5_50/mean_forest_basins$basin_sum)-1)*100

mean_forest_basins$change_70_SSP2 <- ((mean_forest_basins$basin_sum_SSP2_70/mean_forest_basins$basin_sum)-1)*100
mean_forest_basins$change_70_SSP5 <- ((mean_forest_basins$basin_sum_SSP5_70/mean_forest_basins$basin_sum)-1)*100

write.csv(mean_forest_basins,"/mnt/brunner/brunner/Species_results/binary_dist_CO/env/basins_forest_change.csv")



basins_forest_changes7 <- ggplot()+
  geom_point(mean_forest_basins, mapping=aes(x=reorder(species,mean_pcnt),change_50_SSP2, color="blue"), pch=20)+
  geom_point(mean_forest_basins, mapping=aes(x=reorder(species,mean_pcnt),change_50_SSP5, color="#636363"), pch=20)+
  scale_color_identity(name="projection",
                       breaks = c( "blue", "#636363"),
                       labels = c("SSP2 2050", "SSP5 2050"),
                       guide = "legend")+
  labs(y="Change of suitable basin per species (%)", x="Species sorted by forest cover of the base year (ascending)")+
  scale_y_continuous(limits = c(-100, 1000), breaks = seq(-100, 1000, by = 100))+
  x11();plot(basins_forest_changes7)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/basins_forest_change_2050.png", width = 15, height = 10, units = "cm")

basins_forest_changes8 <- ggplot()+
  geom_point(mean_forest_basins, mapping=aes(x=reorder(species,mean_pcnt),change_70_SSP2, color="blue"), pch=20)+
  geom_point(mean_forest_basins, mapping=aes(x=reorder(species,mean_pcnt),change_70_SSP5, color="#636363"),pch=20)+
  scale_color_identity(name="projection",
                       breaks = c( "blue", "#636363"),
                       labels = c("SSP2 2070", "SSP5 2070"),
                       guide = "legend")+
  labs(y="Change of suitable basin per species (%)", x="Species sorted by forest cover of the base year (ascending)")+
  scale_y_continuous(limits = c(-100, 1000), breaks = seq(-100, 1000, by = 100))+
  x11();plot(basins_forest_changes8)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/basins_forest_change_2070.png", width = 15, height =10, units = "cm")




forest_change <- cbind(mean_forest[,1:2],mean_temp=mean_temp[,2],mean_elev=mean_elev[,2], mean_flow=mean_flow_accum[,2], mean_forest_SSP2_50=mean_forest_SSP2_50[,2], mean_forest_SSP5_50= mean_forest_SSP5_50[,2],
                       mean_forest_SSP2_70=mean_forest_SSP2_70[,2], mean_forest_SSP5_70= mean_forest_SSP5_70[,2])
forest_change$SSP2_50 <- (forest_change$mean_forest_SSP2_50 - forest_change$mean)*100
forest_change$SSP5_50 <- (forest_change$mean_forest_SSP5_50 - forest_change$mean)*100
forest_change$SSP2_70 <- (forest_change$mean_forest_SSP2_70 - forest_change$mean)*100
forest_change$SSP5_70 <- (forest_change$mean_forest_SSP5_70 - forest_change$mean)*100
write.csv("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/basins_forest_change_SSP2_1.png")
nr_positv <- count_if(gt(0), forest_change$SSP2_50)
nr_negativ <- count_if(lt(0), forest_change$SSP2_50)

nr_positv5 <- count_if(gt(0), forest_change$SSP5_50)
nr_negativ5 <- count_if(lt(0), forest_change$SSP5_50)




species_forest_change3 <- ggplot()+
  geom_point(forest_change, mapping=aes(x=reorder(species, mean_elev), SSP2_50, color="blue"), pch=20)+
  geom_point(forest_change, mapping=aes(x=reorder(species, mean_elev), SSP5_50, color= "#636363"), pch=20)+
  scale_color_identity(name="projection",
                       breaks = c( "blue", "#636363"),
                       labels = c("SSP2 2050", "SSP5 2050"),
                       guide = "legend")+
  scale_y_continuous(limits = c(-45, 50), breaks = seq(-40, 50, by = 10))+
  labs(y="Change in mean forest cover per species (%) ", x="Species sorted by base year elevation (ascending)")+
  x11(); plot(species_forest_change3)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/forest_elevation_change_2050.png", width = 15, height = 10, units = "cm")



species_forest_change4 <- ggplot()+
  geom_point(forest_change, mapping=aes(x=reorder(species, mean_elev), SSP2_70, color="blue"), pch=20)+
  geom_point(forest_change, mapping=aes(x=reorder(species, mean_elev), SSP5_70, color= "#636363"), pch=20)+
  scale_color_identity(name="projection",
                       breaks = c( "blue", "#636363"),
                       labels = c("SSP2 2070", "SSP5 2070"),
                       guide = "legend")+
  scale_y_continuous(limits = c(-45, 50), breaks = seq(-40, 50, by = 10))+
  labs(y="Change in mean forest cover per species (%) ", x="Species sorted by base year elevation (ascending)")+
  x11(); plot(species_forest_change4)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/forest_elevation_change_2070.png", width = 15, height = 10, units = "cm")

species_forest_change5 <- ggplot()+
  geom_point(forest_change, mapping=aes(x=reorder(species, mean_temp), SSP2_50, color="blue"), pch=20)+
  geom_point(forest_change, mapping=aes(x=reorder(species, mean_temp), SSP5_50, color= "#636363"), pch=20)+
  scale_color_identity(name="projection",
                       breaks = c( "blue", "#636363"),
                       labels = c("SSP2 2050", "SSP5 2050"),
                       guide = "legend")+
  scale_y_continuous(limits = c(-45, 50), breaks = seq(-40, 50, by = 10))+
  labs(y="Change in mean forest cover per species (%) ", x="Species sorted by base year temperature (ascending)")+
  x11(); plot(species_forest_change5)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/forest_temperature_change_2050.png", width = 15, height = 10, units = "cm")



species_forest_change6 <- ggplot()+
  geom_point(forest_change, mapping=aes(x=reorder(species, mean_temp), SSP2_70, color="blue"), pch=20)+
  geom_point(forest_change, mapping=aes(x=reorder(species, mean_temp), SSP5_70, color= "#636363"), pch=20)+
  scale_color_identity(name="projection",
                       breaks = c( "blue", "#636363"),
                       labels = c("SSP2 2070", "SSP5 2070"),
                       guide = "legend")+
  scale_y_continuous(limits = c(-45, 50), breaks = seq(-40, 50, by = 10))+
  labs(y="Change in mean forest cover per species (%) ", x="Species sorted by base year temperature (ascending)")+
  x11(); plot(species_forest_change6)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/forest_temperature_change_2070.png", width = 15, height = 10, units = "cm")

species_forest_change7 <- ggplot()+
  geom_point(forest_change, mapping=aes(x=reorder(species, mean_flow), SSP2_50, color="blue"), pch=20)+
  geom_point(forest_change, mapping=aes(x=reorder(species, mean_flow), SSP5_50, color= "#636363"), pch=20)+
  scale_color_identity(name="projection",
                       breaks = c( "blue", "#636363"),
                       labels = c("SSP2 2050", "SSP5 2050"),
                       guide = "legend")+
  scale_y_continuous(limits = c(-45, 50), breaks = seq(-40, 50, by = 10))+
  labs(y="Change in mean forest cover per species (%) ", x="Species sorted by base year flow accumulation (ascending)")+
  x11(); plot(species_forest_change7)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/forest_flow accumulation_change_2050.png", width = 15, height = 10, units = "cm")



species_forest_change8 <- ggplot()+
  geom_point(forest_change, mapping=aes(x=reorder(species, mean_flow), SSP2_70, color="blue"), pch=20)+
  geom_point(forest_change, mapping=aes(x=reorder(species, mean_flow), SSP5_70, color= "#636363"), pch=20)+
  scale_color_identity(name="projection",
                       breaks = c( "blue", "#636363"),
                       labels = c("SSP2 2070", "SSP5 2070"),
                       guide = "legend")+
  scale_y_continuous(limits = c(-45, 50), breaks = seq(-40, 50, by = 10))+
  labs(y="Change in mean forest cover per species (%) ", x="Species sorted by base year flow accumulation (ascending)")+
  x11(); plot(species_forest_change8)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_CO/env/forest_flow accumulation_change_2070.png", width = 15, height = 10, units = "cm")




