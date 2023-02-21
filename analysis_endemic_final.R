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

present_richness_table  <-read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/present_richness_table_bin_dist_endemic.csv")

SSP1_2050_richness_table <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP1_2050_richness_table_bin_dist_endemic.csv")
SSP2_2050_richness_table <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP2_2050_richness_table_bin_dist_endemic.csv")
SSP3_2050_richness_table <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP3_2050_richness_table_bin_dist_endemic.csv")
SSP4_2050_richness_table <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP4_2050_richness_table_bin_dist_endemic.csv")
SSP5_2050_richness_table <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP5_2050_richness_table_bin_dist_endemic.csv")

SSP1_2070_richness_table <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP1_2070_richness_table_bin_dist_endemic.csv")
SSP2_2070_richness_table <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP2_2070_richness_table_bin_dist_endemic.csv")
SSP3_2070_richness_table <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP3_2070_richness_table_bin_dist_endemic.csv")
SSP4_2070_richness_table <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP4_2070_richness_table_bin_dist_endemic.csv")
SSP5_2070_richness_table <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP5_2070_richness_table_bin_dist_endemic.csv")

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
present_binary_dist_all <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/present_richness_all_bin_dist_endemic.csv")

SSP1_2050_binary_dist_all <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP1_2050_all_means_bin_dist_endemic.csv")
SSP2_2050_binary_dist_all <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP2_2050_all_means_bin_dist_endemic.csv")
SSP3_2050_binary_dist_all <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP3_2050_all_means_bin_dist_endemic.csv")
SSP4_2050_binary_dist_all <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP4_2050_all_means_bin_dist_endemic.csv")
SSP5_2050_binary_dist_all <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP5_2050_all_means_bin_dist_endemic.csv")

SSP1_2070_binary_dist_all <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP1_2070_all_means_bin_dist_endemic.csv")
SSP2_2070_binary_dist_all <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP2_2070_all_means_bin_dist_endemic.csv")
SSP3_2070_binary_dist_all <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP3_2070_all_means_bin_dist_endemic.csv")
SSP4_2070_binary_dist_all <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP4_2070_all_means_bin_dist_endemic.csv")
SSP5_2070_binary_dist_all <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP5_2070_all_means_bin_dist_endemic.csv")

#================================= count total number of basin with probability
library(expss)
nr_basins_present <- count_if(gt(0), present_richness_table$present_richness)

nr_basins_SSP1_2050 <- count_if(gt(0), SSP1_2050_richness_table$richness_SSP1_2050) 
nr_basins_SSP2_2050 <- count_if(gt(0), SSP2_2050_richness_table$richness_SSP2_2050) 
nr_basins_SSP3_2050 <- count_if(gt(0), SSP3_2050_richness_table$richness_SSP3_2050)
nr_basins_SSP4_2050 <- count_if(gt(0), SSP4_2050_richness_table$richness_SSP4_2050) 
nr_basins_SSP5_2050 <- count_if(gt(0), SSP5_2050_richness_table$richness_SSP5_2050)


nr_basins_SSP1_2070 <- count_if(gt(0), SSP1_2070_richness_table$richness_SSP1_2070) 
nr_basins_SSP2_2070 <- count_if(gt(0), SSP2_2070_richness_table$richness_SSP2_2070) 
nr_basins_SSP3_2070 <- count_if(gt(0), SSP3_2070_richness_table$richness_SSP3_2070) 
nr_basins_SSP4_2070 <- count_if(gt(0), SSP4_2070_richness_table$richness_SSP4_2070) 
nr_basins_SSP5_2070 <- count_if(gt(0), SSP5_2070_richness_table$richness_SSP5_2070) 


nr_basins <- as.data.frame(rbind(nr_basins_present,nr_basins_SSP1_2050,  nr_basins_SSP2_2050, nr_basins_SSP3_2050, nr_basins_SSP4_2050, nr_basins_SSP5_2050,
                                 nr_basins_SSP1_2070, nr_basins_SSP2_2070, nr_basins_SSP3_2070, nr_basins_SSP4_2070, nr_basins_SSP5_2070))

write.csv(nr_basins, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/nr_basins_bin_dist_endemic.csv")



#================================ which species are missing
extinct_SSP1_2050 <- names(which(colSums(SSP1_2050_binary_dist_all ==1) == 0))
write.csv(extinct_SSP1_2050, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/gone_species/extinct_SSP1_2050.csv", row.names = F)
extinct_SSP2_2050 <- names(which(colSums(SSP2_2050_binary_dist_all ==1) == 0))
write.csv(extinct_SSP2_2050, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/gone_species/extinct_SSP2_2050.csv", row.names = F)
extinct_SSP3_2050 <- names(which(colSums(SSP3_2050_binary_dist_all ==1) == 0))
write.csv(extinct_SSP3_2050, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/gone_species/extinct_SSP3_2050.csv", row.names = F)
extinct_SSP4_2050 <- names(which(colSums(SSP4_2050_binary_dist_all ==1) == 0))
write.csv(extinct_SSP4_2050, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/gone_species/extinct_SSP4_2050.csv", row.names = F)
extinct_SSP5_2050 <- names(which(colSums(SSP5_2050_binary_dist_all ==1) == 0))
write.csv(extinct_SSP5_2050, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/gone_species/extinct_SSP5_2050.csv", row.names = F)

extinct_SSP1_2070 <- names(which(colSums(SSP1_2070_binary_dist_all ==1) == 0))
write.csv(extinct_SSP1_2070, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/gone_species/extinct_SSP1_2070.csv", row.names = F)
extinct_SSP2_2070 <- names(which(colSums(SSP2_2070_binary_dist_all ==1) == 0))
write.csv(extinct_SSP2_2070, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/gone_species/extinct_SSP2_2070.csv", row.names = F)
extinct_SSP3_2070 <- names(which(colSums(SSP3_2070_binary_dist_all ==1) == 0))
write.csv(extinct_SSP3_2070, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/gone_species/extinct_SSP3_2070.csv", row.names = F)
extinct_SSP4_2070 <- names(which(colSums(SSP4_2070_binary_dist_all ==1) == 0))
write.csv(extinct_SSP4_2070, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/gone_species/extinct_SSP4_2070.csv", row.names = F)
extinct_SSP5_2070 <- names(which(colSums(SSP5_2070_binary_dist_all ==1) == 0))
write.csv(extinct_SSP5_2070, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/gone_species/extinct_SSP5_2070.csv", row.names = F)

#================================ stats
richness_lists <- Reduce(merge, list(present_richness_table, SSP1_2050_richness_table,SSP2_2050_richness_table,SSP3_2050_richness_table,SSP4_2050_richness_table,SSP5_2050_richness_table,
                                     SSP1_2070_richness_table,SSP2_2070_richness_table,SSP3_2070_richness_table,SSP4_2070_richness_table, SSP5_2070_richness_table))
richness_lists<- richness_lists[-c(1)]
x <- as.data.frame(richness_lists)

stats <-sapply(x, function(x) c( "Stand dev" = sd(x), 
                         "Mean"= mean(x,na.rm=TRUE),
                         "n" = length(x),
                         "Median" = median(x),
                         "CoeffofVariation" = sd(x)/mean(x,na.rm=TRUE),
                         "Minimum" = min(x),
                         "Maximun" = max(x),
                         "Upper Quantile" = quantile(x,1),
                         "LowerQuartile" = quantile(x,0)
)
)

write.csv(stats, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/richness_stats.csv", row.names = T)

#==========Temperature

#density of temp
ann_mean_temp_present_SSP2_SSP5_50_70_density <-ggplot()+
  geom_density(present_richness_env,  mapping= aes(x=mean_ann_mean_temp,  color="#004242"))+
  geom_density(SSP2_2050_richness_env, mapping=aes(x=mean_ann_mean_temp, color="blue"))+
  geom_density(SSP5_2050_richness_env, mapping=aes(x=mean_ann_mean_temp, color="green"))+
  geom_density(SSP2_2070_richness_env , mapping=aes(x=mean_ann_mean_temp, color="#00c2c2"))+
  geom_density(SSP5_2070_richness_env , mapping=aes(x=mean_ann_mean_temp, color="#0f5b2f"))+
  scale_color_identity(name="projection",
                       breaks = c("#004242", "blue", "green", "#00c2c2", "#0f5b2f"),
                       labels = c("present", "SSP2_2050", "SSP5_2050", "SSP2_2070", "SSP5_2070"),
                       guide = "legend")+
  labs(y="density", x="ann mean temp")+     
  ggtitle("Temperature density")+
  #x11();plot(ann_mean_temp_present_SSP2_SSP5_50_70_density)

ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/ann_mean_temp_present_SSP2_SSP5_50_70_density.png", width = 40, height = 20, units = "cm")

#mean temp per species
mean_temp <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_temp_species_present.csv")

env_species_plot1 <-ggplot(mean_temp,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#57bfff")+geom_point()+
  geom_hline(aes(yintercept = mean(present_richness_env$mean_ann_mean_temp)),color="blue")+
   coord_cartesian(ylim = c(90,320))+
  labs(y="mean temperature (C*10)", x="Species")+
  ggtitle("Mean Temperature of species occurences- present")
#x11();plot(env_species_plot1)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/temp_species_present.png", width = 20, height = 20, units = "cm")


mean_temp_SSP2_50 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_temp_species_SSP2_50.csv")

env_species_plot2 <-ggplot(mean_temp_SSP2_50,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#57bfff")+geom_point()+
  geom_hline(aes(yintercept = mean(SSP2_2050_richness_env$mean_ann_mean_temp)),color="blue")+
  coord_cartesian(ylim = c(90,320))+
  labs(y="mean temperature (C*10)", x="Species")+
  ggtitle("Mean Temperature of species occurences-SSP2_2050")
#x11();plot(env_species_plot2)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/temp_species_SSP2_50.png", width = 20, height = 20, units = "cm")

mean_temp_SSP5_50 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_temp_species_SSP5_50.csv")

env_species_plot3 <-ggplot(mean_temp_SSP5_50,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#57bfff")+geom_point()+
  geom_hline(aes(yintercept = mean(SSP5_2050_richness_env$mean_ann_mean_temp)),color="blue")+
  coord_cartesian(ylim = c(90,320))+
  labs(y="mean temperature (C*10)", x="Species")+
  ggtitle("Mean Temperature of species occurences-SSP5 2050")
#x11();plot(env_species_plot3)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/temp_species_SSP5_50.png", width = 20, height = 20, units = "cm")


mean_temp_SSP2_70 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_temp_species_SSP2_70.csv")

env_species_plot4 <-ggplot(mean_temp_SSP2_70,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#57bfff")+geom_point()+
  geom_hline(aes(yintercept = mean(SSP2_2070_richness_env$mean_ann_mean_temp)),color="blue")+
  coord_cartesian(ylim = c(90,320))+
  labs(y="mean temperature (C*10)", x="Species")+
  ggtitle("Mean Temperature of species occurences-SSP2 2070")
#x11();plot(env_species_plot4)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/temp_species_SSP2_70.png", width = 20, height = 20, units = "cm")

mean_temp_SSP5_70 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_temp_species_SSP5_70.csv")

env_species_plot5 <-ggplot(mean_temp_SSP5_70,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#57bfff")+geom_point()+
  geom_hline(aes(yintercept = mean(SSP5_2070_richness_env$mean_ann_mean_temp)),color="blue")+
  coord_cartesian(ylim = c(90,320))+
  labs(y="mean temperature (C*10)", x="Species")+
  ggtitle("Mean Temperature of species occurences-SSP5 2070")
x11();plot(env_species_plot5)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/temp_species_SSP5_70.png", width = 20, height = 20, units = "cm")


# change basins according to temperature
basin_sum_present <-data.frame(species=names(present_binary_dist_all), basin_sum= colSums(present_binary_dist_all))

SSP2_50_sum_present <-data.frame(species=names(SSP2_2050_binary_dist_all), basin_sum_SSP2_50= colSums(SSP2_2050_binary_dist_all))
SSP5_50_sum_present <-data.frame(species=names(SSP5_2050_binary_dist_all), basin_sum_SSP5_50= colSums(SSP5_2050_binary_dist_all))

SSP2_70_sum_present <-data.frame(species=names(SSP2_2070_binary_dist_all), basin_sum_SSP2_70= colSums(SSP2_2070_binary_dist_all))
SSP5_70_sum_present <-data.frame(species=names(SSP5_2070_binary_dist_all), basin_sum_SSP5_70= colSums(SSP5_2070_binary_dist_all))


mean_temp_basins <- Reduce(merge, list(mean_temp, basin_sum_present, SSP2_50_sum_present, SSP5_50_sum_present, SSP2_70_sum_present, SSP5_70_sum_present))
mean_temp_basins <- mean_temp_basins[,-3]
mean_temp_basins$change_50_SSP2 <- (mean_temp_basins$basin_sum_SSP2_50/mean_temp_basins$basin_sum)
mean_temp_basins$change_50_SSP5 <- (mean_temp_basins$basin_sum_SSP5_50/mean_temp_basins$basin_sum)

mean_temp_basins$change_70_SSP2 <- (mean_temp_basins$basin_sum_SSP2_70/mean_temp_basins$basin_sum)
mean_temp_basins$change_70_SSP5 <- (mean_temp_basins$basin_sum_SSP5_70/mean_temp_basins$basin_sum)

write.csv(mean_temp_basins,"/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/basins_temp_change.csv")

basins_temp_changes <- ggplot(mean_temp_basins, aes(x=reorder(species,mean), change_50_SSP2))+geom_point()
x11();plot(basins_temp_changes)
                              
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/basins_temp_change_SSP2_50.png", width = 30, height = 10, units = "cm")




# Elevation 

env_present_not_used <- read.csv("/mnt/brunner/brunner/SDM_data/tables/all_var/topo_bio_lc_stream_present.csv")
#ELEVATION
present_richness_elev <-merge(present_richness_table, env_present_not_used[,c("basin_id","mean_elevation_global")], by = "basin_id")
SSP2_2050_richness_elev <-merge(SSP2_2050_richness_table, env_present_not_used[,c("basin_id","mean_elevation_global")], by = "basin_id")
SSP2_2070_richness_elev <-merge(SSP2_2070_richness_table, env_present_not_used[,c("basin_id","mean_elevation_global")], by = "basin_id")
SSP5_2050_richness_elev <-merge(SSP5_2050_richness_table, env_present_not_used[,c("basin_id","mean_elevation_global")], by = "basin_id")
SSP5_2070_richness_elev <-merge(SSP5_2070_richness_table, env_present_not_used[,c("basin_id","mean_elevation_global")], by = "basin_id")

elev_test <- subset(SSP2_2070_richness_elev, richness_SSP2_2070 > 0)
min(elev_test)

mean_elevation_present_SSP2_50_70_density <-ggplot()+
  geom_density(present_richness_elev,  mapping= aes(x=mean_elevation_global, color="#004242"))+
  geom_density(SSP2_2050_richness_elev, mapping=aes(x=mean_elevation_global, color="blue"))+
  geom_density(SSP5_2050_richness_elev, mapping=aes(x=mean_elevation_global, color="green"))+
  geom_density(SSP2_2070_richness_elev , mapping=aes(x=mean_elevation_global, color="#00c2c2"))+
  geom_density(SSP2_2070_richness_elev , mapping=aes(x=mean_elevation_global, color="#0f5b2f"))+
  scale_color_identity(name="projection",
                       breaks = c("#004242", "blue", "green", "#00c2c2", "#0f5b2f"),
                       labels = c("present", "SSP2_2050", "SSP5_2050" ,"SSP2_2070", "SSP5_2070"),
                       guide = "legend")+
  labs(y="density", x=" mean_elevation")+     
  ggtitle("mean_elevation density")+
  x11();plot(mean_elevation_present_SSP2_50_70_density)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_elevation_present_SSP2_50_70_density.png", width = 20, height = 20, units = "cm")

#mean temp per species
mean_elev <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_elev_species_present.csv")

elev_species_plot1 <-ggplot(mean_elev,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#57bfff")+geom_point()+
  geom_hline(aes(yintercept = mean(present_richness_elev$mean_elevation_global)),color="blue")+
  coord_cartesian(ylim = c(0,4000))+
  labs(y="mean elevation (m)", x="Species")+
  ggtitle("Mean elevation of species occurences- present")
#x11();plot(elev_species_plot1)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/elev_species_present.png", width = 20, height = 20, units = "cm")


mean_elev_SSP2_50 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_elev_species_SSP2_50.csv")

elev_species_plot2 <-ggplot(mean_elev_SSP2_50,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#57bfff")+geom_point()+
  geom_hline(aes(yintercept = mean(SSP2_2050_richness_elev$mean_elevation_global)),color="blue")+
  coord_cartesian(ylim = c(0,4000))+
  labs(y="mean elevation (m)", x="Species")+
  ggtitle("Mean elevation of species occurences-SSP2_2050")
#x11();plot(elev_species_plot2)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/elev_species_SSP2_50.png", width = 20, height = 20, units = "cm")

mean_elev_SSP5_50 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_elev_species_SSP5_50.csv")

elev_species_plot3 <-ggplot(mean_elev_SSP5_50,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#57bfff")+geom_point()+
  geom_hline(aes(yintercept = mean(SSP5_2050_richness_elev$mean_elevation_global)),color="blue")+
  coord_cartesian(ylim = c(0,4000))+
  labs(y="mean elevation (m)", x="Species")+
  ggtitle("Mean elevation of species occurences-SSP5_2050")
#x11();plot(elev_species_plot3)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/elev_species_SSP5_50.png", width = 20, height = 20, units = "cm")


mean_elev_SSP2_70 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_elev_species_SSP2_70.csv")

elev_species_plot4 <-ggplot(mean_elev_SSP2_70,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#57bfff")+geom_point()+
  geom_hline(aes(yintercept = mean(SSP2_2070_richness_elev$mean_elevation_global)),color="blue")+
  coord_cartesian(ylim = c(0,4000))+
  labs(y="mean elevation (m)", x="Species")+
  ggtitle("Mean elevation of species occurences-SSP2_2070")
#x11();plot(elev_species_plot4)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/elev_species_SSP2_70.png", width = 20, height = 20, units = "cm")

mean_elev_SSP5_70 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_elev_species_SSP5_70.csv")

elev_species_plot5 <-ggplot(mean_elev_SSP5_70,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#57bfff")+geom_point()+
  geom_hline(aes(yintercept = mean(SSP5_2070_richness_elev$mean_elevation_global)),color="blue")+
  coord_cartesian(ylim = c(0,4000))+
  labs(y="mean elevation (m)", x="Species")+
  ggtitle("Mean elevation of species occurences-SSP5_2070")
x11();plot(elev_species_plot5)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/elev_species_SSP5_70.png", width = 20, height = 20, units = "cm")




#strahler
median_strahler <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/median_sd_strahler_species_present.csv")
nr_median_present <- median_strahler %>% group_by(median) %>% tally()
write.csv(nr_median_present, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/nr_median_present.csv", row.names = F)
median_strahler_SSP2_50 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/median_sd_strahler_species_SSP2_50.csv")
nr_median_SSP2_50 <- median_strahler_SSP2_50 %>% group_by(median) %>% tally()
write.csv(nr_median_SSP2_50, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/nr_median_SSP2_50.csv", row.names = F)
median_strahler_SSP5_50 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/median_sd_strahler_species_SSP5_50.csv")
nr_median_SSP5_50 <- median_strahler_SSP5_50 %>% group_by(median) %>% tally()
write.csv(nr_median_SSP5_50, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/nr_median_SSP5_50.csv", row.names = F)
median_strahler_SSP2_70 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/median_sd_strahler_species_SSP2_70.csv")
nr_median_SSP2_70 <- median_strahler_SSP2_70 %>% group_by(median) %>% tally()
write.csv(nr_median_SSP2_70, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/nr_median_SSP2_70.csv", row.names = F)
median_strahler_SSP5_70 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/median_sd_strahler_species_SSP5_70.csv")
nr_median_SSP5_70 <- median_strahler_SSP5_70 %>% group_by(median) %>% tally()
write.csv(nr_median_SSP5_70, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/nr_median_SSP5_70.csv", row.names = F)

#forest
forest_present_SSP2_50_70_density <-ggplot()+
  geom_density(present_richness_env,  mapping= aes(x=forest, color="#004242"))+
  geom_density(SSP2_2050_richness_env, mapping=aes(x=forest, color="blue"))+
  geom_density(SSP5_2050_richness_env, mapping=aes(x=forest, color="green"))+
  geom_density(SSP2_2070_richness_env , mapping=aes(x=forest, color="#00c2c2"))+
  geom_density(SSP2_2070_richness_env , mapping=aes(x=forest, color="#0f5b2f"))+
  scale_color_identity(name="projection",
                       breaks = c("#004242", "blue", "green", "#00c2c2", "#0f5b2f"),
                       labels = c("present", "SSP2_2050", "SSP5_2050" ,"SSP2_2070", "SSP5_2070"),
                       guide = "legend")+
  labs(y="density", x="percentage of forest")+     
  ggtitle("forest density")+
  #x11();plot(forest_present_SSP2_50_70_density)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/forest_present_SSP2_50_70_density.png", width = 20, height = 20, units = "cm")

mean_forest <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_forest_species_present.csv")

forest_species_plot1 <-ggplot(mean_forest,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#099a1a")+geom_point()+
  geom_hline(aes(yintercept = mean(present_richness_env$forest)),color="blue")+
  scale_y_continuous(limits = c(-0.15, 1.1), breaks = seq(0, 1, by = 0.1))+
  labs(y="forest cover (%)", x="Species")+
  ggtitle("Forest cover of species occurences- present")
#x11();plot(forest_species_plot1)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/forest_species_present.png", width = 20, height = 20, units = "cm")


mean_forest_SSP2_50 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_forest_species_SSP2_50.csv")

forest_species_plot2 <-ggplot(mean_forest_SSP2_50,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#099a1a")+geom_point()+
  geom_hline(aes(yintercept = mean(SSP2_2050_richness_env$forest)),color="blue")+
  scale_y_continuous(limits = c(-0.15, 1.1), breaks = seq(0, 1, by = 0.1))+
  labs(y="forest cover (%)", x="Species")+
  ggtitle("forest cover of species occurences-SSP2_2050")
#x11();plot(forest_species_plot2)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/forest_species_SSP2_50.png", width = 20, height = 20, units = "cm")

mean_forest_SSP5_50 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_forest_species_SSP5_50.csv")

forest_species_plot3 <-ggplot(mean_forest_SSP5_50,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#099a1a")+geom_point()+
  geom_hline(aes(yintercept = mean(SSP5_2050_richness_env$forest)),color="blue")+
  scale_y_continuous(limits = c(-0.15, 1.1), breaks = seq(0, 1, by = 0.1))+
  labs(y="forest cover (%)", x="Species")+
  ggtitle("forest cover of species occurences-SSP5 2050")
#x11();plot(forest_species_plot3)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/forest_species_SSP5_50.png", width = 20, height = 20, units = "cm")


mean_forest_SSP2_70 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_forest_species_SSP2_70.csv")

forest_species_plot4 <-ggplot(mean_forest_SSP2_70,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#099a1a")+geom_point()+
  geom_hline(aes(yintercept = mean(SSP2_2070_richness_env$forest)),color="blue")+
  scale_y_continuous(limits = c(-0.15, 1.1), breaks = seq(0, 1, by = 0.1))+
  labs(y="forest cover (%)", x="Species")+
  ggtitle("forest cover of species occurences-SSP2 2070")
#x11();plot(forest_species_plot4)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/forest_species_SSP2_70.png", width = 20, height = 20, units = "cm")

mean_forest_SSP5_70 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_forest_species_SSP5_70.csv")

forest_species_plot5 <-ggplot(mean_forest_SSP5_70,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#099a1a")+geom_point()+
  geom_hline(aes(yintercept = mean(SSP5_2070_richness_env$forest)),color="blue")+
  scale_y_continuous(limits = c(-0.15, 1.1), breaks = seq(0, 1, by = 0.1))+
  labs(y="forest cover (%)", x="Species")+
  ggtitle("forest cover of species occurences-SSP5 2070")
#x11();plot(forest_species_plot5)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/forest_species_SSP5_70.png", width = 20, height = 20, units = "cm")


# change basins according to forest
basin_sum_present <-data.frame(species=names(present_binary_dist_all), basin_sum= colSums(present_binary_dist_all))

SSP2_50_sum_present <-data.frame(species=names(SSP2_2050_binary_dist_all), basin_sum_SSP2_50= colSums(SSP2_2050_binary_dist_all))
SSP5_50_sum_present <-data.frame(species=names(SSP5_2050_binary_dist_all), basin_sum_SSP5_50= colSums(SSP5_2050_binary_dist_all))

SSP2_70_sum_present <-data.frame(species=names(SSP2_2070_binary_dist_all), basin_sum_SSP2_70= colSums(SSP2_2070_binary_dist_all))
SSP5_70_sum_present <-data.frame(species=names(SSP5_2070_binary_dist_all), basin_sum_SSP5_70= colSums(SSP5_2070_binary_dist_all))


mean_forest_basins <- Reduce(merge, list(mean_forest, basin_sum_present, SSP2_50_sum_present, SSP5_50_sum_present, SSP2_70_sum_present, SSP5_70_sum_present))
mean_forest_basins <- mean_forest_basins[,-3]
mean_forest_basins$change_50_SSP2 <- (mean_forest_basins$basin_sum_SSP2_50/mean_forest_basins$basin_sum)
mean_forest_basins$change_50_SSP5 <- (mean_forest_basins$basin_sum_SSP5_50/mean_forest_basins$basin_sum)

mean_forest_basins$change_70_SSP2 <- (mean_forest_basins$basin_sum_SSP2_70/mean_forest_basins$basin_sum)
mean_forest_basins$change_70_SSP5 <- (mean_forest_basins$basin_sum_SSP5_70/mean_forest_basins$basin_sum)

write.csv(mean_forest_basins,"/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/basins_forest_change.csv")

basins_forest_changes <- ggplot(mean_forest_basins, aes(x=reorder(species,mean), change_50_SSP2))+geom_point()
x11();plot(basins_forest_changes)

ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/basins_forest_change_SSP2_50.png", width = 30, height = 10, units = "cm")

 #grouping
mean_forest$group<-  cut(mean_forest$mean, breaks= seq(0,1, by=0.2), right=T)
nr_forest_groups <- mean_forest %>% group_by(group) %>% tally

mean_forest_SSP2_50$group<-  cut(mean_forest_SSP2_50$mean, breaks= seq(0,1, by=0.2), right=T)
nr_forest_groups_SSP2_50 <- mean_forest_SSP2_50 %>% group_by(group) %>% tally()
mean_forest_SSP5_50$group<-  cut(mean_forest_SSP5_50$mean, breaks= seq(0,1, by=0.2), right=T)
nr_forest_groups_SSP5_50 <- mean_forest_SSP5_50 %>% group_by(group) %>% tally

mean_forest_SSP2_70$group<-  cut(mean_forest_SSP2_70$mean, breaks= seq(0,1, by=0.2), right=T)
nr_forest_groups_SSP2_70 <- mean_forest_SSP2_70 %>% group_by(group) %>% tally
mean_forest_SSP5_70$group<-  cut(mean_forest_SSP5_70$mean, breaks= seq(0,1, by=0.2), right=T)
nr_forest_groups_SSP5_70 <- mean_forest_SSP5_70 %>% group_by(group) %>% tally
 
all_groups_forest <- cbind(nr_forest_groups, nr_forest_groups_SSP2_50[1:5,], nr_forest_groups_SSP5_50[1:5,],nr_forest_groups_SSP2_70[1:5,],nr_forest_groups_SSP5_70[1:5,])
write.csv(all_groups_forest,"/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/groups_forest.csv" )

#crops
crops_present_SSP2_50_70_density <-ggplot()+
  geom_density(present_richness_env,  mapping= aes(x=crops, color="#004242"))+
  geom_density(SSP2_2050_richness_env, mapping=aes(x=crops, color="blue"))+
  geom_density(SSP5_2050_richness_env, mapping=aes(x=crops, color="green"))+
  geom_density(SSP2_2070_richness_env , mapping=aes(x=crops, color="#00c2c2"))+
  geom_density(SSP2_2070_richness_env , mapping=aes(x=crops, color="#0f5b2f"))+
  scale_color_identity(name="projection",
                       breaks = c("#004242", "blue", "green", "#00c2c2", "#0f5b2f"),
                       labels = c("present", "SSP2_2050", "SSP5_2050" ,"SSP2_2070", "SSP5_2070"),
                       guide = "legend")+
  labs(y="density", x="percentage of crops")+     
  ggtitle("crops density")+
 # x11();plot(crops_present_SSP2_50_70_density)
ggsave("/mnt/brunner/brunner/Species_results/graphs_bin/crops_present_SSP2_50_70_density.png", width = 20, height = 20, units = "cm")

mean_crops <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_crops_species_present.csv")

crops_species_plot1 <-ggplot(mean_crops,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#dbbf1b")+geom_point()+
  geom_hline(aes(yintercept = mean(present_richness_env$crops)),color="blue")+
  scale_y_continuous(limits = c(-0.15, 1.1), breaks = seq(0, 1, by = 0.1))+
  labs(y="crops cover (%)", x="Species")+
  ggtitle("crops cover of species occurences- present")
#x11();plot(crops_species_plot1)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/crops_species_present.png", width = 20, height = 20, units = "cm")


mean_crops_SSP2_50 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_crops_species_SSP2_50.csv")

crops_species_plot2 <-ggplot(mean_crops_SSP2_50,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#dbbf1b")+geom_point()+
  geom_hline(aes(yintercept = mean(SSP2_2050_richness_env$crops)),color="blue")+
  scale_y_continuous(limits = c(-0.15, 1.1), breaks = seq(0, 1, by = 0.1))+
  labs(y="crops cover (%)", x="Species")+
  ggtitle("crops cover of species occurences-SSP2_2050")
#x11();plot(crops_species_plot2)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/crops_species_SSP2_50.png", width = 20, height = 20, units = "cm")

mean_crops_SSP5_50 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_crops_species_SSP5_50.csv")

crops_species_plot3 <-ggplot(mean_crops_SSP5_50,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#dbbf1b")+geom_point()+
  geom_hline(aes(yintercept = mean(SSP5_2050_richness_env$crops)),color="blue")+
  scale_y_continuous(limits = c(-0.15, 1.1), breaks = seq(0, 1, by = 0.1))+
  labs(y="crops cover (%)", x="Species")+
  ggtitle("crops cover of species occurences-SSP5 2050")
#x11();plot(crops_species_plot3)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/crops_species_SSP5_50.png", width = 20, height = 20, units = "cm")


mean_crops_SSP2_70 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_crops_species_SSP2_70.csv")

crops_species_plot4 <-ggplot(mean_crops_SSP2_70,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#dbbf1b")+geom_point()+
  geom_hline(aes(yintercept = mean(SSP2_2070_richness_env$crops)),color="blue")+
  scale_y_continuous(limits = c(-0.15, 1.1), breaks = seq(0, 1, by = 0.1))+
  labs(y="crops cover (%)", x="Species")+
  ggtitle("crops cover of species occurences-SSP2 2070")
#x11();plot(crops_species_plot4)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/crops_species_SSP2_70.png", width = 20, height = 20, units = "cm")

mean_crops_SSP5_70 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_crops_species_SSP5_70.csv")

crops_species_plot5 <-ggplot(mean_crops_SSP5_70,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#dbbf1b")+geom_point()+
  geom_hline(aes(yintercept = mean(SSP5_2070_richness_env$crops)),color="blue")+
  scale_y_continuous(limits = c(-0.15, 1.1), breaks = seq(0, 1, by = 0.1))+
  labs(y="crops cover (%)", x="Species")+
  ggtitle("crops cover of species occurences-SSP5 2070")
#x11();plot(crops_species_plot5)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/crops_species_SSP5_70.png", width = 20, height = 20, units = "cm")

#grouping
mean_crops$group<-  cut(mean_crops$mean, breaks= seq(0,0.8, by=0.05), right=T)
nr_crops_groups <- mean_crops %>% group_by(group) %>% tally

mean_crops_SSP2_50$group<-  cut(mean_crops_SSP2_50$mean, breaks= seq(0,0.8, by=0.05), right=T)
nr_crops_groups_SSP2_50 <- mean_crops_SSP2_50 %>% group_by(group) %>% tally()
mean_crops_SSP5_50$group<-  cut(mean_crops_SSP5_50$mean, breaks= seq(0,0.8, by=0.05), right=T)
nr_crops_groups_SSP5_50 <- mean_crops_SSP5_50 %>% group_by(group) %>% tally

mean_crops_SSP2_70$group<-  cut(mean_crops_SSP2_70$mean, breaks= seq(0,0.8, by=0.05), right=T)
nr_crops_groups_SSP2_70 <- mean_crops_SSP2_70 %>% group_by(group) %>% tally
mean_crops_SSP5_70$group<-  cut(mean_crops_SSP5_70$mean, breaks= seq(0,0.8, by=0.05), right=T)
nr_crops_groups_SSP5_70 <- mean_crops_SSP5_70 %>% group_by(group) %>% tally

myls <- list(nr_crops_groups, nr_crops_groups_SSP2_50, nr_crops_groups_SSP5_50,nr_crops_groups_SSP2_70,nr_crops_groups_SSP5_70)
#maximum number of rows
max.rows <- max(nrow(nr_crops_groups), nrow(nr_crops_groups_SSP2_50), nrow(nr_crops_groups_SSP5_50), nrow(nr_crops_groups_SSP2_70), nrow(nr_crops_groups_SSP5_70))
#insert the needed `NA`s to each dataframe
new_myls <- lapply(myls, function(x) { x[1:max.rows,] })
test <-new_myls[[2]]
#create  wanted dataframe
all_groups_crops <- do.call(cbind, lapply(new_myls, `[`, "n"))
all_groups_crops$group <- test$group

write.csv(all_groups_crops,"/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/groups_crops.csv" )

#pasture
pasture_present_SSP2_50_70_density <-ggplot()+
  geom_density(present_richness_env,  mapping= aes(x=pasture, color="#004242"))+
  geom_density(SSP2_2050_richness_env, mapping=aes(x=pasture, color="blue"))+
  geom_density(SSP5_2050_richness_env, mapping=aes(x=pasture, color="green"))+
  geom_density(SSP2_2070_richness_env , mapping=aes(x=pasture, color="#00c2c2"))+
  geom_density(SSP2_2070_richness_env , mapping=aes(x=pasture, color="#0f5b2f"))+
  scale_color_identity(name="projection",
                       breaks = c("#004242", "blue", "green", "#00c2c2", "#0f5b2f"),
                       labels = c("present", "SSP2_2050", "SSP5_2050" ,"SSP2_2070", "SSP5_2070"),
                       guide = "legend")+
  labs(y="density", x="percentage of pasture")+     
  ggtitle("pasture density")+
  x11();plot(pasture_present_SSP2_50_70_density)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/pasture_present_SSP2_50_70_density.png", width = 20, height = 20, units = "cm")

mean_pasture <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_pasture_species_present.csv")

pasture_species_plot1 <-ggplot(mean_pasture,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#71e35f")+geom_point()+
  geom_hline(aes(yintercept = mean(present_richness_env$pasture)),color="blue")+
 scale_y_continuous(limits = c(-0.15, 1.1), breaks = seq(0, 1, by = 0.1))+
    labs(y="pasture cover (%)", x="Species")+
  ggtitle("pasture cover of species occurences- present")
x11();plot(pasture_species_plot1)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/pasture_species_present.png", width = 20, height = 20, units = "cm")


mean_pasture_SSP2_50 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_pasture_species_SSP2_50.csv")

pasture_species_plot2 <-ggplot(mean_pasture_SSP2_50,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#71e35f")+geom_point()+
  geom_hline(aes(yintercept = mean(SSP2_2050_richness_env$pasture)),color="blue")+
  scale_y_continuous(limits = c(-0.15, 1.1), breaks = seq(0, 1, by = 0.1))+
  labs(y="pasture cover (%)", x="Species")+
  ggtitle("pasture cover of species occurences-SSP2_2050")
x11();plot(pasture_species_plot2)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/pasture_species_SSP2_50.png", width = 20, height = 20, units = "cm")

mean_pasture_SSP5_50 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_pasture_species_SSP5_50.csv")

pasture_species_plot3 <-ggplot(mean_pasture_SSP5_50,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#71e35f")+geom_point()+
  geom_hline(aes(yintercept = mean(SSP5_2050_richness_env$pasture)),color="blue")+
  scale_y_continuous(limits = c(-0.15, 1.1), breaks = seq(0, 1, by = 0.1))+
  labs(y="pasture cover (%)", x="Species")+
  ggtitle("pasture cover of species occurences-SSP5 2050")
x11();plot(pasture_species_plot3)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/pasture_species_SSP5_50.png", width = 20, height = 20, units = "cm")


mean_pasture_SSP2_70 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_pasture_species_SSP2_70.csv")

pasture_species_plot4 <-ggplot(mean_pasture_SSP2_70,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#71e35f")+geom_point()+
  geom_hline(aes(yintercept = mean(SSP2_2070_richness_env$pasture)),color="blue")+
  scale_y_continuous(limits = c(-0.15, 1.1), breaks = seq(0, 1, by = 0.1))+
  labs(y="pasture cover (%)", x="Species")+
  ggtitle("pasture cover of species occurences-SSP2 2070")
x11();plot(pasture_species_plot4)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/pasture_species_SSP2_70.png", width = 20, height = 20, units = "cm")

mean_pasture_SSP5_70 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_pasture_species_SSP5_70.csv")

pasture_species_plot5 <-ggplot(mean_pasture_SSP5_70,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#71e35f")+geom_point()+
  geom_hline(aes(yintercept = mean(SSP5_2070_richness_env$pasture)),color="blue")+
  scale_y_continuous(limits = c(-0.15, 1.1), breaks = seq(0, 1, by = 0.1))+
  labs(y="pasture cover (%)", x="Species")+
  ggtitle("pasture cover of species occurences-SSP5 2070")
x11();plot(pasture_species_plot5)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/pasture_species_SSP5_70.png", width = 20, height = 20, units = "cm")

#grouping
mean_pasture$group<-  cut(mean_pasture$mean, breaks= seq(0,0.8, by=0.05), right=T)
nr_pasture_groups <- mean_pasture %>% group_by(group) %>% tally

mean_pasture_SSP2_50$group<-  cut(mean_pasture_SSP2_50$mean, breaks= seq(0,0.8, by=0.05), right=T)
nr_pasture_groups_SSP2_50 <- mean_pasture_SSP2_50 %>% group_by(group) %>% tally()
mean_pasture_SSP5_50$group<-  cut(mean_pasture_SSP5_50$mean, breaks= seq(0,0.8, by=0.05), right=T)
nr_pasture_groups_SSP5_50 <- mean_pasture_SSP5_50 %>% group_by(group) %>% tally

mean_pasture_SSP2_70$group<-  cut(mean_pasture_SSP2_70$mean, breaks= seq(0,0.8, by=0.05), right=T)
nr_pasture_groups_SSP2_70 <- mean_pasture_SSP2_70 %>% group_by(group) %>% tally
mean_pasture_SSP5_70$group<-  cut(mean_pasture_SSP5_70$mean, breaks= seq(0,0.8, by=0.05), right=T)
nr_pasture_groups_SSP5_70 <- mean_pasture_SSP5_70 %>% group_by(group) %>% tally

all_groups_pasture <- cbind(nr_pasture_groups, nr_pasture_groups_SSP2_50, nr_pasture_groups_SSP5_50,nr_pasture_groups_SSP2_70,nr_pasture_groups_SSP5_70)

write.csv(all_groups_pasture,"/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/groups_pasture.csv" )

#urban
urban_present_SSP2_50_70_density <-ggplot()+
  geom_density(present_richness_env,  mapping= aes(x=urban, color="#004242"))+
  geom_density(SSP2_2050_richness_env, mapping=aes(x=urban, color="blue"))+
  geom_density(SSP5_2050_richness_env, mapping=aes(x=urban, color="green"))+
  geom_density(SSP2_2070_richness_env , mapping=aes(x=urban, color="#00c2c2"))+
  geom_density(SSP2_2070_richness_env , mapping=aes(x=urban, color="#0f5b2f"))+
  scale_color_identity(name="projection",
                       breaks = c("#004242", "blue", "green", "#00c2c2", "#0f5b2f"),
                       labels = c("present", "SSP2_2050", "SSP5_2050" ,"SSP2_2070", "SSP5_2070"),
                       guide = "legend")+
  labs(y="density", x="percentage of urban")+     
  ggtitle("urban density")+
  x11();plot(urban_present_SSP2_50_70_density)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/urban_present_SSP2_50_70_density.png", width = 20, height = 20, units = "cm")

mean_urban <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_urban_species_present.csv")

urban_species_plot1 <-ggplot(mean_urban,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#5e3636")+geom_point()+
  geom_hline(aes(yintercept = mean(present_richness_env$urban)),color="blue")+
  scale_y_continuous(limits = c(-0.05, 0.2), breaks = seq(0, 0.2, by = 0.05))+
  labs(y="urban cover (%)", x="Species")+
  ggtitle("urban cover of species occurences- present")
x11();plot(urban_species_plot1)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/urban_species_present.png", width = 20, height = 20, units = "cm")


mean_urban_SSP2_50 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_urban_species_SSP2_50.csv")

urban_species_plot2 <-ggplot(mean_urban_SSP2_50,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#5e3636")+geom_point()+
  geom_hline(aes(yintercept = mean(SSP2_2050_richness_env$urban)),color="blue")+
 scale_y_continuous(limits = c(-0.05, 0.2), breaks = seq(0, 0.2, by = 0.05))+
  labs(y="urban cover (%)", x="Species")+
  ggtitle("urban cover of species occurences-SSP2_2050")
x11();plot(urban_species_plot2)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/urban_species_SSP2_50.png", width = 20, height = 20, units = "cm")

mean_urban_SSP5_50 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_urban_species_SSP5_50.csv")

urban_species_plot3 <-ggplot(mean_urban_SSP5_50,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#5e3636")+geom_point()+
  geom_hline(aes(yintercept = mean(SSP5_2050_richness_env$urban)),color="blue")+
 scale_y_continuous(limits = c(-0.05, 0.2), breaks = seq(0, 0.2, by = 0.05))+
  labs(y="urban cover (%)", x="Species")+
  ggtitle("urban cover of species occurences-SSP5 2050")
x11();plot(urban_species_plot3)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/urban_species_SSP5_50.png", width = 20, height = 20, units = "cm")


mean_urban_SSP2_70 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_urban_species_SSP2_70.csv")

urban_species_plot4 <-ggplot(mean_urban_SSP2_70,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#5e3636")+geom_point()+
  geom_hline(aes(yintercept = mean(SSP2_2070_richness_env$urban)),color="blue")+
 scale_y_continuous(limits = c(-0.05, 0.2), breaks = seq(0, 0.2, by = 0.05))+
  labs(y="urban cover (%)", x="Species")+
  ggtitle("urban cover of species occurences-SSP2 2070")
x11();plot(urban_species_plot4)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/urban_species_SSP2_70.png", width = 20, height = 20, units = "cm")

mean_urban_SSP5_70 <- read.csv("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/mean_sd_urban_species_SSP5_70.csv")

urban_species_plot5 <-ggplot(mean_urban_SSP5_70,aes(x=reorder(species,mean),mean))+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2, color="#5e3636")+geom_point()+
  geom_hline(aes(yintercept = mean(SSP5_2070_richness_env$urban)),color="blue")+
 scale_y_continuous(limits = c(-0.05, 0.2), breaks = seq(0, 0.2, by = 0.05))+
  labs(y="urban cover (%)", x="Species")+
  ggtitle("urban cover of species occurences-SSP5 2070")
x11();plot(urban_species_plot5)
ggsave("/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/urban_species_SSP5_70.png", width = 20, height = 20, units = "cm")

#grouping
mean_urban$group<-  cut(mean_urban$mean, breaks= seq(0,0.2, by=0.01), right=T)
nr_urban_groups <- mean_urban %>% group_by(group) %>% tally

mean_urban_SSP2_50$group<-  cut(mean_urban_SSP2_50$mean, breaks= seq(0,0.2, by=0.01), right=T)
nr_urban_groups_SSP2_50 <- mean_urban_SSP2_50 %>% group_by(group) %>% tally()
mean_urban_SSP5_50$group<-  cut(mean_urban_SSP5_50$mean, breaks= seq(0,0.2, by=0.01), right=T)
nr_urban_groups_SSP5_50 <- mean_urban_SSP5_50 %>% group_by(group) %>% tally

mean_urban_SSP2_70$group<-  cut(mean_urban_SSP2_70$mean, breaks= seq(0,0.2, by=0.01), right=T)
nr_urban_groups_SSP2_70 <- mean_urban_SSP2_70 %>% group_by(group) %>% tally
mean_urban_SSP5_70$group<-  cut(mean_urban_SSP5_70$mean, breaks= seq(0,0.2, by=0.01), right=T)
nr_urban_groups_SSP5_70 <- mean_urban_SSP5_70 %>% group_by(group) %>% tally

myls <- list(nr_urban_groups, nr_urban_groups_SSP2_50, nr_urban_groups_SSP5_50,nr_urban_groups_SSP2_70,nr_urban_groups_SSP5_70)
#maximum number of rows
max.rows <- max(nrow(nr_urban_groups), nrow(nr_urban_groups_SSP2_50), nrow(nr_urban_groups_SSP5_50), nrow(nr_urban_groups_SSP2_70), nrow(nr_urban_groups_SSP5_70))
#insert the needed `NA`s to each dataframe
new_myls <- lapply(myls, function(x) { x[1:max.rows,] })
test <-new_myls[[2]]
#create  wanted dataframe
all_groups_urban <- do.call(cbind, lapply(new_myls, `[`, "n"))
all_groups_urban$group <- test$group

write.csv(all_groups_urban,"/mnt/brunner/brunner/Species_results/binary_dist_endemic/env/groups_urban.csv" )

#================================= see where changes are happening
richness_present <- raster("/mnt/brunner/brunner/Species_results/binary_dist_endemic/present_richness_map_bin_dist_endemic.tif")


richness_SSP1_2050 <- raster("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP1_2050_richness_map_bin_dist_endemic.tif")
richness_SSP2_2050 <- raster("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP2_2050_richness_map_bin_dist_endemic.tif")
richness_SSP3_2050 <- raster("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP3_2050_richness_map_bin_dist_endemic.tif")
richness_SSP4_2050 <- raster("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP4_2050_richness_map_bin_dist_endemic.tif")
richness_SSP5_2050 <- raster("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP5_2050_richness_map_bin_dist_endemic.tif")


richness_SSP1_2070 <- raster("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP1_2070_richness_map_bin_dist_endemic.tif")
richness_SSP2_2070 <- raster("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP2_2070_richness_map_bin_dist_endemic.tif")
richness_SSP3_2070 <- raster("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP3_2070_richness_map_bin_dist_endemic.tif")
richness_SSP4_2070 <- raster("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP4_2070_richness_map_bin_dist_endemic.tif")
richness_SSP5_2070 <- raster("/mnt/brunner/brunner/Species_results/binary_dist_endemic/SSP5_2070_richness_map_bin_dist_endemic.tif")




#=== change present to 2050 and 2070
change_SSP1_2050 <- richness_SSP1_2050 - richness_present 
#x11();plot(change_SSP1_2050)
writeRaster(change_SSP1_2050, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/change_maps/change_SSP1_2050_endemic.tif", overwrite =T)

change_SSP1_2070 <- richness_SSP1_2070 - richness_present 
#x11();plot(change_SSP1_2070)
writeRaster(change_SSP1_2070, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/change_maps/change_SSP1_2070_endemic.tif", overwrite =T)

change_SSP2_2050 <-  richness_SSP2_2050 - richness_present 
#x11();plot(change_SSP2_2050)
writeRaster(change_SSP2_2050, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/change_maps/change_SSP2_2050_endemic.tif", overwrite =T)

change_SSP2_2070 <-  richness_SSP2_2070 - richness_present 
#x11();plot(change_SSP2_2070)
writeRaster(change_SSP2_2070, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/change_maps/change_SSP2_2070_endemic.tif", overwrite =T)

change_SSP3_2050 <-  richness_SSP3_2050 - richness_present 
#x11();plot(change_SSP3_2050)
writeRaster(change_SSP3_2050, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/change_maps/change_SSP3_2050_endemic.tif", overwrite =T)

change_SSP3_2070 <-  richness_SSP3_2070 - richness_present 
#x11();plot(change_SSP3_2070)
writeRaster(change_SSP3_2070, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/change_maps/change_SSP3_2070_endemic.tif", overwrite =T)

change_SSP4_2050 <-  richness_SSP4_2050 - richness_present 
#x11();plot(change_SSP4_2070)
writeRaster(change_SSP4_2050, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/change_maps/change_SSP4_2050_endemic.tif", overwrite =T)

change_SSP4_2070 <-  richness_SSP4_2070 - richness_present 
#x11();plot(change_SSP4_2070)
writeRaster(change_SSP4_2070, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/change_maps/change_SSP4_2070_endemic.tif", overwrite =T)

change_SSP5_2050 <- richness_SSP5_2050 - richness_present 
#x11();plot(change_SSP5_2070)
writeRaster(change_SSP5_2050, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/change_maps/change_SSP5_2050_endemic.tif", overwrite =T)


change_SSP5_2070 <-  richness_SSP5_2070 - richness_present 
#x11();plot(change_SSP5_2070)
writeRaster(change_SSP5_2070, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/change_maps/change_SSP5_2050_endemic.tif", overwrite =T)

#=== change between 2050 and 2070
change_SSP1_2050_2070 <-  richness_SSP1_2070 - richness_SSP1_2050
#x11();plot(change_SSP1_2050_2070)
writeRaster(change_SSP1_2050_2070, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/change_maps/change_SSP1_2050_2070_endemic.tif", overwrite =T)

change_SSP2_2050_2070 <-  richness_SSP2_2070 - richness_SSP2_2050
#x11();plot(change_SSP2_2050_2070)
writeRaster(change_SSP2_2050_2070, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/change_maps/change_SSP2_2050_2070_endemic.tif", overwrite =T)

change_SSP3_2050_2070 <-  richness_SSP3_2070 - richness_SSP3_2050
#x11();plot(change_SSP3_2050_2070)
writeRaster(change_SSP3_2050_2070, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/change_maps/change_SSP3_2050_2070_endemic.tif", overwrite =T)

change_SSP4_2050_2070 <-  richness_SSP4_2070 - richness_SSP4_2050
#x11();plot(change_SSP4_2050_2070)
writeRaster(change_SSP4_2050_2070, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/change_maps/change_SSP4_2050_2070_endemic.tif", overwrite =T)

change_SSP5_2050_2070 <-  richness_SSP5_2070 - richness_SSP5_2050
#x11();plot(change_SSP5_2050_2070)
writeRaster(change_SSP5_2050_2070, "/mnt/brunner/brunner/Species_results/binary_dist_endemic/change_maps/change_SSP5_2050_2070_endemic.tif", overwrite =T)



#====differences between SSPs

change_SSP1_SSP2 <-  richness_SSP1_2050 - richness_SSP2_2050
#x11();plot(change_SSP1_SSP2)
writeRaster(change_SSP1_SSP2 , "/mnt/brunner/brunner/Species_results/binary_dist_endemic/change_maps/change_SSP1_SSP2.tif", overwrite =T)

change_SSP1_SSP3 <-  richness_SSP1_2050 - richness_SSP3_2050
#x11();plot(change_SSP1_SSP3)
writeRaster(change_SSP1_SSP3 , "/mnt/brunner/brunner/Species_results/binary_dist_endemic/change_maps/change_SSP1_SSP3.tif", overwrite =T)

change_SSP1_SSP4 <-  richness_SSP1_2050 - richness_SSP4_2050
#x11();plot(change_SSP1_SSP4)
writeRaster(change_SSP1_SSP4 , "/mnt/brunner/brunner/Species_results/binary_dist_endemic/change_maps/change_SSP1_SSP4.tif", overwrite =T)

change_SSP1_SSP5 <-  richness_SSP1_2050 - richness_SSP5_2050
#x11();plot(change_SSP1_SSP5)
writeRaster(change_SSP1_SSP5 , "/mnt/brunner/brunner/Species_results/binary_dist_endemic/change_maps/change_SSP1_SSP5.tif", overwrite =T)

change_SSP5_SSP2_2050 <-  richness_SSP5_2050 - richness_SSP2_2050 
#x11();plot(change_SSP2_SSP5)
writeRaster(change_SSP5_SSP2_2050 , "/mnt/brunner/brunner/Species_results/binary_dist_endemic/change_maps/change_SSP2_SSP5_2050_endemic.tif", overwrite =T)

change_SSP5_SSP2_2070 <-   richness_SSP5_2070 - richness_SSP2_2070 
#x11();plot(change_SSP2_SSP5)
writeRaster(change_SSP5_SSP2_2070 , "/mnt/brunner/brunner/Species_results/binary_dist_endemic/change_maps/change_SSP2_SSP5_2070_endemic.tif", overwrite =T)

change_SSP5_SSP3 <-  richness_SSP5_2050 - richness_SSP3_2050
#x11();plot(change_SSP5_SSP3)
writeRaster(change_SSP5_SSP3 , "/mnt/brunner/brunner/Species_results/binary_dist_endemic/change_maps/change_SSP5_SSP3.tif", overwrite =T)

change_SSP5_SSP4 <-  richness_SSP5_2050 - richness_SSP4_2050
#x11();plot(change_SSP5_SSP4)
writeRaster(change_SSP5_SSP4 , "/mnt/brunner/brunner/Species_results/binary_dist_endemic/change_maps/change_SSP5_SSP4.tif", overwrite =T)

change_SSP5_SSP2 <-  richness_SSP5_2050 - richness_SSP2_2050
#x11();plot(change_SSP5_SSP2)
writeRaster(change_SSP5_SSP2 , "/mnt/brunner/brunner/Species_results/binary_dist_endemic/change_maps/change_SSP5_SSP2.tif", overwrite =T)