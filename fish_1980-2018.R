library(spatialEco)
library(raster)
library(sp)
library(rgdal)
library(maptools)
library(rgeos)
library(dplyr)
library(devtools)
library(ggord)
library(ggplot2)
library(ggspatial)
library(tidyverse)
library(sf)
#========================================= ALL FISH
#load spatial points of all fish
fish_all <- shapefile("/mnt/brunner/brunner/data/fish/spec_sp/spec_sp.shp")

#load basin
basin <- raster("/mnt/brunner/brunner/stream_network/basins_study_area.tif")

#extract basin_id for spatial points
fish_all<- raster::extract(basin, fish_all, sp = T)
colnames(fish_all@data)[1]<- "species"
colnames(fish_all@data)[10] <- "basin_id"

fish_full_info <- fish_all@data
years <- fish_full_info %>% group_by(Year) %>% tally()
p <- ggplot(fish_full_info , aes(x=Year))+
    geom_histogram(binwidth = 1)
x11();plot(p)
ggsave("/mnt/brunner/brunner/Species_results/years_all_occ.png", width = 40, height = 20, units = "cm")
  write.csv(years, "/mnt/brunner/brunner/data/fish/spec_sp/year_info.csv", row.names = F )

#remove spatial points with no basin
fish_all_2 <- sp.na.omit(fish_all, col.name= "basin_id")
writeOGR(obj=fish_all_2, dsn="/mnt/brunner/brunner/data/fish/spec_sp/" , layer="fish_all_2", driver="ESRI Shapefile", overwrite_layer = T)

#make csv from this dataset
fish_all_info <- fish_all_2@data
write.csv(fish_all_info, "/mnt/brunner/brunner/data/fish/spec_sp/fish_all_info.csv", row.names = F  )
years <- fish_all_info %>% group_by(Year) %>% tally()
species <- fish_all_info %>% group_by(species) %>% tally()

#subset spatial datapoints for year of observation 1980-2018
fish_80 <- subset(fish_all_2, Year > 1979 & Year < 2019)
writeOGR(obj=fish_80, dsn="/mnt/brunner/brunner/data/fish/spec_sp/" , layer="fish_80", driver="ESRI Shapefile", overwrite_layer = T)
#fish_80_test <- shapefile("/mnt/brunner/brunner/data/fish/spec_sp/fish_80.shp")

#subset data 
fish_sdm <- fish_80
fish_sdm@data <- fish_80@data[,c("basin_id", "Longitd", "Latitud", "species")]
writeOGR(obj=fish_sdm , dsn="/mnt/brunner/brunner/data/fish/spec_sp/" , layer="fish_sdm", driver="ESRI Shapefile", overwrite_layer = T)


fish_info <- fish_80@data
write.csv(fish_info,"/mnt/brunner/brunner/data/fish/spec_sp/fish_80_info.csv", row.names = F )

species_80 <- fish_info %>% group_by(Species) %>% tally()

x11(); plot(basin)
points(fish_80, pch=20)

x11(); plot(basin)
points(fish_80_test, pch=20)
