```R
if(Sys.info()[["nodename"]]=="vmdomisch03") DIR ="/data/brunner/Masterarbeit/R" 
setwd(DIR)
R_temp_delete <- paste0(DIR, "/R_temp_delete")
rasterOptions(tmpdir=R_temp_delete)

rasterOptions()

library(raster)
library(sp)
library(rgdal)
library(maptools)
library(rgeos)
library("spatialEco")
library(dplyr)
library(devtools)
library(ggord)
library(ggplot2)
library(ggspatial)
library(tidyverse)
library(sf)
#===============================LOADING MAPS===========================================#
#=========================================BASIN

#load raster of basin
basin <- raster("/data/brunner/Masterarbeit/maps/Colombia/basin.tif")
#x11(); plot(basin)
basin_df <- as.data.frame(basin, xy=T, na.rm=T)
basin_info <- basin_df %>% group_by(basin) %>% tally()

#=========================================LANDCOVER

#load lc raster of study are
LC_list <- list.files("/data/brunner/Masterarbeit/maps/Colombia/" , pattern = "LCCO", full.names = T)
LC_stack <- stack(LC_list)


# LC_1992 <- raster("/data/brunner/Masterarbeit/maps/Colombia/LCCO_1992.tif" )
# LC_1993 <- raster("/data/brunner/Masterarbeit/maps/Colombia/LCCO_1993.tif" )
# LC_1994 <- raster("/data/brunner/Masterarbeit/maps/Colombia/LCCO_1994.tif" )
# LC_1995 <- raster("/data/brunner/Masterarbeit/maps/Colombia/LCCO_1995.tif" )
# LC_1996 <- raster("/data/brunner/Masterarbeit/maps/Colombia/LCCO_1996.tif" )
# LC_1997 <- raster("/data/brunner/Masterarbeit/maps/Colombia/LCCO_1997.tif" )
# LC_1998 <- raster("/data/brunner/Masterarbeit/maps/Colombia/LCCO_1998.tif" )
# LC_1999 <- raster("/data/brunner/Masterarbeit/maps/Colombia/LCCO_1999.tif" )
# LC_2000 <- raster("/data/brunner/Masterarbeit/maps/Colombia/LCCO_2000.tif" )
# LC_2001 <- raster("/data/brunner/Masterarbeit/maps/Colombia/LCCO_2001.tif" )
# LC_2002 <- raster("/data/brunner/Masterarbeit/maps/Colombia/LCCO_2002.tif" )
# LC_2003 <- raster("/data/brunner/Masterarbeit/maps/Colombia/LCCO_2003.tif" )
# LC_2004 <- raster("/data/brunner/Masterarbeit/maps/Colombia/LCCO_2004.tif" )
# LC_2005 <- raster("/data/brunner/Masterarbeit/maps/Colombia/LCCO_2005.tif" )
# LC_2009 <- raster("/data/brunner/Masterarbeit/maps/Colombia/LCCO_2009.tif" )
# LC_2010 <- raster("/data/brunner/Masterarbeit/maps/Colombia/LCCO_2010.tif" )
# LC_2014 <- raster("/data/brunner/Masterarbeit/maps/Colombia/LCCO_2014.tif" )
# LC_2018 <- raster("/data/brunner/Masterarbeit/maps/Colombia/LCCO_2018.tif" )

#========================================ECOREGION

eco <- shapefile("/data/brunner/Masterarbeit/maps/ecoregions/Ecoregion_CO.shp")
#x11();plot (eco)

# =========== Elevation

dem <- raster("/data/brunner/R/land_data/elevation_global.tif")
dem <- crop(dem, basin)
dem <- mask(dem, basin)
#x11(); plot(dem)

# =========== Slope 

slope <- raster("/data/brunner/R/land_data/slope_global.tif")
slope <- crop(slope, basin)
slope<- mask(slope, basin)
# x11(); plot(slope)

# =========== Roughness 

rough <- raster("/data/brunner/R/land_data/roughness.tif")
rough <- crop(rough, basin)
rough<- mask(rough, basin)

# =========== bioclimdata
biolist <- list.files("/data/brunner/R/bioclim" , pattern = "*.tif$", full.names = T)
bioclim <- stack(biolist)
bioclim <- crop(bioclim, extent(basin))
bioclim <- mask(bioclim, basin)
# x11();plot(bioclim)
# Bio1 = Annual Mean Temperature
# Bio2 = Mean Diurnal Range
# Bio3 = Isothermality
# Bio4 = Temperature Seasonality
# Bio5 = Max Temperature of Warmest Month
# Bio6 = Min Temperature of Coldest Month
# Bio7 = Temperature Annual Range
# Bio8 = Mean Temperature of Wettest Quarter
# Bio9 = Mean Temperature of Driest Quarter
# Bio10 = Mean Temperature of Warmest Quarter
# Bio11 = Mean Temperature of Coldest Quarter
# Bio12 = Annual Precipitation
# Bio13 = Precipitation of Wettest Month
# Bio14 = Precipitation of Driest Month
# Bio15 = Precipitation Seasonality
# Bio16 = Precipitation of Wettest Quarter
# Bio17 = Precipitation of Driest Quarter
# Bio18 = Precipitation of Warmest Quarter
# Bio19 = Precipitation of Coldest Quarter

#============  flow 1960- 2015

flow <- raster("/data/brunner/R/Abfluss/flow.tif")
flow <- crop(flow, basin)
flow <- mask(flow, basin)


#========================================= ALL FISH
#load spatial points of all fish
fish_all <- shapefile("/data/brunner/Masterarbeit/data/fish/spec_sp/spec_sp.shp")


#x11(); plot(basin)
#points(fish_all, pch = 20)

#===========================ADD info to all fish data
# extract info of basins
fish_all<- raster::extract(basin, fish_all, sp = T)

# extract info of Landcovers
fish_all<- raster::extract(LC_stack, fish_all, sp =T)
# fish_all<- raster::extract(LC_1992, fish_all, sp = T)
# fish_all<- raster::extract(LC_1996, fish_all, sp = T)
# fish_all<- raster::extract(LC_2000, fish_all, sp = T)
# fish_all<- raster::extract(LC_2001, fish_all, sp = T)
# fish_all<- raster::extract(LC_2005, fish_all, sp = T)
# fish_all<- raster::extract(LC_2009, fish_all, sp = T)
# fish_all<- raster::extract(LC_2010, fish_all, sp = T)
# fish_all<- raster::extract(LC_2014, fish_all, sp = T)
# fish_all<- raster::extract(LC_2018, fish_all, sp = T)

# extract info of ecoregions (< because ecoregions are a shp extract is not possible > "over" creates a df with infos)
info_eco <- fish_all %over% eco

# extract info of elevation
fish_all <-raster::extract(dem, fish_all, sp = T)

# extract info of slope
 fish_all <-raster::extract(slope, fish_all, sp = T)

# extract info of roughness
 fish_all <-raster::extract(rough, fish_all, sp = T)
 
 # extract info of bioclim
 fish_all <-raster::extract(bioclim, fish_all, sp = T)
 
#exctract info flow
 fish_all <-raster::extract(flow, fish_all, sp = T)

# make data frame to extract data to table
fish_all_info <- fish_all@data
#make all family names to lowercase letters
fish_all_info$Family <- tolower(fish_all_info$Family)

#add biom info to fish df
fish_all_info$biom <- info_eco$BIOME_NAME
fish_all_info$realm <- info_eco$REALM
fish_all_info$eco_name <- info_eco$ECO_NAME
 write.csv(fish_all_info, "/data/brunner/Masterarbeit/data/fish/fish_all_info_0922.csv")
#===========================check distributions
# Distribution year
years <- fish_all_info %>% group_by(Year) %>% tally()
x11(); barplot(years$n, names.arg = years$Year)
# Distribution families
families <- fish_all_info %>% group_by(Family) %>% tally()
x11(); barplot(families$n, names.arg = families$Family)
write.csv(families, "/data/brunner/Masterarbeit/data/fish/families_all.csv")
# Distribution basins
basins_all<- fish_all_info %>% group_by(basin) %>% tally()

#=========================================FISH FROM 1992
#=========================limit the fish spatialpoints and dataframe beginning 1992
# subset the spatialpoint dataframe beginning with the year 1992-2018
fish_92 <- subset(fish_all, Year > 1991 & Year < 2019)

# make a dataframe frome this cut spatialpoints
fish_92_info <- fish_92@data
#make all family names to lowercase letters
fish_92_info$Family <- tolower(fish_92_info$Family)


# make a cut dataframe from the all fish dataframe
fish_92_df <- subset(fish_all_info, Year > 1991 & Year < 2019)
sample <-fish_92_df
write.csv(fish_92_df, "/data/brunner/Masterarbeit/data/fish/fish_92-18_info.csv", row.names = F)

#=========================check distributions & stats
# Distribution years
years_92 <- fish_92_info %>% group_by(Year) %>% tally()
write.csv(fish_92_df, "/data/brunner/Masterarbeit/data/fish/years_92-18.csv", row.names = F)
x11(); barplot(years_92$n, names.arg = years_92$Year)

# Distribution families
families_92 <- fish_92_info %>% group_by(Family) %>% tally()
#x11(); barplot(families_92$n, names.arg = families_92$Family)
write.csv(families_92, "/data/brunner/Masterarbeit/data/fish/families_92-18.csv", row.names = F)

# Distribution species
species_92 <- fish_92_info %>% group_by(Species) %>% tally()
write.csv(species_92, "/data/brunner/Masterarbeit/data/fish/species_92-18.csv", row.names = F)

# Distribution basins
basins_92 <- fish_92_info %>% group_by(basin) %>% tally()
write.csv(basins_92, "/data/brunner/Masterarbeit/data/fish/basins_92-18.csv", row.names = F)
#distributionEcoregions
bioms_92 <-fish_92_df %>% group_by(biom)%>% tally()
write.csv(bioms_92, "/data/brunner/Masterarbeit/data/fish/biomes_92-18.csv", row.names = F)
realms_92 <-fish_92_df %>% group_by(realm)%>% tally()
ecoregions_92 <-fish_92_df %>% group_by(eco_name)%>% tally()
write.csv(ecoregions_92, "/data/brunner/Masterarbeit/data/fish/econames_92-18.csv", row.names = F)

# add the most common basin for each family
sample_maj <- sample %>%
  # add a column n with count by categories
  add_count(Family, basin) %>%
  # select max or first occurrence by family
  group_by(Family) %>%
  # keep only first TRUE
  mutate(Majority = basin[n == max(n)][1]) %>%
  # do not keep temp var
  select(-n)
# remove the dublicated families
maj_basin <- sample_maj[!duplicated(sample_maj$Family),]
maj_basin <- subset(maj_basin, select = c(Family, Majority,biom))
write.csv(maj_basin, "/data/brunner/Masterarbeit/data/fish/maj_basin_92-18.csv", row.names = F)



# add the most common Landcover per basin
LC_sample_maj <- sample %>%
  # add a column n with count by categories
  add_count(basin, sample[10:36]) %>%
  # select max or first occurrence by basin
  group_by(basin) %>%
  # keep only first TRUE
  mutate(
    Majority92 = LCCO_1992[n == max(n)][1], Majority93 = LCCO_1993[n == max(n)][1], Majority94 = LCCO_1994[n == max(n)][1],
    Majority95 = LCCO_1995[n == max(n)][1], Majority96 = LCCO_1996[n == max(n)][1], Majority97 = LCCO_1997[n == max(n)][1],
    Majority98 = LCCO_1998[n == max(n)][1], Majority99 = LCCO_1999[n == max(n)][1], Majority00 = LCCO_2000[n == max(n)][1],
    Majority01 = LCCO_2001[n == max(n)][1], Majority02 = LCCO_2002[n == max(n)][1], Majority03 = LCCO_2003[n == max(n)][1],
    Majority04 = LCCO_2004[n == max(n)][1], Majority05 = LCCO_2005[n == max(n)][1], Majority06 = LCCO_2006[n == max(n)][1],
    Majority07 = LCCO_2007[n == max(n)][1], Majority08 = LCCO_2008[n == max(n)][1], Majority09 = LCCO_2009[n == max(n)][1],
    Majority10 = LCCO_2010[n == max(n)][1], Majority11 = LCCO_2011[n == max(n)][1], Majority12 = LCCO_2012[n == max(n)][1],
    Majority13 = LCCO_2013[n == max(n)][1], Majority14 = LCCO_2014[n == max(n)][1], Majority15 = LCCO_2015[n == max(n)][1],
    Majority16 = LCCO_2016[n == max(n)][1], Majority17 = LCCO_2010[n == max(n)][1], Majority18 = LCCO_2018[n == max(n)][1]
  ) %>%
  # do not keep temp var
  select(-n)

maj_LC <- LC_sample_maj[!duplicated(LC_sample_maj$basin),]
maj_LC <- maj_LC[c(10 , 64:90)]


#adding then most commen lC of the majority of LCs
modefunc <- function(x){
  tabresult <- tabulate(x)
  themode <- which(tabresult == max(tabresult))
  if(sum(tabresult == max(tabresult))>1) themode <- NA
  return(themode)
}

maj_LC$mode_all <- apply(maj_LC[2:28], 1, modefunc)
write.csv(maj_LC, "/data/brunner/Masterarbeit/data/fish/maj_LC_92-18.csv", row.names = F)
#adding this mode to the dataframe (for plotting)
sample2 <- sample
sample2 <-merge(x = sample2, y = maj_LC[ , c("basin", "mode_all")], by = "basin", all.x=TRUE)

# # add the most common Landcover per basin
# LC_sample_maj <- sample %>%
#   # add a column n with count by categories
#   add_count(basin,LCCO_1992, LCCO_1996, LCCO_2000, LCCO_2001, LCCO_2005, LCCO_2009, LCCO_2010, LCCO_2014, LCCO_2018) %>%
#   # select max or first occurrence by basin
#   group_by(basin) %>%
#   # keep only first TRUE
#   mutate(
#     Majority92 = LCCO_1992[n == max(n)][1], Majority96 = LCCO_1996[n == max(n)][1], Majority00 = LCCO_2000[n == max(n)][1],
#     Majority01 = LCCO_2001[n == max(n)][1], Majority05 = LCCO_2005[n == max(n)][1], Majority09 = LCCO_2009[n == max(n)][1],
#     Majority10 = LCCO_2010[n == max(n)][1], Majority14 = LCCO_2014[n == max(n)][1], Majority18 = LCCO_2018[n == max(n)][1]
#     ) %>%
#   # do not keep temp var
#   select(-n)
# remove the dublicated basins


# distribution of landcovergroups
LC_sample <- sample[11:37]
LC_n_fish_t <-sapply(X = LC_sample,
       FUN = table)
LC_n_fish <- as.data.frame(LC_n_fish_t)
write.csv(LC_n_fish, "/data/brunner/Masterarbeit/data/fish/LC_n_fish_92-18.csv", row.names = F)

# changes between years
LC_92_96 <-  fish_92_df %>% count(LCCO_1992, LCCO_1996)
write.csv(LC_92_96, "/data/brunner/Masterarbeit/data/fish/LC_92_96.csv")
LC_92_01 <-  fish_92_df %>% count(LCCO_1992, LCCO_2001)
write.csv(LC_92_01, "/data/brunner/Masterarbeit/data/fish/LC_92_01.csv")
LC_96_00 <-  fish_92_df %>% count(LCCO_1996, LCCO_2000)
write.csv(LC_96_00, "/data/brunner/Masterarbeit/data/fish/LC_96_00.csv")
LC_96_05 <-  fish_92_df %>% count(LCCO_1996, LCCO_2005)
write.csv(LC_96_05, "/data/brunner/Masterarbeit/data/fish/LC_96_05.csv")
LC_01_05 <-  fish_92_df %>% count(LCCO_2001, LCCO_2005)
write.csv(LC_01_05, "/data/brunner/Masterarbeit/data/fish/LC_01_05.csv")
LC_01_10 <-  fish_92_df %>% count(LCCO_2001, LCCO_2010)
write.csv(LC_01_10, "/data/brunner/Masterarbeit/data/fish/LC_01_10.csv")
LC_05_09 <-  fish_92_df %>% count(LCCO_2005, LCCO_2009)
write.csv(LC_05_09 , "/data/brunner/Masterarbeit/data/fish/LC_05_09.csv")
LC_05_14 <-  fish_92_df %>% count(LCCO_2005, LCCO_2014)
write.csv(LC_05_14, "/data/brunner/Masterarbeit/data/fish/LC_05_14.csv")
LC_09_18 <-  fish_92_df %>% count(LCCO_2009, LCCO_2018)
write.csv(LC_09_18 , "/data/brunner/Masterarbeit/data/fish/LC_09_18.csv")
LC_10_14 <-  fish_92_df %>% count(LCCO_2010, LCCO_2014)
write.csv(LC_10_14, "/data/brunner/Masterarbeit/data/fish/LC_10_14.csv")
LC_92_18 <-  fish_92_df %>% count(LCCO_1992, LCCO_2018)
write.csv(LC_92_18, "/data/brunner/Masterarbeit/data/fish/LC_92_18.csv")


# year ranges for species
years_species <- fish_92_df %>% group_by(Species) %>% summarise(min_year=min(Year), max_year=max(Year))
years_species <-years_species %>% mutate(year=max_year-min_year)
write.csv(years_species, "/data/brunner/Masterarbeit/data/fish/years_species_18.csv")

#find ranges &stats of environmental data for families
slope_range <- sample %>% group_by(Family) %>% summarise(min_slope=min(slope_global), max_slope = max(slope_global))
elevation_range <- sample %>% group_by(Family) %>% summarise(min_elev=min(elevation_global), max_elev = max(elevation_global))
temp_stats <- sample %>% group_by(Family) %>% summarise(min_mean_temp= min(annual_mean_temperature), max_mean_temp=max(annual_mean_temperature),
                                                        mean_max_temp=mean(max_temp_warmest_month), mean_min_temp=mean(min_temp_coldest_month),
                                                        median_max_temp= mean (max_temp_warmest_month), median_min_temp=median(min_temp_coldest_month))
prec_stats <-sample %>% group_by(Family) %>% summarise(min_mean_prec= min(annual_prec), max_mean_prec=max(annual_prec),
                                                       mean_max_prec=mean(prec_wettest_month), mean_min_prec=mean(prec_driest_month),
                                                       median_max_prec= median(prec_wettest_month), median_min_prec=median(prec_driest_month))

#all stats combines
family_stats<-sample %>% group_by(Family) %>% summarise(min_lat=min(Latitud), max_lat=max(Latitud),
                                                        min_long=min(Longitd),max_long=max(Longitd),
                                                        min_slope=min(slope_global), max_slope = max(slope_global), 
                                                        min_elev=min(elevation_global), max_elev = max(elevation_global),
                                                        min_mean_temp= min(annual_mean_temperature, na.rm = T), max_mean_temp=max(annual_mean_temperature, na.rm = T),
                                                        mean_max_temp=mean(max_temp_warmest_month, na.rm = T), mean_min_temp=mean(min_temp_coldest_month, na.rm = T),
                                                        median_max_temp= mean (max_temp_warmest_month, na.rm = T), median_min_temp=median(min_temp_coldest_month, na.rm = T),
                                                        min_mean_prec= min(annual_prec, na.rm = T), max_mean_prec=max(annual_prec, na.rm = T),
                                                        mean_max_prec=mean(prec_wettest_month, na.rm = T), mean_min_prec=mean(prec_driest_month, na.rm = T),
                                                        median_max_prec= median(prec_wettest_month, na.rm = T), median_min_prec=median(prec_driest_month, na.rm = T))
write.csv(family_stats, "/data/brunner/Masterarbeit/data/fish/92-18_family_stats.csv", row.names = F)

#all stats combines
basin_stats<-sample %>% group_by(basin) %>% summarise(min_lat=min(Latitud), max_lat=max(Latitud),
                                                        min_long=min(Longitd),max_long=max(Longitd),
                                                        min_slope=min(slope_global), max_slope = max(slope_global), 
                                                        min_elev=min(elevation_global), max_elev = max(elevation_global),
                                                        min_mean_temp= min(annual_mean_temperature, na.rm = T), max_mean_temp=max(annual_mean_temperature, na.rm = T),
                                                        mean_max_temp=mean(max_temp_warmest_month, na.rm = T), mean_min_temp=mean(min_temp_coldest_month, na.rm = T),
                                                        median_max_temp= mean (max_temp_warmest_month, na.rm = T), median_min_temp=median(min_temp_coldest_month, na.rm = T),
                                                        min_mean_prec= min(annual_prec, na.rm = T), max_mean_prec=max(annual_prec, na.rm = T),
                                                        mean_max_prec=mean(prec_wettest_month, na.rm = T), mean_min_prec=mean(prec_driest_month, na.rm = T),
                                                        median_max_prec= median(prec_wettest_month, na.rm = T), median_min_prec=median(prec_driest_month, na.rm = T))

write.csv(basin_stats, "/data/brunner/Masterarbeit/data/fish/92-18_basin_stats.csv", row.names = F)
#fishplots
x11(width = 12, height = 9); plot (basin)
points(fish_92, pch =20)

p1 <- ggplot(basin_df, aes(x = x, y = y)) +
  geom_point(data =fish_92_info, aes(x= Longitd,y= Latitud, colour = factor(Family)))
x11(); plot (p1)

p2 <- ggplot(basin_df, aes(x = x, y = y)) +
  geom_point(data =fish_92_info, aes(x= Longitd,y= Latitud, colour = factor(Year)))
x11(width = 12, height = 9 ); plot (p2)

p3 <- ggplot(basin_df, aes(x = x, y = y)) +
  geom_polygon(aes(fill = factor(basin) , group = basin),show.legend = F) +
  geom_point(data =fish_92_df, aes(x= Longitd,y= Latitud, colour = factor(Family)))
x11(width = 15, height = 10);p3

p4 <- ggplot(basin_df, aes(x = x, y = y)) +
  geom_polygon(aes(fill = factor(basin) , group = basin),show.legend = F) +
  geom_point(data =sample2, aes(x= Longitd,y= Latitud, colour = factor(mode_all)))
x11(width = 15, height = 10);p4 + scale_colour_brewer(palette = "Set1")

#aggregation of fish stats per basin
agg_92_median <- aggregate(fish_92_df[12:42], by= list (fish_92_df$basin), FUN = median)

#aggregation per subcatchment
fish_92_info <-fread(paste0("/data/brunner/Masterarbeit/data/fish/fish_92_info", ".csv"), h=T)
diversity_92 <- subset(fish_92_info, select = c(basin, Species))
diversity_92 <- as.data.frame(diversity_92)

# number of different Species per basin
richness <- diversity_92 %>%
            group_by(basin) %>%
            summarise(n_distinct(Species))
names(richness) <- c("zone", "richness")
write.csv(richness,  "/data/brunner/Masterarbeit/data/fish/basin_richness_92-18.csv")

# example for zonal <- to get all basin ids (also where there are no species occurences)
dem_zonal <- zonal(dem, basin, "count", na.rm = T)
dem_zonal <- as.data.frame(dem_zonal)
# merge them and keep all data <- basins with no occurences
div_basin <-merge(richness, dem_zonal, by = "zone", all=TRUE)
div_basin <- div_basin[-3]
# make NAs zero &  # change headers, see ?reclassify
div_basin[is.na(div_basin)]<- 0
names(div_basin) <- c("is", "becomes")
basin_richness <- reclassify(basin, div_basin)
x11();plot (basin_richness)
writeRaster(basin_richness, "/data/brunner/Masterarbeit/maps/fish/richness_per_basin.tif")


# number of species per basin
basin92 <-(fish_92_info %>% group_by(basin) %>% tally())
n_fish_basin <- basins_92 
names(n_fish_basin) <- c("zone", "n") 
write.csv(richness,  "/data/brunner/Masterarbeit/data/fish/n_basin_92-18.csv")
fish_per_basin <-merge(n_fish_basin, dem_zonal, by = "zone", all=TRUE)
fish_per_basin <- fish_per_basin[-3]
fish_per_basin[is.na(fish_per_basin)] <- 0
basin_n_species <- reclassify(basin, fish_per_basin)
x11(); plot(basin_n_species)
writeRaster(basin_n_species, "/data/brunner/Masterarbeit/maps/fish/number_per_basin.tif")

```

