
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
if (!require("svglite")) { install.packages("svglite", dependencies = TRUE) ; library(svglite)} 

endemic <- read.table("E:/Uni/MA Arbeit/Species_results/binary_dist_endemic/endemic_species_list.txt")

#temperature
mean_temp <- read.csv("E:/Uni/MA Arbeit/Species_results/binary_dist_CO/env/mean_sd_temp_species_present.csv")
mean_temp$mean_C <- mean_temp$mean/10
mean_temp$sd_C <- mean_temp$sd/10
endemic_temp2 <- mean_temp[mean_temp$species %in% endemic$V1,]


env_species_temp <-ggplot()+
  geom_errorbar(mean_temp,mapping=aes(x=reorder(species,mean_C), ymin=mean_C-sd_C,ymax=mean_C+sd_C),width=0.2, color="#57bfff")+
  geom_hline(mean_temp, mapping=aes(yintercept = mean(mean_C)),color="blue")+
  geom_point(mean_temp,mapping =aes(x=reorder(species,mean_C),mean_C))+  
  geom_point(endemic_temp2,mapping=aes(x=reorder(species,mean_C),mean_C), color="#d0d0d0")+  
  scale_y_continuous(limits = c(9, 32), breaks = seq(0, 30, by = 5))+
  labs(y="°C", x="Species sorted by mean present temperature (ascending) ")+
  x11();plot(env_species_temp)
ggsave(file="G:/Uni/Publikation/figures/env_species_temp_2.svg", plot=env_species_temp, width=10, height=5)


#precipitation
mean_prec <- read.csv("E:/Uni/MA Arbeit/Species_results/binary_dist_CO/env/mean_sd_prec_species_present.csv")
endemic_prec2 <- mean_prec[mean_prec$species %in% endemic$V1,]

env_species_prec <-ggplot()+
  geom_errorbar(mean_prec,mapping=aes(x=reorder(species,mean), ymin=mean-sd,ymax=mean+sd),width=0.2, color="#2c6b92")+
  geom_hline(mean_prec, mapping=aes(yintercept = mean(mean)),color="blue")+
  geom_point(mean_prec,mapping =aes(x=reorder(species,mean),mean))+  
  geom_point(endemic_prec2,mapping=aes(x=reorder(species,mean),mean), color="#d0d0d0")+  
  scale_y_continuous(limits = c(500, 8000), breaks = seq(500, 8000, by = 500))+
  labs(y="mm", x="Species sorted by mean present precipitation (ascending)")+
  x11();plot(env_species_prec)
ggsave(file="G:/Uni/Publikation/figures/env_species_prec2.svg", plot=env_species_prec, width=10, height=5)

#forest
mean_forest <- read.csv("E:/Uni/MA Arbeit/Species_results/binary_dist_CO/env/mean_sd_forest_species_present.csv")
mean_forest$mean_pcnt <- mean_forest$mean*100
mean_forest$sd_pcnt <-mean_forest$sd*100
endemic_forest2 <- mean_forest[mean_forest$species %in% endemic$V1,]

env_species_forest <-ggplot()+
  geom_errorbar(mean_forest,mapping=aes(x=reorder(species,mean_pcnt), ymin=mean_pcnt-sd_pcnt,ymax=mean_pcnt+sd_pcnt),width=0.2, color="#099a1a")+
  geom_hline(mean_forest, mapping=aes(yintercept = mean(mean_pcnt)),color="blue")+
  geom_point(mean_forest,mapping =aes(x=reorder(species,mean_pcnt),mean_pcnt))+  
  geom_point(endemic_forest2,mapping=aes(x=reorder(species,mean_pcnt),mean_pcnt), color="#d0d0d0")+  
  scale_y_continuous(limits = c(-15, 110), breaks = seq(-10, 100, by = 10))+
  labs(y="forest cover (%)", x="Species sorted by mean present forest cover (ascending)")+
  x11();plot(env_species_forest)
ggsave(file="G:/Uni/Publikation/figures/env_species_forest2.svg", plot=env_species_forest, width=10, height=10)

#crops 
mean_crops <- read.csv("E:/Uni/MA Arbeit/Species_results/binary_dist_CO/env/mean_sd_crops_species_present.csv")
mean_crops$mean_pcnt <- mean_crops$mean*100
mean_crops$sd_pcnt <-mean_crops$sd*100
endemic_crops2 <- mean_crops[mean_crops$species %in% endemic$V1,]

env_species_crops <-ggplot()+
  geom_errorbar(mean_crops,mapping=aes(x=reorder(species,mean_pcnt), ymin=mean_pcnt-sd_pcnt,ymax=mean_pcnt+sd_pcnt),width=0.2, color="#dbbf1b")+
  geom_hline(mean_crops, mapping=aes(yintercept = mean(mean_pcnt)),color="blue")+
  geom_point(mean_crops,mapping =aes(x=reorder(species,mean_pcnt),mean_pcnt))+  
  geom_point(endemic_crops2,mapping=aes(x=reorder(species,mean_pcnt),mean_pcnt), color="#d0d0d0")+  
  scale_y_continuous(limits = c(-15, 110), breaks = seq(-10, 100, by = 10))+
  labs(y="crops cover (%)", x="Species sorted by mean present crops cover (ascending)")+
  x11();plot(env_species_crops)
ggsave(file="G:/Uni/Publikation/figures/env_species_crops2.svg", plot=env_species_crops, width=10, height=10)

#pasture "#71e35f"
mean_pasture <- read.csv("E:/Uni/MA Arbeit/Species_results/binary_dist_CO/env/mean_sd_pasture_species_present.csv")
mean_pasture$mean_pcnt <- mean_pasture$mean*100
mean_pasture$sd_pcnt <-mean_pasture$sd*100
endemic_pasture2 <- mean_pasture[mean_pasture$species %in% endemic$V1,]


env_species_pasture <-ggplot()+
  geom_errorbar(mean_pasture,mapping=aes(x=reorder(species,mean_pcnt), ymin=mean_pcnt-sd_pcnt,ymax=mean_pcnt+sd_pcnt),width=0.2, color="#71e35f")+
  geom_hline(mean_pasture, mapping=aes(yintercept = mean(mean_pcnt)),color="blue")+
  geom_point(mean_pasture,mapping =aes(x=reorder(species,mean_pcnt),mean_pcnt))+  
  geom_point(endemic_pasture2,mapping=aes(x=reorder(species,mean_pcnt),mean_pcnt), color="#d0d0d0")+  
  scale_y_continuous(limits = c(-20, 110), breaks = seq(-10, 100, by = 10))+
  labs(y="pasture cover (%)", x="Species sorted by mean present pasture cover (ascending)")+
  x11();plot(env_species_pasture)
ggsave(file="G:/Uni/Publikation/figures/env_species_pasture2.svg", plot=env_species_pasture, width=10, height=10)

#urban
mean_urban <- read.csv("E:/Uni/MA Arbeit/Species_results/binary_dist_CO/env/mean_sd_urban_species_present.csv")
mean_urban$mean_pcnt <- mean_urban$mean*100
mean_urban$sd_pcnt <-mean_urban$sd*100
endemic_urban2 <- mean_urban[mean_urban$species %in% endemic$V1,]

env_species_urban <-ggplot()+
  geom_errorbar(mean_urban,mapping=aes(x=reorder(species,mean_pcnt), ymin=mean_pcnt-sd_pcnt,ymax=mean_pcnt+sd_pcnt),width=0.2, color="#5e3636")+
  geom_hline(mean_urban, mapping=aes(yintercept = mean(mean_pcnt)),color="blue")+
  geom_point(mean_urban,mapping =aes(x=reorder(species,mean_pcnt),mean_pcnt))+  
  geom_point(endemic_urban2,mapping=aes(x=reorder(species,mean_pcnt),mean_pcnt), color="#d0d0d0")+  
  scale_y_continuous(limits = c(-5, 20), breaks = seq(-5, 20, by = 5))+
  labs(y="urban cover (%)", x="Species sorted by mean present urban cover (ascending)")+
  x11();plot(env_species_urban)
ggsave(file="G:/Uni/Publikation/figures/env_species_urban2.svg", plot=env_species_urban, width=10, height=10)
