###======================================================#
###---- Copy pre-compiled data to working directory -----
###======================================================#
export myDIR=/mnt/$USER/$USER/tmp # set working directory in bash
echo $myDIR # check
mkdir $myDIR # create directory

cd $myDIR  # navigate to that directory
#cp -r /data/shared/spdep_SDMs . # copy to current folder




###=======================================#
###---- Set variables for the server -----
###=======================================#
R
### Get current user
USER <- Sys.getenv("USER") # Linux
if(Sys.info()[["sysname"]]=="Linux") DIR=paste0("/mnt/", USER, "/", USER, "/tmp/SDM_data")
#if(Sys.info()[["sysname"]]=="Windows") DIR="C:/temp/spdep_SDMs" # generic Windows path
# if(Sys.getenv()[["USERNAME"]]=="domisch") DIR="D:/codes/spdep_SDMs" # Sami's path
if(Sys.info()[["sysname"]]=="Linux") N_CORES=5
# if(Sys.info()[["sysname"]]=="Windows") N_CORES=3
#dir.create(DIR)
setwd(DIR)
MAXENT_DIR=paste0(DIR, "/MAXENT_DIR")

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
if (!require("biomod2")) { install.packages("BiocManager") ; library(BiocManager)}
if (!require("biomod2")) { BiocManager::install("biomod2") ; library(biomod2)}
if (!require("viridis")) { install.packages("viridis", dependencies = TRUE) ; library(viridis)}
if (!require("RColorBrewer")) { install.packages("RColorBrewer", dependencies = TRUE) ; library(RColorBrewer)}
if (!require("ggplot2")) { install.packages("ggplot2", dependencies = TRUE) ; library(ggplot2)}
if (!require("ggthemes")) { install.packages("ggthemes", dependencies = TRUE) ; library(ggthemes)} # theme_map()
if (!require("caret")) { install.packages("caret", dependencies = TRUE) ; library(caret)} # correlation
if (!require("scales")) { install.packages("scales", dependencies = TRUE) ; library(scales)} # scaling

### Create temporary raster folder
### Make sure to check and empty once in a while --> use /tmp (default), will be flushed after 10 days
dir.create(paste0(DIR,"/R_temp_delete"))
rasterOptions(tmpdir=paste0(DIR, "/R_temp_delete/"))


#' Transform raster to data.table
#'
#' @param x  Rastereco_behav_ object
#' @param row.names  `NULL` or a character vector giving the row names for the data frame. Missing values are not allowed
#' @param optional  logical. If `TRUE`, setting row names and converting column names (to syntactic names: see make.names) is optional
#' @param xy  logical. If `TRUE`, also return the spatial coordinates
#' @param centroids  logical. If TRUE return the centroids instead of all spatial coordinates (only relevant if xy=TRUE)
#' @param sepNA	logical. If TRUE the parts of the spatial objects are separated by lines that are NA (only if xy=TRUE and, for polygons, if centroids=FALSE
#' @param ...	 Additional arguments (none) passed to `raster::as.data.frame`
#'
#' @value returns a data.table object
#' @examples
#' logo <- brick(system.file("external/rlogo.grd", package="raster"))
#' v <- as.data.table(logo)
#' @import

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

###--------------------------#
###--- Load species data ----
###--------------------------#

### Select the size of the stream network (smaller = faster test)
# THRES=40 # ca. 1500 stream reaches
#THRES=110 # ca. 8700 stream reaches
# THRES=20
# THRES=200


### Load the basin -raster file
basin_id <- raster(paste0(DIR, "/basins_study_area.tif")) # do not use basins_cat here!
# basin_id <- raster(paste0(DIR, "/basin_rwatershed_", THRES, ".tif")) # QTP data
names(basin_id) <- "basin_id"


### Species points
#make species points from csv
# species_info <-fread(paste0(DIR, "/fish_cropped_80", ".csv"), h=T)
# sp_points <- SpatialPointsDataFrame(species_info[,c("Longitd", "Latitud")], data = species_info)
# projection(sp_points) <- projection(basin_id)
# #x11(); plot(sp_points2)
# sp_points$Species <- gsub(" ", ".", sp_points$Species) 


#directly with species points > with cut data set (1980-2018 #)
sp_points <- shapefile(paste0(DIR, "/fish_sdm.shp"))
head(sp_points)
sp_points$species <- gsub(" ", ".", sp_points$species)  # remove whitespace, biomod will remove it later...
 #sp_point_info <- as.data.frame(sp_points@data)
 #species_count <- sp_point_info %>% group_by(Species) %>% tally()
#Tetra_count <- sp_point_info %>% group_by(basin_id)%>% summarise(sp_point_info$Species)


# #directly with species points > with cut data set (1980-2018 & less variables)
# sp_points <- shapefile(paste0(DIR, "/fish_sdm.shp"))
# head(sp_points)
# #colnames(sp_points@data)[1] <- "cat"
# sp_points$Species <- gsub(" ", ".", sp_points$Species)  # remove whitespace, biomod will remove it later...


# #directly with species points >all datapoints unchanged, just basin extracted
# sp_points3 <- shapefile(paste0(DIR, "/fish_all_2.shp"))
# head(sp_points3)
# #colnames(sp_points@data)[1] <- "cat"
# sp_points3$Species <- gsub(" ", ".", sp_points3$Species)  # remove whitespace, biomod will remove it later...

### Prepare species data for pseudo-absences
my_species <- c("Acestrocephalus.sardina",
                "Corydoras.septentrionalis",
                "Metynnis.argenteus", 
                "Serrasalmus.rhombeus",
                "Tetragonopterus.argenteus") 

#my_species <- c(spec_name$species)
### Run thorugh all species in a loop 
# sp_dt_for_biomod <- na.omit(sp_dt[,"basin_id"]) # for storing species data in the loop

sp_dt <- as.data.table(basin_id, na.rm=T)
sp_dt <- sp_dt[!duplicated(sp_dt),]

for (i in 1:length(my_species)) {
  cat("Running species", i, "\n")
  tmp <- as.data.table(extract(basin_id, subset(sp_points, species==my_species[i]) , df=T))
  tmp <- tmp[!is.na(tmp$basin_id),]
  tmp <- tmp[!duplicated(tmp$basin_id),] # remove duplicates
  tmp$spec <- 1 # write species occurrence
  names(tmp) <- c("ID", "basin_id", my_species[i] )
  sp_dt <- merge(sp_dt, tmp[,c("basin_id",  my_species[i]), with=F] , by="basin_id", all.x=T)
  rm(tmp)
}


# test <- sp_points[sp_points$Species=="Abramites.equesi"]
# test2 <-base::subset(sp_points, sp_points$Species== "Abramites.eques") 
# test3 <-base::subset(sp_points, sp_points$Species== "Tetragonopterus.argenteus")
#test4 <-base::subset(sp_points2, sp_points2$Species== "Abramites.equesi") 

### Count number of points per species and basin for modelling
data.frame(number_observations=colSums(sp_dt[,-1], na.rm=T))

###--------------------------------#
###--- Load environmental data ----
###--------------------------------#

##=====present data=====##
env_data <- fread(paste0(DIR, "/bio_lc_stream_present", ".csv"), h=T)
### Select the variables
var_list <- names(env_data[,-c("basin_id")])

### Attach species data --> this table now contains all data
all_data <- merge(env_data, sp_dt, by="basin_id", all.x=T) # add species data
head(all_data)
tail(all_data)
str(all_data)

### Remove duplicates (if any)
all_data <- all_data[!duplicated(all_data$basin_id),] # all removed earlier
summary(all_data)
str(all_data)
length(unique(all_data$basin_id))


###-------------------------------#
###--- Run an ensemble model -----
###-------------------------------#

### Create grid_id layer --> needed later for getting the raster predictions
domain_cells <- as.data.table.raster(basin_id)
domain_cells$seq_id <- seq.int(1:nrow(domain_cells))
setkey(domain_cells, seq_id) # set key



### Download Maxent from within R
download.file(url = "https://biodiversityinformatics.amnh.org/open_source/maxent/maxent.php?op=download", destfile = "maxent.zip", mode='wb')
### Unzip
MAXENT_DIR=paste0(DIR, "/MAXENT_DIR")
unzip("maxent.zip", exdir=MAXENT_DIR, overwrite = T)

### Maxent citation
# Steven J. Phillips, Miroslav Dud?k, Robert E. Schapire. [Internet] Maxent software for modeling species niches and distributions (Version 3.4.1). Available from url: http://biodiversityinformatics.amnh.org/open_source/maxent/. Accessed on 2020-5-13.

# wget -O maxent.zip "https://biodiversityinformatics.amnh.org/open_source/maxent/maxent.php?op=download"
# unzip -j maxent.zip





###--- Prepare the data for biomod2----
SPP_NA <- all_data[,..my_species] # get table of species data
EXP <- all_data[,..var_list]
COORD <- cbind(all_data[,"basin_id"], all_data[,"basin_id"]); names(COORD) <- c("x", "y")
# save(all_data, var_list, sp_dt, SPP_NA, EXP, COORD, file=paste0(DIR, "/precooked_species_variables_coordinates.RData"))

# SPP_NA <- all_data[,c("obs.spp.1.NA", "obs.spp.2.NA", "obs.spp.3.NA")]
# EXP <- all_data[,c("slope",  "forest", "urban", "bio1", "bio7", "log_flow")]

#For Future projections
#EXP_PAST <- SSP1_50_CCSM[,..var_list] # for now only dummy data




###--- Select if stream network distance calculations should be done ----
DO_DISTANCE=TRUE
# DO_DISTANCE=FALSE # set to false if it should not be done


### Setup parallel computation
cl <- parallel::makeCluster(N_CORES, outfile = "")
registerDoParallel(cl)

### test without loop: 
# i=2 

###--- Start loop - this runs through all species and saves all results for later lookup ----
# for(i in 1:length(SPP_NA)) {


###--- Run in parallel ----  
foreach(i=1:length(SPP_NA), 
        .errorhandling = c("pass"), 
        .verbose = T, 
        .packages = c("raster", "plyr", "dplyr", "data.table", "biomod2", "scales", "viridis")) %dopar% {
          # .errorhandling = c("stop", "remove", "pass")
          
          
          
          ### Prepare species name and folder
          MY_SPP <- names(SPP_NA[,i,with=F]) # get species name
          MY_SPP_DIR <- paste0(DIR, "/", MY_SPP) # create species folder path
          system(paste0("rm -rf ", MY_SPP_DIR))   # remove folder (previous model run)
          dir.create(MY_SPP_DIR)  # create species folder
          
          ### Copy maxent folder into each species' directory
          file.copy(MAXENT_DIR, MY_SPP_DIR, recursive=TRUE)
          MAXENT_DIR_TMP <- paste0(MY_SPP_DIR, "/MAXENT_DIR/maxent")
          
          cat(">>>>>>  Running models for", MY_SPP, "   <<<<<<", "\n")
          
          
          
          
          ##-----------------------------------------------------------#
          ####--- Calculate stream-network distance between points ----
          ##-----------------------------------------------------------#
          
          if(DO_DISTANCE==TRUE) { 
            
            ### Load stream network raster
            str_net <- raster(paste0(DIR, "/stream_network_110_study_area", ".tif")) 
            names(str_net) <- "basin_id"
            #x11(); plot(str_net)
            #values_stream <-na.omit(unique(values(str_net)))
            
            ### Add a "1" is species occurs in the basin that corresponds to the stream reach (no "snapping" needed)
            str_net_tmp <- as.data.table(str_net, na.rm=F) # get data.table
           
            
            # str_id <- na.omit(str_net_tmp)
            # str_id <- str_id[!duplicated(str_id$basin_id),] 
            # write.csv(str_id, "/mnt/brunner/brunner/stream_network/basin_ids_stream_network.csv", row.names = F)
            
            str_net_tmp$seq_id <- seq.int(1:nrow(str_net_tmp)) # add ID
            setkey(str_net_tmp, seq_id)
            
            ### Prepare species data
            tmp_spp <- sp_dt[,c("basin_id", MY_SPP), with=F] #?sp_dt only 1 or NA?#
            tmp_spp <- na.omit(tmp_spp)
            set_these_to_1 <- sort(unique(tmp_spp$basin_id)) # get those basins that have species observations
            
            ### Set the values in the stream network
            str_net_tmp$for_distance <- ifelse(is.na(str_net_tmp$basin_id), NA, 0) # set all stream cells to zero, land remains NA
            
            # str_distance <-str_net_tmp$for_distance
            # str_distance <-as.data.frame(str_distance)  
            # str_distance <- na.omit(str_distance)
            
            ### Get index of "set_these_to_1" basins, and repalce those with "1"
            my_index <- which(str_net_tmp$basin_id  %in% set_these_to_1) #alle rasterzellen in basins wo arten vorkommen
            str_net_tmp$for_distance <- replace(str_net_tmp$for_distance, my_index, 1) # write "1" in those/all pixels of basins that have species
            str_net_tmp <- as.data.table(arrange(str_net_tmp, seq_id))
            #distance_count <- str_net_tmp %>% group_by(for_distance)%>% tally()
             # subset(str_net_tmp, str_net_tmp$basin_id==8019) # check
            
            #  str_distance2 <-str_net_tmp$for_distance
            #  str_distance2 <-as.data.frame(str_distance2)  
            # str_distance2 <- na.omit(str_distance2)
            
            ### Insert the new values into the raster (0=no species, 1=species observed in that basin/stream)
            str_net_for_dist <- setValues(basin_id, str_net_tmp$for_distance)
            # writeRaster(str_net_for_dist, paste0(DIR, "/temp_distance_file.tif"), overwrite=T) # check 
            
            ### Calculate the within-stream distance between points
            str_net_for_dist <- gridDistance(str_net_for_dist, 1, omit=NA)
            writeRaster(str_net_for_dist, paste0(MY_SPP_DIR, "/", MY_SPP, "_str_net_distance.tif"), options= c("COMPRESS=DEFLATE", "ZLEVEL=9"), overwrite=T)
            
            
            ### Assign the values to the entire basin raster file
            dist_per_bas <- as.data.table(as.data.frame(zonal(str_net_for_dist, basin_id, fun='mean', digits=2, na.rm=T)))
            names(dist_per_bas) <- c("basin_id", "distance")
            
            ### Rescale
            dist_per_bas$distance_scaled <-  rescale(dist_per_bas$distance, to = c(1, 0)) 
            dist_per_bas$distance_scaled <- ifelse(is.na(dist_per_bas$distance_scaled), -1, dist_per_bas$distance_scaled) # -1 if non-connected 
            
            ### Get as raster
            basin_table <- as.data.table(basin_id, na.rm=F) # get basins in a table
            basin_table$seq_id <- seq.int(1:nrow(basin_table))
            setkey(basin_table, seq_id)
            basin_for_dist <- merge(basin_table, dist_per_bas, by="basin_id", all.x=T)
            basin_for_dist <- as.data.table(arrange(basin_for_dist, seq_id))
            basin_for_dist_r <- setValues(basin_id, basin_for_dist$distance_scaled)
            writeRaster(basin_for_dist_r,  paste0(MY_SPP_DIR, "/", MY_SPP, "_basin_distance.tif"), options= c("COMPRESS=DEFLATE", "ZLEVEL=9"), overwrite=T)
            save(basin_for_dist, file=paste0(MY_SPP_DIR, "/", MY_SPP, "_basin_distance_weight_table.RData")) # distance=raw, disctance_scaled= -1...0-1
            
            # ### Attach as a variable for the model --> don't do this other than for testing purpose only... 
            ### all_data <- merge(all_data, dist_per_bas, by="basin_id")
            ### EXP <- cbind(EXP, all_data[,"distance_scaled"])
          }
          
          
          
          
          ###-------------------#
          ###--- Run models ----
          ###-------------------#
          
          ### Initiate
          myBiomodData <- BIOMOD_FormatingData(resp.var = SPP_NA[,i,with=F] ,
                                               expl.var = EXP,
                                               resp.xy =  COORD, # NULL
                                               resp.name = MY_SPP,
                                               PA.nb.rep = 1,
                                               PA.nb.absences = 10000, 
                                               PA.strategy = 'random')
          
          
          ### Options definition
          # myBiomodOption <- BIOMOD_ModelingOptions() # default options
          
          myBiomodOption <- BIOMOD_ModelingOptions(
            RF = list( do.classif = TRUE,
                       ntree = 1000, # default 500, models more stable with 1000 trees?
                       mtry = 'default',
                       nodesize = 5,
                       maxnodes = NULL),
            
            # NOTE: Maxent needs x11 connextion even without GUI...
            MAXENT.Phillips = list(
              path_to_maxent.jar = MAXENT_DIR_TMP, # MAXENT_DIR,
              # path_to_maxent.jar = "C:\temp\spdep_SDMs\MAXENT_DIR", # in windows
              memory_allocated = 4096, # give more RAM
              background_data_dir = 'default',
              maximumbackground = 'default',
              maximumiterations = 200,
              visible = FALSE,
              linear = TRUE,
              quadratic = TRUE,
              product = TRUE,
              threshold = FALSE, # edited
              hinge = FALSE, # edited
              lq2lqptthreshold = 80,
              l2lqthreshold = 10,
              hingethreshold = 15,
              beta_threshold = -1,
              beta_categorical = -1,
              beta_lqp = -1,
              beta_hinge = -1,
              betamultiplier = 1,
              defaultprevalence = 0.5)
          )
          
          
          ### Modelling
          myBiomodModelOut <- BIOMOD_Modeling(
            myBiomodData,
            #models = c('SRE','CTA','RF','MARS','FDA','MAXENT.Phillips','GLM','GAM','GBM','ANN'),
            #models = c('SRE','CTA','RF','FDA','MAXENT.Phillips','GLM','GAM','GBM','ANN'),# Error for MARS: <simpleError in validObject(.Object): invalid class “MARS_biomod2_model” object: FALSE>
            #models = c('SRE','CTA','RF', 'FDA','MAXENT.Phillips', 'GLM'), # fast ones
            # models = c('RF', 'MAXENT.Phillips'), #  Machine learning subset
            # models = c('RF','CTA', 'FDA', 'MAXENT.Phillips'), #  Machine learning and classification subset
            # models = c('RF', 'MAXENT.Phillips'),
            # models = c('MAXENT.Phillips'),
            #models = c('CTA','RF','FDA','MAXENT.Phillips','GLM','GBM'),
            #models = c('CTA','RF','MAXENT.Phillips','GLM','GBM'),
            models = c('CTA','MAXENT.Phillips','GLM','GBM'),
            models.options = myBiomodOption,
            NbRunEval=2, # repetitions: use 2 for testing, 10 for the final run
            DataSplit=70,
            Yweights=NULL,
            VarImport=3,
            models.eval.meth = c('TSS','ROC'),
            SaveObj = TRUE,
            rescal.all.models = TRUE)
          
          gc()
          
          ### Ensemble model
          myBiomodEM <- BIOMOD_EnsembleModeling(
            modeling.output = myBiomodModelOut,
            chosen.models = 'all',  # only final model?
            em.by = 'all',
            eval.metric = c('TSS'),
            eval.metric.quality.threshold = c(0.4), # e.g. 0.4
            prob.mean = F,
            prob.cv = F,
            prob.ci = F,
            prob.ci.alpha = 0.05,
            prob.median = F,
            committee.averaging = F,
            prob.mean.weight = T,
            prob.mean.weight.decay = "proportional" ) # 1.6
          
          ###---- Project on present-day data
          myBiomodProj <- BIOMOD_Projection(
            modeling.output = myBiomodModelOut,
            new.env = EXP,
            proj.name = 'present',
            selected.models = 'all',
            binary.meth= 'TSS', # NULL
            compress = 'xz',
            clamping.mask = T)
          
          ### Do ensemble-models projections on present variables
          # ?BIOMOD_EnsembleForecasting
          myBiomodEF <- BIOMOD_EnsembleForecasting(
            projection.output = myBiomodProj,
            EM.output = myBiomodEM,
            binary.meth = 'TSS', # NULL
            total.consensus = TRUE)
          
          gc() 
          
          
          ### Save the ensemble evaluation metrics separately:
          (x <- get_evaluations(myBiomodEM))
          write.csv(x, paste0(MY_SPP_DIR, "/", MY_SPP, "_ensemble_evaluation.csv")); rm(x)
          
          
          ### Load predictions - probability
          load(paste0(MY_SPP_DIR, "/proj_present/proj_present_", MY_SPP, "_ensemble.RData"))
          tmp_proj_present <- as.data.table(as.data.frame(ef.out)); rm(ef.out)
          names(tmp_proj_present) <- "probs"
          tmp_proj_present <- tmp_proj_present/1000
          tmp_proj_present <- cbind(all_data[,c("basin_id")], tmp_proj_present)
          
          
          ### Load predictions - binary TSS
          load(paste0(MY_SPP_DIR, "/proj_present/proj_present_", MY_SPP, "_ensemble_TSSbin.RData"))
          ef.out_bin <- get(paste0("proj_present_",MY_SPP, "_ensemble_TSSbin"))
          tmp_proj_present_bin <- as.data.table(as.data.frame(ef.out_bin)); rm(ef.out_bin)
          names(tmp_proj_present_bin) <- "binary"
          tmp_proj_present_bin <- tmp_proj_present_bin
          tmp_proj_present_bin <- cbind(all_data[,c("basin_id")], tmp_proj_present_bin)
          
          
          
          
          ### Save all predictions as a table - useful for e.g. Marxan later on..
          all_predictions <- merge(tmp_proj_present, tmp_proj_present_bin, by="basin_id")
          all_predictions$semi_binary <- ifelse(all_predictions$binary==1, all_predictions$probs, 0)
          save(all_predictions, file=paste0(MY_SPP_DIR, "/", MY_SPP, "_all_predictions_basin_table.RData"))
          ## if distance_weights needed: load the separate RData file and attach to this table in another loop...
          
          
          ### Prepare maps
          tmp_proj_present_domain <- merge(domain_cells, all_predictions, by="basin_id", all.x=T)
          ### Sort data to match the spatial configuration
          tmp_proj_present_domain <- as.data.table(arrange(tmp_proj_present_domain, seq_id))
          ### Insert the values into the raster template
          pred_ensemble_r <- setValues(basin_id, tmp_proj_present_domain$probs)
          names(pred_ensemble_r) <- "pred_ensemble"
          
          pred_ensemble_bin_r <- setValues(basin_id, tmp_proj_present_domain$binary)
          names(pred_ensemble_bin_r) <- "pred_ensemble_binary"
          
          pred_ensemble_semi_bin_r <- setValues(basin_id, tmp_proj_present_domain$semi_binary)
          names(pred_ensemble_semi_bin_r) <- "pred_ensemble_semi_binary"
          
          ### Write raster to disk
          writeRaster(pred_ensemble_r, paste0(MY_SPP_DIR, "/", MY_SPP, "_ensemble_prob.tif"), options= c("COMPRESS=DEFLATE", "ZLEVEL=9") , overwrite=T)
          writeRaster(pred_ensemble_bin_r, paste0(MY_SPP_DIR, "/", MY_SPP, "_ensemble_bin.tif"), options= c("COMPRESS=DEFLATE", "ZLEVEL=9") , overwrite=T)
          writeRaster(pred_ensemble_semi_bin_r, paste0(MY_SPP_DIR, "/", MY_SPP, "_ensemble_semi_bin.tif"), options= c("COMPRESS=DEFLATE", "ZLEVEL=9") , overwrite=T)
          
          ### Get species points
          # mypoints <- as.data.table(as.data.frame(subset(sp_points, species==MY_SPP)))
          
          
          ### Save maps as pdf
          pdf(paste0(MY_SPP_DIR, "/", MY_SPP, "_ensemble_prob.pdf"), width=20, height=20)
          plot(pred_ensemble_r, col=viridis(15), zlim=c(0,1))
          points(subset(sp_points, species==MY_SPP), pch=16, col="red")
          dev.off()
          
          pdf(paste0(MY_SPP_DIR, "/", MY_SPP, "_ensemble_bin.pdf"), width=20, height=20)
          plot(pred_ensemble_bin_r, col=viridis(15), zlim=c(0,1))
          points(subset(sp_points, species==MY_SPP), pch=16, col="red")
          dev.off()
          
          pdf(paste0(MY_SPP_DIR, "/", MY_SPP, "_ensemble_semi_bin.pdf"), width=20, height=20)
          plot(pred_ensemble_semi_bin_r, col=viridis(15), zlim=c(0,1))
          points(subset(sp_points, species==MY_SPP), pch=16, col="red")
          dev.off()
          
          
          
          
          
          ###--- Multiply the predictions with the distance map (clipping)----
          if(DO_DISTANCE==TRUE) {
            # pred_ensemble_r_dist <- pred_ensemble_r * basin_for_dist_r
            # pred_ensemble_r_dist[pred_ensemble_r_dist<=0] <- 0
            # ### Save map as pdf
            # pdf(paste0(MY_SPP_DIR, "/", MY_SPP, "_ensemble_distance.pdf"), width=20, height=20)
            # plot(pred_ensemble_r_dist, col=viridis(15), zlim=c(0,1))
            # points(subset(sp_points, species==MY_SPP), pch=16, col="red")
            # dev.off()
            
            ### Export ratser
            # writeRaster(pred_ensemble_r_dist, paste0(MY_SPP_DIR, "/", MY_SPP, "_ensemble_prob_distance.tif"), overwrite=T)
            # rm(pred_ensemble_r_dist) # clean up
            
            ### Save the distance-weighted predictions
            all_predictions <- merge(all_predictions, dist_per_bas, by="basin_id") # attach basin distances
            ### Attach the different distance-corrected 
            all_predictions$tmp <- all_predictions$distance_scaled * all_predictions$probs # multuiply probabilities with distance
            all_predictions$probs_dist <- ifelse(all_predictions$tmp <=0, 0 , all_predictions$tmp) # omit if negative (non-connected basin)
            all_predictions$binary_dist <- ifelse(all_predictions$tmp <=0, 0 , all_predictions$binary) 
            all_predictions$semi_binary_dist <- ifelse(all_predictions$tmp <=0, 0 , all_predictions$semi_binary) 
            all_predictions <- subset(all_predictions, select=-c(tmp))
            save(all_predictions, file=paste0(MY_SPP_DIR, "/", MY_SPP, "_all_predictions_basin_table.RData"))
            
            ### Prepare maps
            tmp_proj_present_domain <- merge(domain_cells, all_predictions, by="basin_id", all.x=T)
            ### Sort data to match the spatial configuration
            tmp_proj_present_domain <- as.data.table(arrange(tmp_proj_present_domain, seq_id))
            ### Insert the values into the raster template
            pred_ensemble_dist_r <- setValues(basin_id, tmp_proj_present_domain$probs_dist)
            names(pred_ensemble_dist_r) <- "pred_ensemble_dist"
            
            pred_ensemble_bin_dist_r <- setValues(basin_id, tmp_proj_present_domain$binary_dist)
            names(pred_ensemble_bin_dist_r) <- "pred_ensemble_binary_dist"
            
            pred_ensemble_semi_bin_dist_r <- setValues(basin_id, tmp_proj_present_domain$semi_binary_dist)
            names(pred_ensemble_semi_bin_dist_r) <- "pred_ensemble_semi_binary_dist"
            
            ### Write raster to disk
            writeRaster(pred_ensemble_dist_r, paste0(MY_SPP_DIR, "/", MY_SPP, "_ensemble_prob_dist.tif"), options= c("COMPRESS=DEFLATE", "ZLEVEL=9") , overwrite=T)
            writeRaster(pred_ensemble_bin_dist_r, paste0(MY_SPP_DIR, "/", MY_SPP, "_ensemble_bin_dist.tif"), options= c("COMPRESS=DEFLATE", "ZLEVEL=9") , overwrite=T)
            writeRaster(pred_ensemble_semi_bin_dist_r, paste0(MY_SPP_DIR, "/", MY_SPP, "_ensemble_semi_bin_dist.tif"), options= c("COMPRESS=DEFLATE", "ZLEVEL=9") , overwrite=T)
            
            ### Get species points
            # mypoints <- as.data.table(as.data.frame(subset(sp_points, species==MY_SPP)))
            
            
            ### Save maps as pdf
            pdf(paste0(MY_SPP_DIR, "/", MY_SPP, "_ensemble_prob_dist.pdf"), width=20, height=20)
            plot(pred_ensemble_dist_r, col=viridis(15), zlim=c(0,1))
            points(subset(sp_points, species==MY_SPP), pch=16, col="red")
            dev.off()
            
            pdf(paste0(MY_SPP_DIR, "/", MY_SPP, "_ensemble_bin_dist.pdf"), width=20, height=20)
            plot(pred_ensemble_bin_dist_r, col=viridis(15), zlim=c(0,1))
            points(subset(sp_points, species==MY_SPP), pch=16, col="red")
            dev.off()
            
            pdf(paste0(MY_SPP_DIR, "/", MY_SPP, "_ensemble_semi_bin_dist.pdf"), width=20, height=20)
            plot(pred_ensemble_semi_bin_dist_r, col=viridis(15), zlim=c(0,1))
            points(subset(sp_points, species==MY_SPP), pch=16, col="red")
            dev.off()
            
          }
          
          
          
          ### Get variable importance
          var_imp <- as.data.frame(variables_importance(myBiomodEM@em.models[[1]], data=EXP , method="full_rand" , nb_rand=3)$mat)
          var_imp$mean <- rowMeans(var_imp)
          tmp <- subset(var_imp, select=c(mean))
          var_imp$percent_contribution <- round(apply(tmp, 1, function(x) x/colSums(tmp) ),2)  # scale from 0-1
          var_imp$variable <- row.names(var_imp)
          write.csv(var_imp, paste0(MY_SPP_DIR, "/", MY_SPP, "_variable_importance.csv"), row.names = F)
          
          ### Write raster to disk
          writeRaster(pred_ensemble_r, paste0(MY_SPP_DIR, "/", MY_SPP, "_ensemble.tif"), options= c("COMPRESS=DEFLATE", "ZLEVEL=9"), overwrite=T)
          writeRaster(pred_ensemble_bin_r, paste0(MY_SPP_DIR, "/", MY_SPP, "_ensemble_bin.tif"), options= c("COMPRESS=DEFLATE", "ZLEVEL=9"), overwrite=T)
          
          gc()
          
          
          
          
          ###--- Check single algorithm results ----
          myCurrentProj <- get_predictions(myBiomodProj)
          ### Get single "full" model maps
          single_full_pred <- as.data.table(as.data.frame(myCurrentProj[,,"Full","PA1"]))
          single_full_pred <- single_full_pred/1000 # rescale to 0-1
          single_full_pred <- cbind(all_data[,c("basin_id")], single_full_pred)
          
          ### Prepare map
          single_full_pred_domain <- merge(domain_cells, single_full_pred, by="basin_id", all.x=T)
          ### Sort data to match the spatial configuration
          single_full_pred_domain <- as.data.table(arrange(single_full_pred_domain, seq_id))
          
          
          ### Create character vector of all possible algorithms  
          check_these <- c("SRE", "CTA","RF","MARS","GLM","GAM","GBM","ANN","MAXENT.Phillips","MAXENT.Phillips.2")
          mystack <- stack() # empty stack for storing layers
          
          ### Run through all columns, check if model exists and write into raster stack
          for (k in check_these) {
            if(k %in% colnames(single_full_pred_domain)) { 
              tmp <- setValues(basin_id, single_full_pred_domain[,get(k)])
              names(tmp) <- k 
              mystack <- addLayer(mystack, tmp) 
            }
          }
          
          ### Save single algorighm map as pdf
          pdf(paste0(MY_SPP_DIR, "/", MY_SPP, "_single_algorithms.pdf"), width=20, height=20)
          plot(mystack, col=viridis(15), zlim=c(0,1))
          dev.off()
          
          
          ### Export the raw, non-scaled and non-normalized variable importances of each algorthm
          ### careful interpretation is "the higher the value, the more influence it has"
          var_imp_single <- as.data.frame(myBiomodModelOut@variables.importances@val[,,"Full", "PA1"])
          var_imp_single$variable <- row.names(var_imp_single)
          write.csv(var_imp_single, paste0(MY_SPP_DIR, "/", MY_SPP, "_single_algorithm_raw_variable_importance.csv"), row.names = F )
          
          ### Clean up
          rm(mystack, single_full_pred, single_full_pred_domain, var_imp_single)
          gc()
          
          
          
          
          # ###----------------------------------------------------------------------------#
          # ###---- Project on past/future data - if applicable --> else comment this block
          # ###----------------------------------------------------------------------------#
          # 
          # 
          
          SSPs <- c("SSP1", "SSP2", "SSP3", "SSP4", "SSP5") #, "SSP2", "SSP3", "SSP4", "SSP5"
          GCMs <- c("CCSM", "IPSL", "MIROC") #, "IPSL", "MIROC"
          YEAR <- c("2050", "2070") #, "2070"
        
          
          for (z in YEAR) {
            
            cat(">>>> Running year", z, "\n")


            
                for (w in SSPs) {

                  cat(">>>> Running SSP", w, "\n")
        
            
            
                      # inner loop
                      for (l in GCMs) {
                        
                        cat(">>>>Running GCM", l, "\n")
                        my_scenario <- paste0(w, "_", l, "_", z)
                            dat <- fread(paste0(DIR, "/SRES_data/", my_scenario, ".csv"))
                            dat_list <- names(dat[,-c("basin_id")])
                            fut_data <- merge(dat, sp_dt, by="basin_id", all.x=T)
                            fut_data <- fut_data[,..dat_list]
  
                            myBiomodProj_future<- BIOMOD_Projection(
                              modeling.output = myBiomodModelOut,
                              new.env = fut_data,
                              proj.name = my_scenario,
                              selected.models = 'all',
                              binary.meth= 'TSS', # NULL
                              compress = 'xz',
                              clamping.mask = F)
                            
                            myBiomodEF_future <- BIOMOD_EnsembleForecasting(
                              projection.output = myBiomodProj_future,
                              EM.output = myBiomodEM,
                              binary.meth = 'TSS', # NULL
                              total.consensus = TRUE)
                            
                            ### Load predictions - probability
                            load(paste0(MY_SPP_DIR, "/proj_", my_scenario, "/proj_", my_scenario, "_", MY_SPP, "_ensemble.RData"))
                            tmp_proj_future <- as.data.table(as.data.frame(ef.out)); rm(ef.out)
                            names(tmp_proj_future) <- "probs"
                            tmp_proj_future <- tmp_proj_future/1000
                            tmp_proj_future <- cbind(all_data[,c("basin_id")], tmp_proj_future)
                            
                            ### Prepare map
                            tmp_proj_future_domain <- merge(domain_cells, tmp_proj_future, by="basin_id", all.x=T)
                            ### Sort data to match the spatial configuration
                            tmp_proj_future_domain <- as.data.table(arrange(tmp_proj_future_domain, seq_id))
                            ### Insert the values into the raster template
                            pred_ensemble_future_r <- setValues(basin_id, tmp_proj_future_domain$probs)
                            names(pred_ensemble_future_r) <- "pred_ensemble"
                            
                            ### Save map as pdf
                            pdf(paste0(MY_SPP_DIR, "/", MY_SPP, "_ensemble_prob_future","_", my_scenario, ".pdf"), width=20, height=20)
                            plot(pred_ensemble_future_r, col=viridis(15), zlim=c(0,1))
                            dev.off()
                            
                            
                            ### Load predictions - binary TSS
                            load(paste0(MY_SPP_DIR, "/proj_", my_scenario, "/proj_", my_scenario, "_", MY_SPP, "_ensemble_TSSbin.RData")) 
                            ef.out_bin <- get(paste0("proj_", my_scenario, "_", MY_SPP, "_ensemble_TSSbin"))
                            
                            tmp_proj_future_bin <- as.data.table(as.data.frame(ef.out_bin)); rm(ef.out_bin)
                            names(tmp_proj_future_bin) <- "binary"
                            tmp_proj_future_bin <- tmp_proj_future_bin
                            tmp_proj_future_bin <- cbind(all_data[,c("basin_id")], tmp_proj_future_bin)
                            
                            ### Prepare map
                            tmp_proj_future_bin_domain <- merge(domain_cells, tmp_proj_future_bin, by="basin_id", all.x=T)
                            ### Sort data to match the spatial configuration
                            tmp_proj_future_bin_domain <- as.data.table(arrange(tmp_proj_future_bin_domain, seq_id))
                            ### Insert the values into the raster template
                            pred_ensemble_future_bin_r <- setValues(basin_id, tmp_proj_future_bin_domain$binary)
                            names(pred_ensemble_future_bin_r) <- "pred_ensemble_binary"
                            
                            
                            ### Save map as pdf
                            pdf(paste0(DIR, "/",MY_SPP, "/", MY_SPP, "_ensemble_binary_future.pdf"), width=20, height=20) 
                            plot(pred_ensemble_future_bin_r, col=viridis(15), zlim=c(0,1))
                            dev.off()
                            
                            ### Write raster to disk
                            writeRaster(pred_ensemble_future_r, paste0(DIR, "/",MY_SPP, "/", MY_SPP, "_ensemble_prob_", my_scenario, ".tif"), options= c("COMPRESS=DEFLATE", "ZLEVEL=9"), overwrite=T) # uncomment for full run
                            writeRaster(pred_ensemble_future_bin_r, paste0(DIR, "/",MY_SPP, "/", MY_SPP, "_ensemble_bin_future" ,"_", my_scenario, ".tif"), options= c("COMPRESS=DEFLATE", "ZLEVEL=9"), overwrite=T) # uncommentfor full run
                            
                            ### Save all predictions as a table - useful for e.g. Marxan later on..
                            all_future_predictions <- merge(tmp_proj_future, tmp_proj_future_bin, by="basin_id")
                            all_future_predictions$semi_binary <- ifelse(all_future_predictions$binary==1, all_future_predictions$probs, 0)
                            
                            
                            ### Save the distance-weighted predictions
                            all_future_predictions <- merge(all_future_predictions, dist_per_bas, by="basin_id") # attach basin distances
                            ### Attach the different distance-corrected
                            all_future_predictions$tmp <- all_future_predictions$distance_scaled * all_future_predictions$probs # multuiply probabilities with distance
                            all_future_predictions$probs_dist <- ifelse(all_future_predictions$tmp <=0, 0 , all_future_predictions$tmp) # omit if negative (non-connected basin)
                            all_future_predictions$binary_dist <- ifelse(all_future_predictions$tmp <=0, 0 , all_future_predictions$binary)
                            all_future_predictions$semi_binary_dist <- ifelse(all_future_predictions$tmp <=0, 0 , all_future_predictions$semi_binary)
                            all_future_predictions <- subset(all_future_predictions, select=-c(tmp))
                            save(all_future_predictions, file=paste0(MY_SPP_DIR, "/", MY_SPP, "_all_future_predictions_basin_table" ,"_", my_scenario, ".RData")) 
                            ### Clean up
                            rm(dat,
                               fut_data,
                               myBiomodProj_future,
                               myBiomodEF_future,
                               tmp_proj_future,
                               pred_ensemble_future_r,
                               tmp_proj_future_bin)
                            
                             } # exit inner loop


                  } # exit outer loop

          } # close YEAR

          
          
          #dat <- fread(paste0(DIR, "/SRES_data/", "SSP1", "_", "CCSM", "_", "2050", ".csv"))
          
                    # 
          #
          # ###---- END past/future projections
          
          
          
          
          
          
          ### Save all files: --> not needed, all files wre written already to disk
          # save.image(paste0(MY_SPP_DIR, "/", MY_SPP, "_all_results.Rdata"))
          # save(myBiomodEF, 
          #      myBiomomodProj, 
          #      myBiomodEM, 
          #      myBiomodModelOut, 
          #      myBiomodData, 
          #      MY_SPP, 
          #     file=paste0(MY_SPP_DIR, "/", MY_SPP, "_all_results.Rdata"))
          
          ### Remove files for next species run
          rm(pred_ensemble_r, 
             pred_ensemble_bin_r, 
             tmp_proj_present, 
             tmp_proj_present_bin, 
             pred_ensemble_dist_r,
             pred_ensemble_bin_dist_r,
             pred_ensemble_semi_bin_dist_r,
             myBiomodEF, 
             myBiomodProj, 
             myBiomodEM, 
             myBiomodModelOut, 
             myBiomodData, 
             MY_SPP,
             all_predictions)
          
          MAXENT_DIR_TMP <- paste0(MY_SPP_DIR, "/MAXENT_DIR")
          system(paste0("rm -rf ", MAXENT_DIR_TMP))   # remove maxent folder
            
          gc() 
          
          #keep the files in a vector before deleting everything
          # keep= c("SPP_NA", "EXP", "DIR", "MAXENT_DIR")
          # remove = ls()
          # 
          # #Delete everything but the object in 'keep'.
          # rm(list=(remove[is.na(match(remove, keep))]))  
          
        } # end loop

  stopCluster(cl) # close the parrallel backend





### Useful functions to explore:
# niche overlap among species tables, https://cran.r-project.org/web/packages/fuzzySim/fuzzySim.pdf for modOverlap() function, see also here https://modtools.wordpress.com/2015/10/30/modoverlap/
#
# https://rdrr.io/cran/biomod2/man/BIOMOD_RangeSize.html for range change estimates











###---- Plot "nice" maps with species occurrences -----

### Species name
MY_SPP="Tetragonopterus.argenteus"

### Create grid_id layer --> needed later for getting the raster predictions
basin_id <- raster(paste0(DIR, "/basins_study_area", ".tif")); names(basin_id) <- "basin_id"
domain_cells <- as.data.table.raster(basin_id)
domain_cells$seq_id <- seq.int(1:nrow(domain_cells))
setkey(domain_cells, seq_id) # set key


### Subset species points
sp_points <- shapefile(paste0(DIR, "/fish_sdm.shp"))
mypoints <- as.data.table(as.data.frame(subset(sp_points, Species==MY_SPP)))


### Load predictions - probability
load(paste0(DIR, "/", MY_SPP, "/proj_present/proj_present_", MY_SPP, "_ensemble.RData"))
tmp_proj_present <- as.data.table(as.data.frame(ef.out)); rm(ef.out)
names(tmp_proj_present) <- "probs"
tmp_proj_present <- tmp_proj_present/1000
tmp_proj_present <- cbind(all_data[,c("basin_id")], tmp_proj_present)

### Prepare map
tmp_proj_present_domain <- merge(domain_cells, tmp_proj_present, by="basin_id", all.x=T)
### Sort data to match the spatial configuration
tmp_proj_present_domain <- as.data.table(arrange(tmp_proj_present_domain, seq_id))
### Insert the values into the raster template
r <- setValues(basin_id, tmp_proj_present_domain$probs)
names(r) <- "pred_ensemble"

### Shortcut if raster was saved earlier
# r <- raster(paste0(DIR, "/",MY_SPP, "/", MY_SPP, "_ensemble.tif"))


mymap <- as.data.table(r, xy=T)
res=ncell(mymap)


### Plot the probability map
p <- ggplot() + 
  geom_tile(data=mymap, aes(x=x, y=y, fill=pred_ensemble)) + 
  scale_fill_viridis(na.value="transparent", limits=c(0,1)) +
  coord_equal() +
  theme_map() +
  theme(legend.key.width=unit(2, "cm")) +
  theme(panel.grid.major = element_blank(), # remove grey plot backgound
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.title.x = element_text(face="bold", size=1),
        axis.text.x  = element_text(size=10)) +
  theme(axis.title.y = element_text(face="bold", size=1),
        axis.text.y  = element_text(size=10)) +
  theme(axis.line = element_line(colour="black", size=1))+
  theme(axis.ticks.length=unit(.15, "cm"))+
  xlab("Longitude") + ylab("Latitude") +
  theme(legend.position="bottom")

### Add species points
# p <- p + geom_point(data=mypoints, aes(x=coords.x1, y=coords.x2), size=2, col="black",  fill=NA) # points
p <- p + geom_point(data=mypoints, aes(x=coords.x1, y=coords.x2), shape = 21, colour = "black", fill = "NA", size = 3, stroke = 1.5) # hollow circles


### Plot as pdf
pdf(paste0(DIR, "/", MY_SPP, "/pred_ensemble.pdf"))
plot(p)
dev.off()

### Plot as svg --> edit in Inkscape or similar
svg(paste0(DIR, "/", MY_SPP, "/pred_ensemble.svg"))
plot(p)
dev.off()







###--- Plot the binary map ----

### Load predictions - binary TSS
load(paste0(DIR, "/", MY_SPP, "/proj_present/proj_present_", MY_SPP, "_ensemble_TSSbin.RData"))
ef.out_bin <- get(paste0("proj_present_",MY_SPP, "_ensemble_TSSbin"))
tmp_proj_present_bin <- as.data.table(as.data.frame(ef.out_bin)); rm(ef.out_bin)
names(tmp_proj_present_bin) <- "binary"
tmp_proj_present_bin <- tmp_proj_present_bin
tmp_proj_present_bin <- cbind(all_data[,c("basin_id")], tmp_proj_present_bin)

### Prepare map
tmp_proj_present_domain <- merge(domain_cells, tmp_proj_present_bin, by="basin_id", all.x=T)
### Sort data to match the spatial configuration
tmp_proj_present_domain <- as.data.table(arrange(tmp_proj_present_domain, seq_id))
### Insert the values into the raster template
r <- setValues(basin_id, tmp_proj_present_domain$binary)
names(r) <- "pred_ensemble_bin"

### Shortcut
# r <- tmp_proj_present_bin
# r <- raster(paste0(DIR, "/",MY_SPP, "/", MY_SPP, "_ensemble_binary.tif"))

res=ncell(r)
mymap <- as.data.table(r, xy=T)


### Plot the binary map
p <- ggplot() + 
  geom_tile(data=mymap, aes(x=x, y=y, fill=pred_ensemble_bin)) + 
  scale_fill_viridis(na.value="transparent", limits=c(0,1)) +
  coord_equal() +
  theme_map() +
  theme(legend.key.width=unit(2, "cm")) +
  theme(panel.grid.major = element_blank(), # remove grey plot backgound
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.title.x = element_text(face="bold", size=1),
        axis.text.x  = element_text(size=10)) +
  theme(axis.title.y = element_text(face="bold", size=1),
        axis.text.y  = element_text(size=10)) +
  theme(axis.line = element_line(colour="black", size=1))+
  theme(axis.ticks.length=unit(.15, "cm"))+
  xlab("Longitude") + ylab("Latitude") +
  theme(legend.position="bottom")

### Add species points
# p <- p + geom_point(data=mypoints, aes(x=coords.x1, y=coords.x2), size=2, col="black",  fill=NA) # points
p <- p + geom_point(data=mypoints, aes(x=coords.x1, y=coords.x2), shape = 21, colour = "black", fill = "NA", size = 3, stroke = 1.5) # hollow circles



### Plot as pdf
pdf(paste0(DIR, "/", MY_SPP, "/pred_ensemble_binary.pdf"))
plot(p)
dev.off()

### Plot as svg --> edit in Inkscape or similar
svg(paste0(DIR, "/", MY_SPP, "/pred_ensemble_binary.svg"))
plot(p)
dev.off()






###--- Get mean / median elevation of species' predicted suitable habitats----

# - merge the table: "basin_id | binary_prediction" to the "all_data" table, subset the "prediction==1" and get mean(all_data$dem) 

###--- Create a bi-plot of e.g. flow and temperature, or other variables













### TEST

# ### Get variable importance manually - final output not yet checked!
# ### --> allows to check the single algorithms performance and predictor contributions
# 
# # range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))} # define functions
#
# load(paste0(DIR, "/", MY_SPP, "/all_results.Rdata")) # get model outputs
# var_imp <- myBiomodModelOut@variables.importances@val[,,"Full", "PA1"]
# 
# ### Change this depending on your selected algorithms:
#   ### "RF", "MAXENT.Phillips"
# mod_eval <-data.frame(algorithm=rep(c("RF","MAXENT.Phillips"), each=2, times=1),
#                       method=rep(c("TSS", "ROC"), times=2),
#                       rbind(myBiomodModelOut@ models.evaluation@val[,, "RF", "Full", "PA1"],
#                             myBiomodModelOut@ models.evaluation@val[,, "MAXENT.Phillips", "Full", "PA1"]) )
# 
# # mod_eval <-data.frame(algorithm=rep(c("RF","CTA"), each=2, times=1),
# #                       method=rep(c("TSS", "ROC"), times=2),
# #                       rbind(myBiomodModelOut@ models.evaluation@val[,, "RF", "Full", "PA1"],
# #                             myBiomodModelOut@ models.evaluation@val[,, "CTA", "Full", "PA1"]) )
# 
# 
# ### "RF", "GBM","MAXENT.Phillips"
# # mod_eval <-data.frame(algorithm=rep(c("RF", "GBM","MAXENT.Phillips"), each=2, times=1),
# #                       method=rep(c("TSS", "ROC"), times=3),
# #                       rbind(myBiomodModelOut@ models.evaluation@val[,, "RF", "Full", "PA1"],
# #                             myBiomodModelOut@ models.evaluation@val[,, "GBM", "Full", "PA1"],
# #                             myBiomodModelOut@ models.evaluation@val[,, "MAXENT.Phillips", "Full", "PA1"]) )
# 
# 
# ### Get the TSS values and assign as weights
# mod_eval_TSS <- subset(mod_eval, method=="TSS")
# mod_eval_TSS <- subset(mod_eval_TSS, select=c(algorithm, Testing.data))
# 
# var_imp_ensemble <- data.frame(var_imp=apply(var_imp, 1, function(x) {weighted.mean(x, mod_eval_TSS$Testing.data) })  )
# var_imp_ensemble$percent_contribution <- round(apply(var_imp_ensemble, 1, function(x) x/colSums(var_imp_ensemble) ),2)  # scale from 0-1
# var_imp_ensemble$variable <- row.names(var_imp_ensemble)
# write.csv(var_imp_ensemble, paste0(DIR, "/", MY_SPP, "/variable_importance.csv"), row.names = F)



### Copy the maxent files
# setwd(paste0(DIR, "/biomod2/"))
# copy_these <- list.files(paste0(path, "/biomod2/"), pattern = "max")
# file.copy(copy_these, modelDIR, overwrite=T)
# setwd(modelDIR)

