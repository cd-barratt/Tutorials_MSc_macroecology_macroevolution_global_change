#### 1. Initial setup and package install ###
#install.packages(c('sdm','raster','rgdal','rgbif','readr','dismo','dplyr','spThin'))
library(devtools)
install_github()
base.dir <- '/Users/cb76kecu/Dropbox (iDiv)/Micro-Macro_course/2021/Tutorial/Part_1_SDMs/'

dir.create(paste0(base.dir))
dir.create(paste0(base.dir,'env_data/'))
dir.create(paste0(base.dir,'env_data/'))
dir.create(paste0(base.dir,'env_data/present_climate'))
dir.create(paste0(base.dir,'env_data/future_climate_2070_RCP_8.5'))
dir.create(paste0(base.dir,'sdm_R/'))
dir.create(paste0(base.dir,'sdm_R/ensembles'))
dir.create(paste0(base.dir,'sdm_R/GBIF_data'))
dir.create(paste0(base.dir,'sdm_R/GBIF_data/thinned'))


#### 2. Download and prepare the environmental data ###
library(raster)
library (rgdal)

setwd <- base.dir

climate <- raster::getData('worldclim', var='bio', res=5)

template <- extent(-19, 55, -35, 38)
climate.crop <- crop(climate, template, snap="out")

for(i in 1:nlayers(climate.crop)){
  writeRaster(climate.crop[[i]], paste0(base.dir, "env_data/present_climate/bioclim_", i), "ascii", overwrite = T)
  cat('Writing bioclim',i, 
      '... \n')
}

env_data <- paste0(base.dir,'env_data/present_climate')                                 
lst <- list.files(path=env_data,pattern='asc$',full.names = T) 
climate <- stack(lst)

par(mar=c(1,1,1,1))
par(mfrow=c(4,5))
for(i in 1:nlayers(climate)){
  plot(climate[[i]], main= paste0("bio",i))
}


#### 3. Download and prepare future climate data ###
future_climate <- raster::getData('CMIP5', var='bio', res=10, rcp=85, model='AC', year=70)

future_climate.crop <- crop(future_climate, template, snap="out")

for(i in 1:nlayers(future_climate.crop)){
  writeRaster(future_climate.crop[[i]], paste0(base.dir, "env_data/future_climate_2070_RCP_8.5/bioclim_", i), "ascii", overwrite = T)
  cat('Writing bioclim',i, 
      '... \n')
}

env_data <- paste0(base.dir,'env_data/future_climate_2070_RCP_8.5/')									# lists of pred variable locations
lst <- list.files(path=env_data,pattern='asc$',full.names = T) 
future_climate <- stack(lst)

par(mar=c(1,1,1,1))
par(mfrow=c(4,5))
for(i in 1:nlayers(future_climate)){
  plot(future_climate[[i]], main= paste0("bio",i))
}

unlink('./cmip5/', recursive=TRUE)
unlink('./wc10/', recursive=TRUE)


#### 4. Obtain species data from GBIF and spatially rarefy the data ###
library(rgbif)
library(dplyr)
species_name <- "Loxodonta_cyclotis"
species_name_for_GBIF <- gsub("_", " ", species_name)
dat <- occ_search(scientificName = species_name_for_GBIF, hasCoordinate = T, limit = 5000)
dat_ne_lapply <- lapply(dat, "as.data.frame")
dat_ne <- bind_rows(dat_ne_lapply)

sp <- dat_ne[,c('decimalLongitude','decimalLatitude')]
sp$species <- 1
sp <- na.omit(sp)
write.csv(sp,'/Users/cb76kecu/Desktop/Elephant_IUCN_data/Loxodonta_cyclotis.csv')

coordinates(sp) <- ~ decimalLongitude + decimalLatitude
par(mar=c(1,1,1,1))
par(mfrow=c(1,1))
plot(climate[[1]], main = paste("Presence records"))
plot(sp, pch=21, cex=0.5, add=T)

library(spThin)
sp_to_thin <- read.csv(paste0(base.dir,'sdm_R/GBIF_data/',species_name,'.csv'))
out.dir <- paste0(base.dir,'sdm_R/GBIF_data/thinned/')

thinned_dataset <- thin( loc.data = sp_to_thin, 
                         lat.col = "decimalLatitude", long.col = "decimalLongitude", 
                         spec.col = "species", 
                         thin.par = 50, reps = 1, 
                         locs.thinned.list.return = TRUE, 
                         write.files = TRUE, 
                         max.files = 5, 
                         out.dir = out.dir, 
                         out.base = paste0(species_name), 
                         write.log.file = TRUE,
                         log.file = paste0(base.dir,'sdm_R/GBIF_data/thinned/',species_name,'_log_file.txt'))

file.rename(paste0(base.dir,'sdm_R/GBIF_data/thinned/',species_name,'_thin1.csv'), paste0(base.dir,'sdm_R/GBIF_data/thinned/',species_name,'_thinned.csv') )

sp_thinned <- read.csv(paste0(base.dir,'sdm_R/GBIF_data/thinned/',species_name,'_thinned.csv'))
coordinates(sp_thinned) <- ~ decimalLongitude + decimalLatitude

par(mar=c(1,1,1,1))
par(mfrow=c(1,2))
plot(climate[[1]], main = paste("GBIF presence records"))
plot(sp, pch=21, cex=0.5, add=T)
plot(climate[[1]], main = paste("Thinned GBIF presence records"))
plot(sp_thinned, pch=21, cex=0.5, add=T)


#### 5. Build SDMs ###
library(sdm)
#installAll()

library(usdm)
uncorrelated_vars <- extract(climate, sp)
uncorrelated_vars <- data.frame(uncorrelated_vars)

vif <- vifstep(uncorrelated_vars)
vif
climate <- exclude(climate, vif)
climate

future_climate <- exclude(future_climate, vif)
future_climate

par(mar=c(1,1,1,1))
par(mfrow=c(1,1))
plot(climate)
plot(future_climate)

setwd(paste0(base.dir,'sdm_R/ensembles/'))

sp <- read.csv(paste0(base.dir,'sdm_R/GBIF_data/thinned/',species_name,'_thinned.csv'))


library(sdm)
library(rgdal)
library(sp)

L_africana <- read.csv('/Users/cb76kecu/Desktop/Elephant_IUCN_data/Loxodonta_africana.csv')
L_africana <- L_africana[,c('decimalLongitude','decimalLatitude')]
L_africana <- na.omit(L_africana)
L_africana$presence <- 1
coordinates(L_africana) <- ~ decimalLongitude + decimalLatitude

L_cyclotis <- read.csv('/Users/cb76kecu/Desktop/Elephant_IUCN_data/Loxodonta_cyclotis.csv')
L_cyclotis <- L_cyclotis[,c('decimalLongitude','decimalLatitude')]
L_cyclotis <- na.omit(L_cyclotis)
L_cyclotis$presence <- 1
coordinates(L_cyclotis) <- ~ decimalLongitude + decimalLatitude

IUCN_shapefile_L_africana <- readOGR('/Users/cb76kecu/Desktop/Elephant_IUCN_data/L_africana/L_africana.shp')
IUCN_shapefile_L_cyclotis <- readOGR('/Users/cb76kecu/Desktop/Elephant_IUCN_data/L_cyclotis//L_cyclotis.shp')

plot(IUCN_shapefile_L_africana, lwd=0.75, col="orange")
plot(IUCN_shapefile_L_cyclotis, lwd=0.75, col="forestgreen", add=T)
plot(L_africana, pch=20, cex=0.5, col="red", add=T)
plot(L_cyclotis, pch=15, cex=0.5, col="black", add=T)

crs(L_africana) <- crs(IUCN_shapefile_L_africana)
crs(L_cyclotis) <- crs(IUCN_shapefile_L_cyclotis)

erroneous_africana <- over(L_africana, IUCN_shapefile_L_cyclotis, fn = mean) 

d <- sdmData(presence~., train=sp, predictors=climate, bg=list(n=250))

d

getmethodNames()

methods <- c('bioclim','brt','gam','rf') 
#methods <- c('bioclim','gam','glm','rpart') 

m <- sdm(presence ~., d, methods=methods, replication='sub', test.p=30, n=3) 

m

roc(m, smooth=T)
#gui(m) 

ensemble_model  <- ensemble(m,climate, setting=list(id=c(4:12), method='weighted', stat='AUC', opt=2)) 

setwd(paste0(base.dir))

colours <- colorRampPalette(c('navy','lightyellow','orange','red'))(50)
par(mar=c(1,1,1,1))
par(mfrow=c(1,1))
plot(ensemble_model, main = "Ensemble SDM - Present", col=colours, zlim=c(0,1))
points(sp_thinned, pch=19, cex=0.5, add=TRUE)
writeRaster(ensemble_model, './sdm_R/ensembles/present.asc', format="ascii", overwrite=TRUE)
file.rename('./sdm_R/ensembles/present.asc', paste0('./sdm_R/ensembles/',species_name,'_ensemble_present.asc'))

lst1 <- list.files(path=getwd(),pattern='grd$',full.names = T) 
lst2 <- list.files(path=getwd(),pattern='gri$',full.names = T) 
file.remove(lst1)
file.remove(lst2)


#### 6. Project SDMs from present to future climate ###
future_ensemble_model  <- ensemble(m,newdata=future_climate, setting=list(id=4:12, method='weighted', stat='AUC', opt=2)) 

colours <- colorRampPalette(c('navy','lightyellow','orange','red'))(50)
par(mar=c(1,1,1,1))
par(mfrow=c(1,1))
plot(future_ensemble_model, main = "Ensemble SDM - 2070, RCP8.5", col=colours, zlim=c(0,1))
points(sp_thinned, pch=19, cex=0.5, add=TRUE)
writeRaster(future_ensemble_model, './sdm_R/ensembles/future_2070_RCP_8.5.asc', format="ascii", overwrite=TRUE) 
file.rename('./sdm_R/ensembles/future_2070_RCP_8.5.asc',paste0('./sdm_R/ensembles/',species_name,'_ensemble_future_2070_RCP_8.5.asc'))

lst1 <- list.files(path=getwd(),pattern='grd$',full.names = T) 
lst2 <- list.files(path=getwd(),pattern='gri$',full.names = T) 
file.remove(lst1)
file.remove(lst2)


#### 7. Now compare present vs future models and evaluate them ###
ensemble_model <- raster(paste0('./sdm_R/ensembles/',species_name,'_ensemble_present.asc'))
future_ensemble_model <- raster(paste0('./sdm_R/ensembles/',species_name,'_ensemble_future_2070_RCP_8.5.asc'))
difference = ensemble_model-future_ensemble_model

difference_colours <- colorRampPalette(c('green','lightyellow','orange','red'))(50)
par(mar=c(5.1, 4.1, 4.1, 2.1))
par(mfrow=c(1,3))
plot(ensemble_model, main = "Ensemble SDM - Present", col=colours, zlim=c(0,1))
plot(future_ensemble_model, main = "Ensemble SDM - 2070, RCP 8.5", col=colours, zlim=c(0,1))
plot(difference, main="Difference", col=difference_colours, zlim=c(-0.3,0.7))
dev.copy(tiff, './comparison.tiff', width=2900,height=1600, res=300)
dev.off()
file.rename('./comparison.tiff', paste0('./',species_name,'_Present_vs_2070_RCP8.5_comparison.tiff'))

Australia_protected_areas <- readOGR('./Protected_areas/Australia_protected_areas.shp')
par(mar=c(5.1, 4.1, 4.1, 2.1))
par(mfrow=c(1,1))
plot(difference, main="Protected_areas", col=difference_colours, zlim=c(-0.3,0.7))
lines(Australia_protected_areas)
dev.copy(tiff, './Protected_Areas_comparison_PAs.tiff', width=2900,height=1600, res=300)
dev.off()
file.rename('./Protected_Areas_comparison_PAs.tiff', paste0('./',species_name,'_Protected_Areas_Present_vs_2070_RCP8.5_comparison.tiff'))
dev.off()


