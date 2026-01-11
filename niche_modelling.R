# Generate ecological niche models
# Elvira Castillo Almansa (castilloalmansae@gmail.com)

# Install and load necessary packages
library(biomod2) 
library(ecospat)  
library(raster) 
library(terra) 


# Set working directory
setwd("/route")

# Create folders to store outputs
dir.create("SpeciesXY/")
dir.create("Values/")
dir.create("Geotif/")
dir.create("Plots/")
dir.create("Images/")

# Define climate and occurrence data folders
ClimateFolder <- "/covariates" # Climate files should follow this structure: "bio_5.tif"
SpeciesFolder <- "/species"

# Load environmental variables as raster stacks
IndVar.names <- rast(list.files(ClimateFolder, full.names = TRUE))

## Generate background

 # Create mask
shapefile <- terra::vect("/shapefile.shp")
extent <- ext(c(-10, 80, 20, 80))
Mask <- crop(shapefile, extent)
plot(Mask)

 # Sample random background points.
Points.random <- spatSample(Mask, 10000)
Background.xy <- as.data.frame(geom(Points.random)[, c("x", "y")])
plot(Background.xy)

## Presence data

 # Define species file name.
myRespName <- "species.csv" 
SpeciesName.2 <- gsub(".csv","", myRespName, fixed=TRUE)

 # Load, clean and save occurrence data
Species.data <- read.csv(paste(SpeciesFolder,myRespName,sep=""), sep = ";")
head(Species.data)
Species.data$Longitude <- as.numeric(gsub(",", ".", Species.data$Longitude))
Species.data$Latitude <- as.numeric(gsub(",", ".", Species.data$Latitude))
Species.data<- Species.data[c("Longitude", "Latitude")]

Species.data.temp <- extract(IndVar.names,Species.data)
Species.data.temp2 <- cbind(Species.data, Species.data.temp)
Species.data.temp3  <- na.omit(Species.data.temp2)
XY <- Species.data.temp3[1:2]
colnames(XY) <- c("x", "y")
head(XY)
write.csv(XY, paste0("SpeciesXY/",SpeciesName.2,".csv"), row.names = F)

 # Prepare data for biomod2.
myResp.xy <- rbind(XY,Background.xy);names(myResp.xy)<-c("x","y");row.names(myResp.xy)<-c(1:nrow(myResp.xy)) 
myResp <- data.frame(c(rep(1,nrow(XY)),rep(NA,nrow(Background.xy))));names(myResp)<-"pa"; row.names(myResp)<-c(1:nrow(myResp.xy)) 
myExpl <- data.frame (extract (IndVar.names, myResp.xy))

## Modelling for current scenarios

 # Format input for biomod2.
myBiomodData <- BIOMOD_FormatingData (resp.var = myResp, resp.xy = myResp.xy, expl.var = IndVar.names, resp.name = SpeciesName.2, PA.nb.rep = 1, PA.nb.absences = nrow(Background.xy), PA.strategy = 'random') 

 # Run models (10 reps for 3 algorithms).
myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData, modeling.id = 'AllModels', models = c('GBM','RF', 'GLM'), CV.strategy = "random", CV.nb.rep=10, CV.perc=0.8, weights=NULL, var.import=3, metric.eval = c('ROC','TSS'), scale.models = FALSE, do.progress = TRUE, prevalence = 0.5, seed.val = 42)

 # Build ensemble using models with ROC > 0.8.
myBiomodEM.ROC  <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                           models.chosen = 'all',
                                           em.by = 'all',
                                           em.algo = c("EMwmean"),
                                           metric.select = c('ROC'),
                                           metric.select.thresh = c(0.8),
                                           var.import = 3,
                                           metric.eval = c('ROC', "TSS", "KAPPA"),
                                           seed.val = 42)

 # Project models to study area.
myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                  new.env = IndVar.names,
                                  proj.name = "Current",
                                  models.chosen = 'all')

# Build ensemble forecast.
BIOMOD_EnsembleForecasting(bm.em = myBiomodEM.ROC, 
                           bm.proj = myBiomodProj,
                           models.chosen = 'all',
                           metric.binary = 'all',
                           metric.filter = 'all') 

# Export raster
Pred.ROC <- rast(paste(SpeciesName.2, "/proj_Current/proj_Current_", SpeciesName.2, "_ensemble.tif", sep=""))
terra::writeRaster (Pred.ROC, paste("Geotif/",SpeciesName.2,".Current.tif",sep=""), overwrite=TRUE)

# Compute evaluation metrics (mean & SD)

## ROC
myEMeval.replicas.ROC <- get_evaluations(myBiomodModelOut, metric.eval= "ROC")
myEMeval.replicas.ROC.2  <- myEMeval.replicas.ROC[1:30, 10]
myEMeval.replicas.ROC.2.mean <- round(mean(myEMeval.replicas.ROC.2), digits = 3)
myEMeval.replicas.ROC.2.sd <- round(sd(myEMeval.replicas.ROC.2), digits = 3)

## TSS
myEMeval.replicas.TSS <- get_evaluations(myBiomodModelOut, metric.eval= "TSS")
myEMeval.replicas.TSS.2  <- myEMeval.replicas.TSS[1:30, 10]
myEMeval.replicas.TSS.2.mean <- round(mean(myEMeval.replicas.TSS.2), digits = 3)
myEMeval.replicas.TSS.2.sd <- round(sd(myEMeval.replicas.TSS.2), digits = 3)

# Consensus model evaluation.
myEMeval.Consenso <- get_evaluations(myBiomodEM.ROC)

# Save evalutation results. 
myEMeval.replicas <- get_evaluations(myBiomodModelOut)
save(myEMeval.replicas,file=paste("Values/",SpeciesName.2,"_replicas.RData",sep=""))
write.table(myEMeval.replicas,file=paste("Values/",SpeciesName.2,"_replicas.txt",sep=""),row.names=T,col.names=T,sep=",")
save(myEMeval.Consenso,file=paste("Values/",SpeciesName.2,"_consenso.RData",sep=""))

# Plot and export map.
pdf(paste("Plots/",SpeciesName.2,"_Current.pdf",sep=""),width=11.7, height=8.3)
plot(Pred.ROC, legend = TRUE, col = rev(terrain.colors(50)))
points(XY,pch=20,bg="black",cex=1)
title(main = paste(SpeciesName.2, "   Current", sep = ""), line = 3)
mtext(side = 1, line = 2, text = paste("Ocurrences = ", nrow(XY),"  | ROC (Validation mean)= ", myEMeval.replicas.ROC.2.mean, sep = ""))
par(cex=.7);
dev.off()

## Save full session.
aaa <- paste("Images/",SpeciesName.2,".RData" ,sep=""); save.image (file = aaa)
