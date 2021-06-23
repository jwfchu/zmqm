###################################################
##
## MaxEnt models for ZM/QM CSAS
## Jackson W.F. Chu - Fisheries and Oceans Canada
##
###################################################
#
# Set local  working directory
setwd("zmqm")

# Local locally installed libraries
library(dismo)                  
library(PresenceAbsence)
library(raster) 
library(usdm)
library(rJava)
library(parallel)
library(reshape2)
library(MASS)

## Initial settings

spp <- "zm"                    # species
numFolds <- 5                  # Number of cross validation folds
type <- "spatial"              # random or spatially blocked cross validation
nboot <- 100                   # Number of bootstrap samples for uncertainty calculation
calbg <- "no"                  # No to sampling background points (add code later if targeted background sampling is done)
nbg <- 10000                   # number of background sampling points (n=10,000, Liu et al)

# Categorical variables
facVars = "FW_mask"          # c('','') if more than one

# Data directories
species.dir <- "SpeciesDat"  # Location of species data (presence.csv and absence.csv)
env.dir <- "EnvLayers"       # Location of environmental rasters (.tif)

# sub directory for results
results.dir <- "Results"
dir.create(file.path(getwd(),results.dir, spp, "CV"), recursive = T)

# Working projection

NAcrs <- 
  
# LOAD DATA
  
env.list <- list.files(path=env.dir,pattern=".tif$", full.names = T)
envstack <- stack(env.list)

# remove environmental variables not used for modelling from stack

# number of env variables used in model
nvars <- length(names(modstack))

# Background point sampling
if (calbg == "yes"){
  
  pstack <- mask(modstack, pres, inverse=TRUE)                                  # Convert raster cells overlapping with presence points to NA
  bg <- sampleRandom( pstack, size=nbg, na.rm=T, xy=TRUE )                      # Randomly sample background points from pstack 
  abs <- as.data.frame(bg)
  
  coordinates(abs) <- ~x+y                                                      # Convert to spatial
  abs <- SpatialPoints(abs)
  proj4string(abs) <- NAcrs
  
  dat <- rbind(as.data.frame(pres), as.data.frame(abs))                         # create spdat from pres and abs
  dat$sp <-  c( rep(1,length(pres)), rep(0,length(abs)) )
  
  spdat <- dat[order(dat$y),]                                                   # Order by y for spatial folds
  coordinates(spdat) <- ~x+y                                                    # Convert to spatial points
  proj4string(spdat) <- NAcrs                                                   # Set coordinate reference
}

# Use FW mask 

# Use the FW mask as an environmental variable in the model
depth_mask <- as.factor(depth_mask)

## ENV variable extraction at presence-background points

envdata <- as.data.frame(extract( x=modstack, y=spdat ) )
if (facVars %in% names(envdata))  
  envdata[[facVars]] <- as.factor(envdata[[facVars]])                           # Set categorical variables as factors
obsdat <- spdat@data$sp                                                         # Convert spdat to obs vector


### CROSS VALIDATION

CV <- function( folds ){
  cuts <- cut( seq(1,length(obsdat)), breaks=folds, labels=FALSE )              # Create equally size cuts
    if( type == "spatial"){                                                     # Randomize order or spatial folds
    fold <- cuts
  } else if ( type == "random" ){
        set.seed(42)                                                           # set seed for reproducible results
    fold <- cuts[sample(length(cuts))]
  } else {
    stop( "CV type is not 'spatial' or 'random'")
  }
  cv <- list()
  for ( f in 1:numFolds ){
    foldname <- paste0("fold",f)
    cv[[foldname]][['train']] <- which(fold!=f)
    cv[[foldname]][['test']] <- which(fold==f)
  }
  
  # Plot folds
  spdat@data$cv <- fold
  foldplot <- spplot( spdat, "cv", scales=list(draw=T), cex=.5,
                      legendEntries = as.character(1:folds), 
                      key.space=list(x=0.9,y=0.9,corner=c(0,1)))
  print(foldplot)
  #return CV list
  return(cv)
}

# run CV function
cv <- CV( folds = numFolds)

## RUN MAXENT MODEL

runmaxent <- function ( fold ){
  # require
  require(dismo)
  require(rJava)
  
  # Get cv fold datasets
  train <- cv[[fold]][["train"]]
  train.envdat <- envdata[ train, ]
  train.spdat <- obsdat[ train ]
  
  # Create directory
  outdir <- file.path(results.dir,spp,"CV",paste("fold",fold, sep="_"))
  dir.create(outdir)
  # run model
  trainmodel <- maxent(x=train.envdat,
                       p=train.spdat,
                       args=c("-P","-J","writeplotdata"), 
                       path=outdir)
  # Save the model object
  save(trainmodel, file=file.path(getwd(), outdir, "maxent.model.object.RData"))
  
  # return 
  return (trainmodel)
}