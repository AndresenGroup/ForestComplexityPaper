# by Jacob May
# last updated 6/8/22

# this file is basically a library of all the forest canopy structural metrics
# used in the paper Unraveling Forest Complexity: Resource Use Efficiency, 
# Disturbance, and the Structure-Function Relationship
# https://doi.org/10.1029/2021JG006748
# also includes some utility functions

# last tested with R 4.1.2 and lidR 4.0.1
# references are induced in the above paper

###############################################################################

# other notes:
# some functions contain footRaster input,
# we tested using a raster mask instead of cropping in some cases
# you can use a "footprint" or mask raster where 1 is a pixel you want and 0
# for pixels you don't want
# however this wasn't used in the final analysis and is only coded for some
# of the functions

##############################################################################

#fill NA pixels in a raster with the average of surrounding pixels
#if there are no surrounding pixels then it loops until the surroundings
#are filled in, which makes it not good for rasters missing large areas
#input: raster
#output: raster

na.fill <- function(raster){
  
  while(any(is.na(raster[]))){
    for(i in 1:length(raster@data@values)){
      
      
      if(is.na(raster[i])){
        
        sur <- vector()
        j <- 1
        
        if(i > raster@ncols){
          sur[j] <- raster[i - raster@ncols]
          j <- j + 1
        }
        if(i > 1){
          if((i - 1) %% raster@ncols != 0){
            sur[j] <- raster[i - 1]
            j <- j + 1
          }
        }
        if(i %% raster@ncols !=0){
          sur[j] <- raster[i + 1]
          j <- j + 1
        }
        if(i <= length(raster) - raster@ncols){
          sur[j] <- raster[i + raster@ncols]
          j <- j + 1
        }
        
        raster[i] <- mean(na.exclude(sur))
        
        
      }
    }
  }
  
  return(raster)
}

################################################################

#rumple metric
#simple code to find rumple
#input: las, gridRes = resolution for CHM, plotting = make a plot
#output: single number, the rumple index

rumple <- function(las, gridRes = 1, plotting = FALSE){
  
  chm <- grid_canopy(las, res = gridRes, pitfree())
  chm <- na.fill(chm)
  if(plotting){plot(chm, col = magma(100))}
  
  return(rumple_index(chm))
}

###############################################################

#vertical point distribution
#get the vertical distribution of point in a point cloud
#input: las, dz = vertical bin size, gndCutoff = understory removal height,
#maxZ = maximum height - if not specified will take max of cloud
#so if running on multiple sites this should be set so they are the same
#output: data.frame with first column bin height and second column fraction of points in the bin

verticalDist <- function(las, dz = 1, gndCutoff = 0, maxZ = NULL){
  
  if(is.null(maxZ)){
    maxZ <- max(las@data[["Z"]])
  }
  
  pointCount <- vector()
  height <- vector()
  count <- 1

  for(i in seq(gndCutoff, maxZ, dz)){
    pointCount[count] <- length(filter_poi(las, Z >= i & Z < i + dz)@data[["Z"]])
    height[count] <- i
    count <- count + 1
  }
  
  pointCount <- pointCount/max(pointCount)
  return(data.frame(height = height, pointCountNorm = pointCount))
}

######################################################################

#vertical distribution maximum metric
#just gives the height bin of maximum points in the vertical distribution
#inputs: same as previous, gndCutoff needs to be set high though
#output: single value, height of max point density
#if both this and the distribution data frame are needed just use the data frame one
#to reduce compute time

verticalDistMax <- function(las, dz = 1, gndCutoff = 5, maxZ = NULL){
  
  if(is.null(maxZ)){
    maxZ <- max(las@data[["Z"]])
  }
  
  pointCount <- vector()
  height <- vector()
  count <- 1
  
  for(i in seq(gndCutoff, maxZ, dz)){
    pointCount[count] <- length(filter_poi(las, Z >= i & Z < i + dz)@data[["Z"]])
    height[count] <- i
    count <- count + 1
  }
  
  return(height[which.max(pointCount)])
}


verticalDistGridHelper <- function(Z){
  
  if(!is.null(Z)){
    test <- hist(Z, breaks = c(min(Z), 50, 1), plot = FALSE)
    maxDen <- test$mids[which.max(test$counts)]
    
  } else {
    maxDen <- NA
  }
  
  
  return(maxDen)
  
}


verticalDistGrid <- function(las, gridRes, gndCutoff = 0, rasterResut = FALSE){
  
  
  lasCut <- filter_poi(las, Z >= gndCutoff)
  VDGRast <- grid_metrics(lasCut, ~verticalDistGridHelper(Z), gridRes)
  
  if(rasterResut){
    return(VDGRast)
  } else {
    VDM_sd <- sd(na.exclude(VDGRast@data@values))
    VDM_mean <- mean(na.exclude(VDGRast@data@values))
    return(c(VDM_mean, VDM_sd))
  }
  
}




########################################################################

#max height metrics
#input: las, resolution, climatology footprint to crop with, raster result or numeric result
#output: raster or vector with [1]=mean and [2]=sd

cr_rasthmax <- function(las, gridRes = 1, footRaster = NULL, rasterResult = FALSE){
  
  maxRast <- grid_metrics(las, ~max(Z), gridRes)
  if(!is.null(footRaster)){
    maxRast <- maxRast*footRaster
    maxRast[maxRast == 0] <- NA
  }
  
  if(rasterResult){
    return(maxRast)
  } else {
    maxZ_sd <- sd(na.exclude(maxRast@data@values))
    maxZ_mean <- mean(na.exclude(maxRast@data@values))
    return(c(maxZ_mean, maxZ_sd))
  }
}

#######################################################################

#height variance metrics
#same as previous with sd instead of max
#includes a ground cutoff since result would be heavily effected by ground points

hvar <- function(las, gridRes = 1, groundCut = 0.5, footRaster = NULL, rasterResult = FALSE){
  
  varRast <- grid_metrics(filter_poi(las, Z > groundCut), ~sd(Z), gridRes)
  if(!is.null(footRaster)){
    varRast <- varRast*footRaster
    varRast[varRast == 0] <- NA
  }
  
  if(rasterResult){
    return(varRast)
  } else {
    sdZ_sd <- sd(na.exclude(varRast@data@values))
    sdZ_mean <- mean(na.exclude(varRast@data@values))
    return(c(sdZ_mean, sdZ_sd))
  }
}

########################################################################

hmean <- function(las, gridRes = 1, groundCut = 0.5, footRaster = NULL, rasterResult = FALSE){
  
  meanRast <- grid_metrics(filter_poi(las, Z > groundCut), ~mean(Z), gridRes)
  if(!is.null(footRaster)){
    meanRast <- meanRast*footRaster
    meanRast[meanRast == 0] <- NA
  }
  
  if(rasterResult){
    return(meanRast)
  } else {
    meanZ_sd <- sd(na.exclude(meanRast@data@values))
    meanZ_mean <- mean(na.exclude(meanRast@data@values))
    return(c(meanZ_mean, meanZ_sd))
  }
}

##########################################################################
#density metrics
#density mean is mostly useless for characterization but included for completeness

density.metric <- function(las, gridRes = 1, footRaster = NULL, rasterResult = FALSE){
  
  grid_density <- grid_metrics(las, ~length(Z)/(gridRes*gridRes), gridRes)
  if(!is.null(footRaster)){
    grid_density <- grid_density*footRaster
    grid_density[grid_density == 0] <- NA
  }
  
  if(rasterResult){
    return(grid_density)
  } else {
    density_mean <- mean(na.exclude(grid_density@data@values))
    density_sd <- sd(na.exclude(grid_density@data@values))
    return(c(density_mean, density_sd))
  }
}

##########################################################################
#Gap Fraction, or deep gap fraction or 1 - cover fraction
#can't be cropped the way the other raster functions are cause it relies on counting empty pixels

gap.fraction <- function(las, gridRes = 1, groundCut = NULL, footRaster = NULL, rasterResult = FALSE){
  
  if(is.null(groundCut)){
    groundCut <- gridRes
  }
  
  ground_density <- grid_metrics(filter_poi(las, Z < groundCut), ~length(Z), gridRes)
  
  if(!is.null(footRaster)){
    ground_density[is.na(ground_density[])] <- Inf
    ground_density <- ground_density*footRaster
    ground_density[ground_density == 0] <- NA
    ground_density[is.infinite(ground_density[])] <- 0
    gapFraction <- length(ground_density[ground_density])/length(na.exclude(ground_density@data@values))
  } else {
    gapFraction <- length(na.exclude(ground_density@data@values))/length(ground_density@data@values)
  }
  
  if(rasterResult){
    return(entropyRast)
  } else {
    return(gapFraction)
  }
}

##########################################################################

#normalized Shannon's vertical complexity index

shannon.index <- function(las, 
                          gridRes = 1, 
                          dz = NULL, 
                          groundCut = 0, 
                          footRaster = NULL, 
                          rasterResult = FALSE){
  
  if(is.null(dz)){
    dz <- gridRes
  }
  
  entropyRast <- grid_metrics(filter_poi(las, Z > groundCut), ~entropy(Z, by = gridRes), gridRes)
  
  if(!is.null(footRaster)){
    entropyRast <- entropyRast*footRaster
    entropyRast[entropyRast == 0] <- NA
  }
  
  if(rasterResult){
    return(entropyRast)
  } else {
    entropy_sd <- sd(na.exclude(entropyRast@data@values))
    entropy_mean <- mean(na.exclude(entropyRast@data@values))
    return(c(entropy_mean, entropy_sd))
  }
}

######################################################################
#Vertical complexity index (VCI), how does the vertical dist differ from random?

VCImetric <- function(las, 
                gridRes = 1, 
                dz = 1, 
                groundCut = 0.05, 
                footRaster = NULL, 
                rasterResult = FALSE){
  
  if(is.null(dz)){
    dz <- gridRes
  }
  
  maxZ <- max(las@data[["Z"]])
  VCIRast <- grid_metrics(filter_poi(las, Z > groundCut), ~VCI_FIX(Z, by = 0.1, zmax = 100), gridRes)
  #VCIRast <- grid_metrics(filter_poi(las, Z > groundCut), ~VCI(Z, by = dz, zmax = grid_metrics(las, ~max(Z),)), gridRes)
  
  
  if(!is.null(footRaster)){
    VCIRast <- VCIRast*footRaster
    VCIRast[VCIRast == 0] <- NA
  }
  
  if(rasterResult){
    return(VCIRast)
  } else {
    VCI_sd <- sd(na.exclude(VCIRast@data@values))
    VCI_mean <- mean(na.exclude(VCIRast@data@values))
    return(c(VCI_mean, VCI_sd))
  }
}

#########################################
#VCI was not working so here it is copied from
#the lidr package and renamed to fix function name errors
VCI_FIX = function(z, zmax, by = 1)
{
  z <- z[z < zmax]
  return(entropy(z, by, zmax))
}


#########################################################################
#relative height function, can get the height of specified quartiles of cloud
#or used for each pixel in a raster
#easy to add more quartiles here
#ground cut is very important, should be small but still important

Relative.Height.Helper <- function(Z){
  
  zSort <- sort(Z)
  RH25 <- zSort[ceiling(length(Z)*0.25)]
  RH50 <- zSort[ceiling(length(Z)*0.50)]
  RH75 <- zSort[ceiling(length(Z)*0.75)]
  #if(ceiling(length(Z)*0.95) > length(Z)){
  #  RH95 <- zSort[length(Z)]
  #} else {
  RH95 <- zSort[ceiling(length(Z)*0.95)]
  #}

  
  list(RH25 = RH25, RH50 = RH50, RH75 = RH75, RH95 = RH95)
}

##########################################################################
#Relative height rasters, only using mean to declutter
#requires "Relative.Height.Helper" above
#has major issues with low grid resolution, 1m max recommended for now


Relative.Height <- function(las,
                            gridRes = 1,
                            groundCut = 0,
                            footRaster = NULL,
                            rasterResult = FALSE){
  
  RHrast <- grid_metrics(filter_poi(las, Z > groundCut),
                         ~Relative.Height.Helper(Z),
                         gridRes)
  
  CRrast <- (subset(RHrast, 4) - subset(RHrast, 1))/subset(RHrast, 4)
  
  if(rasterResult){
    plot(CRrast, col = viridis(100))
    plot(RHrast, col = height.colors(100))
    #return(RHrast)
  } else {
    RH25_mean <- mean(subset(RHrast, 1)@data@values, na.rm = TRUE)
    RH50_mean <- mean(subset(RHrast, 2)@data@values, na.rm = TRUE)
    RH75_mean <- mean(subset(RHrast, 3)@data@values, na.rm = TRUE)
    RH95_mean <- mean(subset(RHrast, 4)@data@values, na.rm = TRUE)
    CR_mean <- mean(CRrast@data@values, na.rm = TRUE)
    
    return(c(RH25_mean, RH50_mean, RH75_mean, RH95_mean, CR_mean))
  }
  
}

########################################################################
#Canopy Ratio

# just added as fifth return from relative heights, saves computation


#########################################################################
#LAI from forestR method, for individual pixels
#corrected to fit LAI2000 field measurement sites from Ting Zheng

LAI.Helper <- function(Z, 
                dZ = 1, 
                gndCut = 8, 
                VAImax = 8){
  
  LAD <- vector()
  b <- 1 - exp(-0.5*VAImax)
  returns <- length(Z)
  
  if(ceiling(max(Z)) < gndCut + dZ){
    return(0)
  }
  
  VAI_obs <- -0.5*log(1 - (sum(Z > gndCut)/returns)*b)
  
  for(i in seq(ceiling(max(Z)), gndCut + dZ, dZ*-1)){
    if(sum(Z < i-dZ) == 0) {
      LAD[i] <- 0
    } else {
      LAD[i] <- VAI_obs * sum(Z < i & Z > i-dZ)/returns
    }
  }
  #return(sum(na.exclude(LAD)))
  #based on linear fit from Ting's field sites
  VAI_adj <- 1.7598*VAI_obs + 0.4704 #old values from forestr cor, using for baileys paper
  #VAI_adj <- VAI_obs
  #commented out adjustment to do correlation testing
  # VAI_adj <- 0.5252*VAI_obs + 1.2616 #new values from almeida, using for ting/general use
  
  return(VAI_adj)
}

#########################################################################
#LAI metrics
#BROKEN FOR FOOTPRINTS

LAI <- function(las, 
                gridRes = 1, 
                dz = 1, 
                groundCutoff = 1, 
                VAI_max = 8,
                footRaster = NULL, 
                rasterResult = FALSE){
  
  LAIRast <- grid_metrics(las, ~LAI.Helper(Z), gridRes)
  
  if(!is.null(footRaster)){
    #wont work for LAI, need to use the inf method
    LAIRast <- LAIRast*footRaster
    LAIRast[LAIRast == 0] <- NA
  }
  
  if(rasterResult){
    return(LAIRast)
  } else {
    LAI_sd <- sd(na.exclude(LAIRast@data@values))
    LAI_mean <- mean(na.exclude(LAIRast@data@values))
    return(c(LAI_mean, LAI_sd))
  }
}

#######################################################################
#####ignore for now
LAD.Helper <- function(Z, 
                       dZ = 1, 
                       gndCut = 1, 
                       VAImax = 8){
  
  LAD <- vector()
  height <- vector()
  b <- 1 - exp(-0.5*VAImax)
  returns <- length(Z)
  
  if(ceiling(max(Z)) < gndCut + dZ){
    return(0)
  }
  
  VAI_obs <- -0.5*log(1 - (sum(Z > gndCut)/returns)*b)
  
  #VAI_adj <- 1.7598*VAI_obs + 0.4704
  VAI_adj <- VAI_obs
  
  count <- 1
  for(i in seq(ceiling(max(Z)), gndCut + dZ, dZ*-1)){
    if(sum(Z < i-dZ) == 0) {
      LAD[i] <- 0
    } else {
      LAD[i] <- VAI_adj * sum(Z < i & Z > i-dZ)/returns
    }
    height[count] <- i
    count <- count + 1
  }
  return(sum(na.exclude(LAD)))
  #based on linear fit from Ting's field sites
  #VAI_adj <- 1.7598*VAI_obs + 0.4704
  
  #return(height[which.max(LAD)])
  return(VAI_adj)
}

#commented out the vai_adj lines to do lai testing on correlation with lai2000

###############################################################################
# LAD max was supposed to return the height of the max leaf area density
# in a pixel, this may not be working correctly and was not used in the paper

LAD.max <- function(las, 
                gridRes = 1, 
                dz = 1, 
                groundCutoff = 6, 
                VAI_max = 8,
                footRaster = NULL, 
                rasterResult = FALSE){
  
  LADmaxRast <- grid_metrics(filter_poi(las, Z >= groundCutoff), ~LAD.Helper(Z), gridRes)
  
  if(!is.null(footRaster)){
    #wont work for LAI, need to use the inf method
    LADmaxRast <- LADmaxRast*footRaster
    LADmaxRast[LADmaxRast == 0] <- NA
  }
  
  if(rasterResult){
    return(LADmaxRast)
  } else {
    LAD_sd <- sd(na.exclude(LADmaxRast@data@values))
    LAD_mean <- mean(na.exclude(LADmaxRast@data@values))
    return(c(LAD_mean, LAD_sd))
  }
}


