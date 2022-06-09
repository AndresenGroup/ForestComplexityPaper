# by Jacob May
# Last updated 6/8/22
# last tested with R 4.1.2 and lidR 4.0.1

# This file runs the metrics from allmetrics.R on a list of lidar datasets
# the lidar datasets must be already cleaned and normalized
# I used lidR and Cloud Compare to do this, see their documentation if needed

# obviously, file paths and naming conventions need to be changed
# the settings file I used will be included but could be substituted for
# just a list of the different lidar dataset names

###############################################################################

#libraries
library(lidR)
library(future)
library(data.table)
library(ggplot2)
library(tidyverse)

library(sf)
library(sp)
library(raster)
library(rgdal) 

#colors
library(RColorBrewer)
library(wesanderson)
library(viridis)
library(ggsci)


#future options
plan(multisession)
options(future.rng.onMisue = "ignore")
# set to number of cores on your pc, doesn't matter for most functions
set_lidr_threads(8) 

#script for getting metrics for all sites
start <- Sys.time()

src <- "C:/Users/andresenlab/Documents/Jake/scripts"
setwd(src)
source("allmetrics.R")


projectDir <- "C:/Users/andresenlab/Documents/Jake/sites/"
setwd(projectDir)
settings <- read.csv(file = "metrics_settings.csv")
df_data <- data.frame()
gridRes <- 100

for(i in 1:nrow(settings)){
  
  #i <- 1
  
  #Load the normalized and cropped las file
  plotName <- settings[i, 1]
  print(plotName)
  setwd(file.path(projectDir, plotName))
  las <- readLAS(paste(plotName, "_crop.laz", sep = ""), select = "xyzrn")
  
  # ended up not using footprint masks so this section was unused
  
  #check for footprint extent overlap issues
  # if(settings[i, 4] == FALSE){
  #   foot <- raster(paste(plotName, "_footRast.tif", sep = ""))
  #   area <- sum(footCropBool1@data@values, na.rm = TRUE)
  #   hmax <- grid_metrics(las, ~max(Z), gridRes)
  #   footCropBool <- resample(footCropBool1, hmax, method = "bilinear")
  # } else {
  #   foot <- NULL
  #   area <- (las@bbox[1, 2] - las@bbox[1, 1])*(las@bbox[2, 2] - las@bbox[2, 1])
  # }
  
  
  
  src <- "C:/Users/andresenlab/Documents/Jake/scripts"
  setwd(src)
  source("allmetrics.R")
  
  plot(hmax(las, gridRes, rasterResult = TRUE), col = height.colors(200))
  
  # compute the metrics from allmetrics.R
  rumpleMetric <- rumple(las, gridRes)
  print("rumple")
  verticalDistMax_result <- verticalDistMax(las, gridRes)
  print("vdistmax")
  hmax_result <- hmax(las, gridRes)
  print("hmax")
  hvar_result <- hvar(las, gridRes)
  print("hvar")
  hmean_result <- hmean(las, gridRes, groundCut = 0.5)
  print("hmean")
  density <- density.metric(las, gridRes)
  print("density")
  gap <- gap.fraction(las, gridRes, groundCut = 0.25)
  print("gap frac")
  shannon <- shannon.index(las, gridRes, groundCut = 0.05)
  print("shannon")
  VCI_result <- VCImetric(las, gridRes, groundCut = 0.05)
  print("VCI")
  LAI_result <- LAI(las, gridRes, dz = gridRes)
  print("LAI")
  RH <- Relative.Height(las, gridRes)
  print("Relative Heights")
  
  # add the metrics for the dataset to a data frame
  df_site <- data.frame(
    site = plotName,
    #area = area,
    rumple = rumpleMetric,
    verticalDistMax = verticalDistMax_result,
    maxZ_mean = hmax_result[1],
    maxZ_sd = hmax_result[2],
    sdZ_mean = hvar_result[1],
    sdZ_sd = hvar_result[2],
    meanZ_mean = hmean_result[1],
    meanZ_sd = hmean_result[2],
    density_mean = density[1],
    density_sd = density[2],
    gap_fraction = gap,
    shannon_mean = shannon[1],
    shannon_sd = shannon[2],
    VCI_mean = VCI_result[1],
    VCI_sd = VCI_result[2],
    LAI_mean = LAI_result[1],
    LAI_sd = LAI_result[2],
    RH25 <- RH[1],
    RH50 <- RH[2],
    RH75 <- RH[3],
    RH95 <- RH[4],
    canopy_ratio <- RH[5]
    )
  
  df_site
  
  #combine data frame with previous datasets'
  df_data <- rbind(df_data, df_site)
  
  #gc()
  
}
# save the data
write.csv(df_data, "metrics_3_8_100m.csv", row.names = FALSE)
totaltime <- Sys.time() - start
print(totaltime)
