########################
## Metric calibration ##
########################

setwd("C:/Users/sarah.bauduin/Documents/GitHub/appendix_wolfOccupancy")

library(raster)
library(rgeos)
library(colorRamps)
library(plot.matrix)
library(ggplot2)
library(viridis)

# Two parameters to estimate:
# The threshold for the occupancy probability above which a cell is considered "occupied"
# The buffer around an "occupied" cell to consider surrounding cells also "occupied"

# To calibrate these two parameters, we aim to reproduce three patterns:
# Having "occupied" cells inside the pack permanent presence areas (PPAs)
# Having "occupied" cells inside the non-pack PPAs
# Having non-occupied cells outside the PPAs

# An area is defined as a PPA if it is occupied two consecutive winters
# A PPA is a circle of radius 9km (represents a pack territory)
load("data/PPA.RData") # listPPA, listPackPPA, listNonPackPPA
# PPA data are available from winter 2009/2010 to winter 2019/2020
yearStartPPA <- 2009


##########################
# Occupancy maps and PPA #
##########################

# Load the estimated occupancy maps
load("modelOutputs/noScaleShortDisp_noCull_effAl_7500_rasterYears.RData")
# The first map is of winter 1993/1994
load("data/franceShape.RData")

for(i in 1:length(rasterYears[17:27])){ # from winter 2009/2010 to winter 2019/2020
  
  plot(rasterYears[17:27][[i]], main = paste0("Winter ", yearStartPPA + i - 1, "-", yearStartPPA + i))
  plot(franceShape, add = TRUE)
  plot(listPackPPA[[i]], add = TRUE, lwd = 2, border = "red")
  plot(listNonPackPPA[[i]], add = TRUE, lwd = 2, border = "black")
  legend("topleft", legend = c("Pack PPA", "Non-pack PPA"), fill = c("red", "black"))
}

# Function to create PPA
# Combine two consecutive occupancy maps 
# Depends on the two parameters to adjust
occ2winters <- function(occWinter1, occWinter2, probThresh, bufferkm){
  
  # Maps in rasterYears do not have coordinate system attached
  occWinter1@crs <- CRS("+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
  occWinter2@crs <- CRS("+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
  
  occWinterBoth <- occWinter1
  occWinterBoth[] <- 0
  
  # Raster 1st winter - Define the "occupied" cells
  occWinter1[which(is.na(raster::values(occWinter1)) | raster::values(occWinter1) < probThresh)] <- 0
  occWinter1[which(raster::values(occWinter1) >= probThresh)] <- 1
  occWinter1[raster::values(occWinter1) == 0] <- NA
  if(bufferkm != 0){
    occWinter1 <- buffer(occWinter1, width = bufferkm * 1000)
  }
  # Raster 2nd winter - Define the "occupied" cells
  occWinter2[which(is.na(raster::values(occWinter2)) | raster::values(occWinter2) < probThresh)] <- 0
  occWinter2[which(raster::values(occWinter2) >= probThresh)] <- 1
  occWinter2[raster::values(occWinter2) == 0] <- NA
  if(bufferkm != 0){
    occWinter2 <- buffer(occWinter2, width = bufferkm * 1000)
  }
  
  # Combine the two maps
  occWinterBoth[] <- raster::values(occWinter1) + raster::values(occWinter2)
  occWinterBoth[is.na(raster::values(occWinterBoth)) | raster::values(occWinterBoth) != 2] <- 0
  occWinterBoth[raster::values(occWinterBoth) == 2] <- 1
  return(occWinterBoth)
}


# Function to compare the raster maps made of two consecutive winters with the PPAs
# Only for the cells that have been sampled
validationPPA <- function(raster2winter, PPA, packPPA, nonPackPPA, sampledCells){
  
  raster2winter@crs <- CRS("+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
  
  # Keep only the data in the sampled areas
  raster2winter <- crop(raster2winter, sampledCells)
  PPA <- crop(PPA, sampledCells)
  packPPA <- crop(packPPA, sampledCells)
  nonPackPPA <- crop(nonPackPPA, sampledCells)
  
  # Proportion of "occupied" cells inside the pack PPAs polygons
  inPackPPA <- unlist(extract(raster2winter, packPPA))
  propInPackPPA <- sum(inPackPPA, na.rm = TRUE) / length(inPackPPA)
  
  # Proportion of "occupied" cells inside the non-pack PPAs polygons
  inNonPackPPA <- unlist(extract(raster2winter, nonPackPPA))
  propInNonPackPPA <- sum(inNonPackPPA, na.rm = TRUE) / length(inNonPackPPA)
  
  # Proportion of "non occupied" cells outside the PPAs polygons
  inPPA <- do.call("rbind", extract(raster2winter, PPA, cellnumbers = TRUE))
  outPPA <- setdiff(1:length(raster2winter), unique(inPPA[,1]))
  propOutPPA <- sum(raster2winter[outPPA] == 0) / length(outPPA)
  
  return(list(propInPackPPA, propInNonPackPPA, propOutPPA))
}


# Calibrating parameters
# Computing the three patterns with probThresh from 0 to 1 (each 0.01)
# and for the bufferkm with 0, 10, 15 km

# Define for each year the cells that were sampled
cellSampled <- list()
for(i in 1:length(rasterYears)){
  rasterOfTheYear <- rasterYears[[i]]
  rasterOfTheYear[!is.na(rasterOfTheYear)] <- 1
  cellSampled[[i]] <- rasterToPolygons(rasterOfTheYear, dissolve = TRUE)
  cellSampled[[i]]@proj4string <- CRS("+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
}

rasterYearCal <- rasterYears[16:27] # Winter 2008/2009 to winter 2019/2020
cellSampledCal <- cellSampled[16:27]

listResPropInPackPPA <- list()
listResPropInNonPackPPA <- list()
listResPropOutPPA <- list()

for(year in 1:(length(rasterYearCal) - 1)){
  
  ResPropInPackPPA <- matrix(nrow = 3, ncol = length(seq(0.01, 1, by = 0.01)))
  ResPropInNonPackPPA <- matrix(nrow = 3, ncol = length(seq(0.01, 1, by = 0.01)))
  ResPropOutPPA <- matrix(nrow = 3, ncol = length(seq(0.01, 1, by = 0.01)))
  
  j <- 1
  for(probThresh in seq(0.01, 1, by = 0.01)){
    i <- 1
    for(bufferkm in c(0, 10, 15)){
      
      # Build a map of PPA with two consecutive winters
      raster2winter <- occ2winters(occWinter1 = rasterYearCal[[year]], 
                                   occWinter2 = rasterYearCal[[year + 1]],
                                   probThresh = probThresh, bufferkm = bufferkm)
      
      # Compare the created PPA with the ones from field data to compute the three patterns
      valid2winter <- validationPPA(raster2winter = raster2winter, PPA = listPPA[[year]], 
                                    packPPA = listPackPPA[[year]], nonPackPPA = listNonPackPPA[[year]],
                                    sampledCells = cellSampledCal[[year]])
      
      ResPropInPackPPA[i,j] <- valid2winter[[1]]
      ResPropInNonPackPPA[i,j] <- valid2winter[[2]]
      ResPropOutPPA[i,j] <- valid2winter[[3]]
      
      i <- i+1
      print(paste0("bufferkm  = ", bufferkm, " km"))
    }
    print(paste0("probThresh = ", probThresh))
    j <- j+1
  }
  
  listResPropInPackPPA[[year]] <- ResPropInPackPPA 
  listResPropInNonPackPPA[[year]] <- ResPropInNonPackPPA
  listResPropOutPPA[[year]] <- ResPropOutPPA
  
  print(paste0("year = ", year))
}

save(listResPropInPackPPA, listResPropInNonPackPPA, listResPropOutPPA, 
     rasterYearCal, cellSampledCal,
     file = "calibrationOutputs/noScaleShortDisp_noCull_effAl_7500_validPPA.RData")


########################
# Calibrate parameters #
########################

load("calibrationOutputs/validPPA.RData")

# Compute the mean pattern values over all years
meanPropInPackPPA <- apply(simplify2array(listResPropInPackPPA), 1:2, mean)
meanPropInNonPackPPA <- apply(simplify2array(listResPropInNonPackPPA), 1:2, mean)
meanPropOutPPA <- apply(simplify2array(listResPropOutPPA), 1:2, mean)

# Rescale between 0 and 1 with 0 the worst combination and 1 the best combination
meanPropInPackPPARescale = (meanPropInPackPPA - min(meanPropInPackPPA)) / (max(meanPropInPackPPA) - min(meanPropInPackPPA))
meanPropInNonPackPPARescale = (meanPropInNonPackPPA - min(meanPropInNonPackPPA)) / (max(meanPropInNonPackPPA) - min(meanPropInNonPackPPA))
meanPropOutPPARescale = (meanPropOutPPA - min(meanPropOutPPA)) / (max(meanPropOutPPA) - min(meanPropOutPPA))

# Reorder the lines of the matrix to have results for the highest buffer first
meanPropInPackPPARescale2 <- rbind(meanPropInPackPPARescale[3,], meanPropInPackPPARescale[2,], meanPropInPackPPARescale[1,])
meanPropInNonPackPPARescale2 <- rbind(meanPropInNonPackPPARescale[3,], meanPropInNonPackPPARescale[2,], meanPropInNonPackPPARescale[1,])
meanPropOutPPARescale2 <- rbind(meanPropOutPPARescale[3,], meanPropOutPPARescale[2,], meanPropOutPPARescale[1,])

plot(meanPropInPackPPARescale2, border = NA, col = viridis(20),
     key = list(side = 3, cex.axis = 0.75), xlab = "Occupancy probability threshold ", 
     ylab = "Buffer size around 'occupied' cells", axis.col = NULL, axis.row = NULL,
     #main = "Proportion of 'occupied' cells in pack PPAs")
     main = "")
axis(1, at = 1:ncol(meanPropInPackPPARescale2), labels = seq(0.01, 1, by = 0.01))
axis(2, at = 1:nrow(meanPropInPackPPARescale2), labels = c(0, 10, 15))

plot(meanPropInNonPackPPARescale2, border = NA, col = viridis(20),
     key = list(side = 3, cex.axis = 0.75), xlab = "Occupancy probability threshold ", 
     ylab = "Buffer size around 'occupied' cells", axis.col = NULL, axis.row = NULL,
     #main = "Proportion of 'occupied' cells in non-pack PPAs")
     main = "")
axis(1, at = 1:ncol(meanPropInNonPackPPARescale2), labels = seq(0.01, 1, by = 0.01))
axis(2, at = 1:nrow(meanPropInNonPackPPARescale2), labels = c(0, 10, 15))

plot(meanPropOutPPARescale2, border = NA, col = viridis(20),
     key = list(side = 3, cex.axis = 0.75), xlab = "Occupancy probability threshold ", 
     ylab = "Buffer size around 'occupied' cells", axis.col = NULL, axis.row = NULL,
     #main = "Proportion of 'non occupied' cells outisde of PPAs")
     main = "")
axis(1, at = 1:ncol(meanPropOutPPARescale2), labels = seq(0.01, 1, by = 0.01))
axis(2, at = 1:nrow(meanPropOutPPARescale2), labels = c(0, 10, 15))

# Combining the results for the three patterns to find the calibration that best represent them all
# Weight to weight presence in PPAs as much as absence outside of PPAs
allPattern <- (meanPropInPackPPARescale2 * 0.5) + (meanPropInNonPackPPARescale2 * 0.5) + meanPropOutPPARescale2
allPattern = (allPattern - min(allPattern)) / (max(allPattern) - min(allPattern))
plot(allPattern, border = NA, col = viridis(20),
     key = list(side = 3, cex.axis = 0.75), xlab = "Occupancy probability threshold ", 
     ylab = "Buffer size around 'occupied' cells", axis.col = NULL, axis.row = NULL,
     #main = "Three paterrns combined")
     main = "")
axis(1, at = 1:ncol(allPattern), labels = seq(0.01, 1, by = 0.01))
axis(2, at = 1:nrow(allPattern), labels = c(0, 10, 15))


##############
# Validation #
##############

# Find the best calibration (where the three patterns combined = 1)
probThresh = 0.51
bufferkm = 10

# Compare the produce maps with these parameter values with the real PPAs
ResPropInPackPPA <- rep(NA, (length(rasterYearCal) - 1))
ResPropInNonPackPPA <- rep(NA, (length(rasterYearCal) - 1))
ResPropOutPPA <- rep(NA, (length(rasterYearCal) - 1))
yearStart <- 2008

for(year in 1:(length(rasterYearCal) - 1)){
  
  # Build a map of PPA with two consecutive winters
  raster2winter <- occ2winters(occWinter1 = rasterYearCal[[year]], 
                               occWinter2 = rasterYearCal[[year + 1]],
                               probThresh = probThresh, bufferkm = bufferkm)
  
  
  # Compare the created PPA with the ones from field data to compute the three patterns
  valid2winter <- validationPPA(raster2winter = raster2winter, PPA = listPPA[[year]], 
                                packPPA = listPackPPA[[year]], nonPackPPA = listNonPackPPA[[year]],
                                sampledCells = cellSampledCal[[year]])
  
  ResPropInPackPPA[year] <- valid2winter[[1]]
  ResPropInNonPackPPA[year] <- valid2winter[[2]]
  ResPropOutPPA[year] <- valid2winter[[3]]
  
  plot(raster2winter, main = paste0("Winters ", yearStart + year - 1, "-", yearStart + year,
                                   " / ", yearStart + year, "-", yearStart + year + 1), legend = FALSE)
  plot(franceShape, add = TRUE)
  plot(listPackPPA[[year]], add = TRUE, lwd = 2, border = "red")
  plot(listNonPackPPA[[year]], add = TRUE, lwd = 2, border = "black")
  legend("topleft", legend = c("'Occupied' cells" ,"Pack PPAs", "Non-pack PPAs"), 
         fill = c("forestgreen","red", "black"))  

}

round((mean(ResPropInPackPPA) * 100)) # % of cells correctly defined as 'occupied' inside pack PPAs
round((mean(ResPropInNonPackPPA) * 100)) # % of cells correctly defined as 'occupied' inside non-pack PPAs
round((mean(ResPropOutPPA) * 100)) # % of cells correctly defined as 'non-occupied' outside of PPAs


# Compute how well the occupancy maps reproduce PPA appareance or disapperance
validDisAppearPPA <- function(raster2winter, PPA, sampledCells){
  
  raster2winter@crs <- CRS("+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
  
  # Crop for the concerned areas
  raster2winter <- crop(raster2winter, sampledCells)
  
  # Proportion of 'occupied' cells in PPAs
  if(!is.null(PPA)){
    PPA <- crop(PPA, sampledCells)
    inPPA <- unlist(extract(raster2winter, PPA))
    propInPPA <- sum(inPPA, na.rm = TRUE) / length(inPPA)
  } else {
    propInPPA <- NA
  }
  return(propInPPA)
}

# PPA when they firstly appear
load("data/newPPA.RData") # listNewPPA, listNewPackPPA, listNewNonPackPPA
# PPA when they firstly disappear
load("data/lostPPA.RData") # listLostPPA, listLostPackPPA, listLostNonPackPPA

resPropInNewPPA <- rep(NA, (length(rasterYearCal) - 1))
resPropInLostPPA <- rep(NA, (length(rasterYearCal) - 1))

for(year in 2:(length(rasterYearCal) - 1)){

  # Build a map of PPA with two consecutive winters
  raster2winter <- occ2winters(occWinter1 = rasterYearCal[[year]], 
                               occWinter2 = rasterYearCal[[year + 1]],
                               probThresh = probThresh, bufferkm = bufferkm)
  
  # How well new PPA appearing are represented by 'occupied' cells
  resPropInNewPPA[year] <- validDisAppearPPA(raster2winter = raster2winter, 
                                             PPA = listNewPPA[[year]], 
                                             sampledCells = cellSampledCal[[year]])
  
  # How well PPA disappearing are represented by 'non-occupied' cells
  resPropInLostPPA[year] <- 1 - validDisAppearPPA(raster2winter = raster2winter, 
                                              PPA = listLostPPA[[year - 1]], 
                                              sampledCells = cellSampledCal[[year]])
}

round((mean(resPropInNewPPA, na.rm = TRUE) * 100)) # % of cells correctly defined as 'occupied' in new PPA (that appear this year)
round((mean(resPropInLostPPA, na.rm = TRUE) * 100)) # % of cells correctly defined as 'non-occupied' in areas where PPA disappear (disappear this year)

# Comparing the growth rates over time from the occupancy maps with those derived from population size computed with CMR models
estimateCMR <- c(17.1, 35.4, 47.7, 25.1, 62.6, 47.9, 81.7, 110.5, 102.7, 135.9, 132.6, 
                 101.7, 130.3, 141.4, 141.5, 175.5, 210.3, 174.5, 353.6, 280.2, 376.7,
                 561.2, 571.9, 682.4, 645.7, 783.8) # Winter 1995/1996 to 2020/2021

cellOcc <- c()
# How many cells are 'occupied' each year
for(year in 1:length(rasterYears[3:28])){ # Winter 1995/1996 to 2020/2021
  
  occWinter <- rasterYears[3:28][[year]]
  occWinter@crs <- CRS("+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
  
  occWinter[which(is.na(raster::values(occWinter)) | raster::values(occWinter) < probThresh)] <- 0
  occWinter[which(raster::values(occWinter) >= probThresh)] <- 1
  occWinter[raster::values(occWinter) == 0] <- NA
  if(bufferkm != 0){
    occWinter <- buffer(occWinter, width = bufferkm * 1000)
  }
  cellOcc <- c(cellOcc, sum(raster::values(occWinter), na.rm = TRUE))
}
growthCellOcc <- cellOcc[2:length(cellOcc)] / cellOcc[1:(length(cellOcc) - 1)]
growthRateCMR <- estimateCMR[2:length(estimateCMR)]/estimateCMR[1:(length(estimateCMR) - 1)]
diffGowthCellOcc <- abs(growthCellOcc - growthRateCMR)
round(mean(diffGowthCellOcc, na.rm = TRUE), digits = 3)
round(sd(diffGowthCellOcc, na.rm = TRUE), digits = 3)

winterYears <- c("96/97", "97/98", "98/99", "99/00", "00/01", "01/02", "02/03", "03/04", "04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21")
plot(x = 1996:2020, y = growthRateCMR, xlab = "Winter", ylab = "Growth rate", type = "l",lwd = 2, xaxt = 'n')
lines(x = 1996:2020, y = growthCellOcc, col = "forestgreen", lwd = 2)
abline(h = 1, lty = 2)
axis(1, at = 1996:2020, labels=winterYears)
legend("topright", legend = c("CMR", "Occupancy"), col = c("black", "forestgreen"), lty = c(1, 1), lwd = c(2, 2))

# Look at the growth compared to t0
occGrowth <-((cellOcc - cellOcc[1]) / cellOcc[1]) * 100
occ0025_ <- occ0025[3:28] # load occ0025 and occ0975 using 4.UseIndicator.R
occ0975_ <- occ0975[3:28]
occGrowth0025 <-((occ0025_ - occ0025_[1]) / occ0025_[1]) * 100
occGrowth0975 <-((occ0975_ - occ0975_[1]) / occ0975_[1]) * 100
CMRGrowth <-((estimateCMR - estimateCMR[1]) / estimateCMR[1]) * 100
cmr25 <- c(1.8, 20.8, 25.9, 10, 32.9, 28.3, 49.8, 75.2, 71.7, 96.9, 94.7, 73.7, 95, 100.1, 105.5,
           132, 161.3, 132.1, 278.6, 217.1, 285.9, 432.2, 460.9, 539.3, 508.6, 641.5)
CMRGrowth0025 <-((cmr25 - cmr25[1]) / cmr25[1]) * 100
cmr975 <- c(33.9, 60.5, 73.2, 59.2, 137.6, 89, 150.1, 217.6, 201.8, 281.9, 258.5, 202.2, 277.6, 288,
            285.5, 358.8, 407.7, 352.5, 706.5, 583.2, 755.4, 883.6, 904.2, 968.6, 934.7, 978.1)
CMRGrowth0975 <-((cmr975 - cmr975[1]) / cmr975[1]) * 100

plot(x = 1995:2020, y = occGrowth, ylim = c(0, max(CMRGrowth)),
     main = "", xlab = "Years", ylab = "Growth compared to t0 in %", type = "l")
# polygon(c(1995:2020, rev(1995:2020)),
#         c(occGrowth0975, rev(occGrowth0025)), col = "gray85", lty = 0)
lines(x = 1995:2020, y = CMRGrowth, lwd = 2)
# polygon(c(1996:2021, rev(1996:2021)), 
#         c(CMRGrowth0975, rev(CMRGrowth0025)), col = "gray85", lty = 0)
