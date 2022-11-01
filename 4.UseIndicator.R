##########################
## Use of the indicator ##
##########################

setwd("C:/Users/sarah.bauduin/Documents/GitHub/appendix_wolfOccupancy")

library(raster)

# Once parameters are calibrated, we can define cells as 'occupied' and 'non-occupied' each year
# and then compare them from year to year to follow trends in population

# Load the estimated occupancy maps
load("modelOutputs/rasterYears.RData")
load("data/franceShape.RData") # shapefile of France

# The bast calibration from the metric calibration
probThresh = 0.72
bufferkm = 15

# Maps into 'occupied' and 'non-occupied'
rasterYearsTr <- list()

for(i in 1:length(rasterYears)){
  
  occWinter <- rasterYears[[i]]
  occWinter[which(is.na(values(occWinter)) | values(occWinter) < probThresh)] <- 0
  occWinter[which(values(occWinter) >= probThresh)] <- 1
  occWinter[values(occWinter) == 0] <- NA
  if(bufferkm != 0){
    occWinter <- buffer(occWinter, width = bufferkm * 1000)
  }
  occWinter[is.na(occWinter)] <- 0
  
  rasterYearsTr[[i]] <- occWinter
}

# Compare each resulting map with the previous one
yearStart <- 1993

for(i in 1:(length(rasterYearsTr) - 1)){
  
  breakpoints <- c(-1.1, -0.1, 0.1, 1.1) # For raster colors
  colors <- c("red","white","forestgreen")
  
  rasterYearsTrDiff <- rasterYearsTr[[i]]
  rasterYearsTrDiff[] <-  values(rasterYearsTr[[i + 1]]) -  values(rasterYearsTr[[i]])
  plot(rasterYearsTrDiff, main = paste0("Difference between winter ", yearStart + i - 1, "-", yearStart + i,
                                        " \nand winter ", yearStart + i, "-", yearStart + i + 1),
       breaks = breakpoints, col = colors, legend = FALSE)
  plot(franceShape, add = TRUE)
  legend("topleft", legend = c("Newly 'occupied' cells" ,"Newly 'non-occupied' cells"), 
         fill = c("forestgreen","red"))
  
  nCellDiff <- sum(values(rasterYearsTr[[i + 1]]), na.rm = TRUE) - sum(values(rasterYearsTr[[i]]), na.rm = TRUE)
  nCellChange <- sum(values(rasterYearsTrDiff) != 0)
  nCellBecome1 <- sum(values(rasterYearsTrDiff) == -1)
  nCellBecome0 <- sum(values(rasterYearsTrDiff) == 1)
  
  print(paste0("Between winter ", yearStart + i - 1, "-", yearStart + i, 
               " and winter ", yearStart + i, "-", yearStart + i + 1, 
               ", there was ", ifelse(nCellDiff >= 0, "an increase", "a decrease"), 
               " of ", abs(nCellDiff), " 10 x 10 km 'occupied' cells in total over France. Locally, ", 
               nCellChange, " cells changed status (i.e., ", nCellBecome0, " cells were 'occupied' and became 'non-occupied' and ", 
               nCellBecome1, " cells were 'non-occupied' and became 'occupied')."))
  
}


# Changes in occupancy over time
# To include a confidence interval, use the occupancy maps defined with the quantiles
load("modelOutputs/rasterYearsQ0025.RData")
load("modelOutputs/rasterYearsQ0975.RData")

occ0025 <- numeric()
occ0975 <- numeric()
occMean <- numeric()

for(i in 1:length(rasterYearsQ0025)){
  
  # 2.5% confidence interval
  occWinter <- rasterYearsQ0025[[i]]
  occWinter@crs <- CRS("+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
  
  occWinter[which(is.na(values(occWinter)) | values(occWinter) < probThresh)] <- 0
  occWinter[which(values(occWinter) >= probThresh)] <- 1
  occWinter[values(occWinter) == 0] <- NA
  if(bufferkm != 0){
    occWinter <- buffer(occWinter, width = bufferkm * 1000)
  }
  
  occ0025 <- c(occ0025, sum(values(occWinter), na.rm = TRUE))
  
  # 97.5% confidence interval
  occWinter <- rasterYearsQ0975[[i]]
  occWinter@crs <- CRS("+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
  
  occWinter[which(is.na(values(occWinter)) | values(occWinter) < probThresh)] <- 0
  occWinter[which(values(occWinter) >= probThresh)] <- 1
  occWinter[values(occWinter) == 0] <- NA
  if(bufferkm != 0){
    occWinter <- buffer(occWinter, width = bufferkm * 1000)
  }
  
  occ0975 <- c(occ0975, sum(values(occWinter), na.rm = TRUE))
  
  # Mean
  occWinter <- rasterYears[[i]]
  occWinter@crs <- CRS("+proj=lcc +lat_0=46.5 +lon_0=3 +lat_1=49 +lat_2=44 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
  
  occWinter[which(is.na(values(occWinter)) | values(occWinter) < probThresh)] <- 0
  occWinter[which(values(occWinter) >= probThresh)] <- 1
  occWinter[values(occWinter) == 0] <- NA
  if(bufferkm != 0){
    occWinter <- buffer(occWinter, width = bufferkm * 1000)
  }
  
  occMean <- c(occMean, sum(values(occWinter), na.rm = TRUE))
}

plot(x = yearStart:(yearStart + length(occMean) - 1), y = occMean, ylim = c(0, max(occ0975)),
     main = "", xlab = "Years", ylab = "Number of occupied cells", type = "l")
polygon(c(yearStart:(yearStart + length(occMean) - 1), rev(yearStart:(yearStart + length(occMean) - 1))), 
        c(occ0975, rev(occ0025)), col = "gray85", lty = 0)
lines(x = yearStart:(yearStart + length(occMean) - 1), y = occMean, lwd = 2)

