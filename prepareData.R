##################
## Prepare data ##
##################

setwd("C:/Users/sarah.bauduin/Documents/GitHub/appendix_wolfOccupancy")

library(sp)
library(raster)
library(sf)
library(dplyr)
library(rgeos)

# Grid for the analysis
load("data/gridFr.RData") # 5753 cells


#############
# Wolf data #
#############

wolfData <- read.csv("data/wolfData/wolfDataExtract.csv", header = TRUE, sep = ";")
# Select validated data with a date and coordinates
wolves <- wolfData[wolfData$Fiabilite == "R" & wolfData$date != "" & !is.na(wolfData$X) & !is.na(wolfData$Y), ] # 29719 observations
# Coordinate system for the locations is RGF93/lambert 93
# +proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs
wolvesSP <- SpatialPointsDataFrame(coords = wolves[,c("X", "Y")], data = wolves, 
                                   proj4string = CRS(as.character("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")),
                                   match.ID = TRUE, bbox = NULL)
wolvesSP$date <- as.Date(wolvesSP$date, format = "%Y-%m-%d")

# Fill the grid (which cells had at least one observation)
# for each year and each month of November, December, January, February, March
yearStart <- 1993 # start of the wolf monitoring
yearEnd <- 2020 # winter 2020/2021 is the last available with validated data
y <- array(0, c(length(gridFr), 5, yearEnd-yearStart+1))

for(year in yearStart:yearEnd){
  
  # Extract wolf observations for the different time period
  novYear <- wolvesSP[wolvesSP$date >= as.Date(paste0(year, "-11-01"), format = "%Y-%m-%d") & 
                        wolvesSP$date < as.Date(paste0(year, "-12-01"), format = "%Y-%m-%d"),]
  decYear <- wolvesSP[wolvesSP$date >= as.Date(paste0(year, "-12-01"), format = "%Y-%m-%d") & 
                        wolvesSP$date < as.Date(paste0(year + 1, "-01-01"), format = "%Y-%m-%d"),]
  janYearPlus1 <- wolvesSP[wolvesSP$date >= as.Date(paste0(year + 1, "-01-01"), format = "%Y-%m-%d") & 
                             wolvesSP$date < as.Date(paste0(year + 1, "-02-01"), format = "%Y-%m-%d"),] 
  febYearPlus1 <- wolvesSP[wolvesSP$date >= as.Date(paste0(year + 1, "-02-01"), format = "%Y-%m-%d") & 
                             wolvesSP$date < as.Date(paste0(year + 1, "-03-01"), format = "%Y-%m-%d"),]
  marYearPlus1 <- wolvesSP[wolvesSP$date >= as.Date(paste0(year + 1, "-03-01"), format = "%Y-%m-%d") & 
                             wolvesSP$date < as.Date(paste0(year + 1, "-04-01"), format = "%Y-%m-%d"),]
  
  # Find in which cells observations were made
  if(length(novYear) != 0){
    novCells <- over(novYear, gridFr)
    y[novCells$ID, 1, year - yearStart + 1] <- 1
  }
  if(length(decYear) != 0){
    decCells <- over(decYear, gridFr)
    y[decCells$ID, 2, year - yearStart + 1] <- 1
  }
  if(length(janYearPlus1) != 0){
    janCells <- over(janYearPlus1, gridFr)
    y[janCells$ID, 3, year - yearStart + 1] <- 1
  }
  if(length(febYearPlus1) != 0){
    febCells <- over(febYearPlus1, gridFr)
    y[febCells$ID, 4, year - yearStart + 1] <- 1
  }
  if(length(marYearPlus1) != 0){
    marCells <- over(marYearPlus1, gridFr)
    y[marCells$ID, 5, year - yearStart + 1] <- 1
  }
  
  print(year)
}
save(y, file = "data/y.RData")


##############
# Covariates #
##############

# Sampling effort
load("data/effort.RData")
# df ncol = 28 from winter 1993-1994 to 2020-2021 (5 month considered : November, December, January, February, March)
# Each column has 5753 values = number of grid cells


# CLC from 2000, 2006, 2012 and 2018
clc2000 <- shapefile("data/CLC/CLC2000/CLC00_FR_RGF.shp")
clc2006 <- shapefile("data/CLC/CLC2006/CLC06_FR_RGF.shp")
clc2012 <- shapefile("data/CLC/CLC2012/CLC12_FR_RGF.shp")
clc2018 <- shapefile("data/CLC/CLC2018/clc2018.shp")

# Package sf
gridFr_sf <- st_as_sf(gridFr) %>%
  st_transform(crs = st_crs(st_as_sf(clc2000)))

# Proportion of forest cover (mix, coniferous and deciduous forest)
# 2000
forest2000 <- st_as_sf(clc2000) %>%
  filter(CODE_00 %in% c("311", "312", "313"))
cellForest2000 <- st_intersection(gridFr_sf , forest2000) %>%
  mutate(areaInter = st_area(.)) %>%
  group_by(ID) %>%
  summarise(areaforestCell = sum(areaInter)) %>%
  mutate(propForestCell2000 = areaforestCell / 1e+08) %>%
  as_tibble() %>%
  select(ID, propForestCell2000)
# 2006
forest2006 <- st_as_sf(clc2006) %>%
  filter(CODE_06 %in% c("311", "312", "313"))
cellForest2006 <- st_intersection(gridFr_sf , forest2006) %>%
  mutate(areaInter = st_area(.)) %>%
  group_by(ID) %>%
  summarise(areaforestCell = sum(areaInter)) %>%
  mutate(propForestCell2006 = areaforestCell / 1e+08) %>%
  as_tibble() %>%
  select(ID, propForestCell2006)
# 2012
forest2012 <- st_as_sf(clc2012)%>%
  filter(CODE_12 %in% c("311", "312", "313"))
cellForest2012 <- st_intersection(gridFr_sf , forest2012) %>%
  mutate(areaInter = st_area(.)) %>%
  group_by(ID) %>%
  summarise(areaforestCell = sum(areaInter)) %>%
  mutate(propForestCell2012 = areaforestCell / 1e+08) %>%
  as_tibble() %>%
  select(ID, propForestCell2012)
# 2018
forest2018 <- st_as_sf(clc2018)%>%
  filter(code_18 %in% c("311", "312", "313"))
cellForest2018 <- st_intersection(gridFr_sf , forest2018) %>%
  mutate(areaInter = st_area(.)) %>%
  group_by(ID) %>%
  summarise(areaforestCell = sum(areaInter)) %>%
  mutate(propForestCell2018 = areaforestCell / 1e+08) %>%
  as_tibble() %>%
  select(ID, propForestCell2018)

forestGrid <- data.frame(ID = 1:length(gridFr))
forestGrid <- merge(forestGrid, as.data.frame(cellForest2000), all = TRUE)
forestGrid <- merge(forestGrid, as.data.frame(cellForest2006), all = TRUE)
forestGrid <- merge(forestGrid, as.data.frame(cellForest2012), all = TRUE)
forestGrid <- merge(forestGrid, as.data.frame(cellForest2018), all = TRUE)
forestGrid[is.na(forestGrid)] <- 0
# Transform in the same format as the effort by replicating columns for the missing year
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}
forestCov <- cbind(rep.col(forestGrid[,2], 2000 - yearStart), rep.col(forestGrid[,3], 2006 - 2000),
                   rep.col(forestGrid[,4], 2012 - 2006), rep.col(forestGrid[,5], yearEnd + 1 - 2012))

save(forestCov, file = "data/forestCov.RData")


# Proportion of farmland cover (highly anthropogenic)
# 2000
farm2000 <- st_as_sf(clc2000) %>%
  filter(CODE_00 %in% c("211", "212", "213", "221", "222", "223", "241", "242", "244"))
cellFarm2000 <- st_intersection(gridFr_sf , farm2000) %>%
  mutate(areaInter = st_area(.)) %>%
  group_by(ID) %>%
  summarise(areafarmCell = sum(areaInter)) %>%
  mutate(propFarmCell2000 = areafarmCell / 1e+08) %>%
  as_tibble() %>%
  select(ID, propFarmCell2000)
# 2006
farm2006 <- st_as_sf(clc2006) %>%
  filter(CODE_06 %in% c("211", "212", "213", "221", "222", "223", "241", "242", "244"))
cellFarm2006 <- st_intersection(gridFr_sf , farm2006) %>%
  mutate(areaInter = st_area(.)) %>%
  group_by(ID) %>%
  summarise(areafarmCell = sum(areaInter)) %>%
  mutate(propFarmCell2006 = areafarmCell / 1e+08) %>%
  as_tibble() %>%
  select(ID, propFarmCell2006)
# 2012
farm2012 <- st_as_sf(clc2012) %>%
  filter(CODE_12 %in% c("211", "212", "213", "221", "222", "223", "241", "242", "244"))
cellFarm2012 <- st_intersection(gridFr_sf , farm2012) %>%
  mutate(areaInter = st_area(.)) %>%
  group_by(ID) %>%
  summarise(areafarmCell = sum(areaInter)) %>%
  mutate(propFarmCell2012 = areafarmCell / 1e+08) %>%
  as_tibble() %>%
  select(ID, propFarmCell2012)
# 2018
farm2018 <- st_as_sf(clc2018) %>%
  filter(code_18 %in% c("211", "212", "213", "221", "222", "223", "241", "242", "244"))
cellFarm2018 <- st_intersection(gridFr_sf , farm2018) %>%
  mutate(areaInter = st_area(.)) %>%
  group_by(ID) %>%
  summarise(areafarmCell = sum(areaInter)) %>%
  mutate(propFarmCell2018 = areafarmCell / 1e+08) %>%
  as_tibble() %>%
  select(ID, propFarmCell2018)

farmGrid <- data.frame(ID = 1:5753)
farmGrid <- merge(farmGrid, as.data.frame(cellFarm2000), all = TRUE)
farmGrid <- merge(farmGrid, as.data.frame(cellFarm2006), all = TRUE)
farmGrid <- merge(farmGrid, as.data.frame(cellFarm2012), all = TRUE)
farmGrid <- merge(farmGrid, as.data.frame(cellFarm2018), all = TRUE)
farmGrid[is.na(farmGrid)] <- 0
farmCov <- cbind(rep.col(farmGrid[,2], 2000 - yearStart), rep.col(farmGrid[,3], 2006 - 2000),
                 rep.col(farmGrid[,4], 2012 - 2006), rep.col(farmGrid[,5], yearEnd + 1 - 2012))
save(farmCov, file = "data/farmCov.RData")


# Proportion of pasture lands
# 2000
pasture2000 <- st_as_sf(clc2000) %>%
  filter(CODE_00 %in% c("231", "243", "321", "322", "323"))
cellPasture2000 <- st_intersection(gridFr_sf , pasture2000) %>%
  mutate(areaInter = st_area(.)) %>%
  group_by(ID) %>%
  summarise(areaPastureCell = sum(areaInter)) %>%
  mutate(propPastureCell2000 = areaPastureCell / 1e+08) %>%
  as_tibble() %>%
  select(ID, propPastureCell2000)
# 2006
pasture2006 <- st_as_sf(clc2006) %>%
  filter(CODE_06 %in% c("231", "243", "321", "322", "323"))
cellPasture2006 <- st_intersection(gridFr_sf , pasture2006) %>%
  mutate(areaInter = st_area(.)) %>%
  group_by(ID) %>%
  summarise(areaPastureCell = sum(areaInter)) %>%
  mutate(propPastureCell2006 = areaPastureCell / 1e+08) %>%
  as_tibble() %>%
  select(ID, propPastureCell2006)
# 2012
pasture2012 <- st_as_sf(clc2012) %>%
  filter(CODE_12 %in% c("231", "243", "321", "322", "323"))
cellPasture2012 <- st_intersection(gridFr_sf , pasture2012) %>%
  mutate(areaInter = st_area(.)) %>%
  group_by(ID) %>%
  summarise(areaPastureCell = sum(areaInter)) %>%
  mutate(propPastureCell2012 = areaPastureCell / 1e+08) %>%
  as_tibble() %>%
  select(ID, propPastureCell2012)
# 2018
pasture2018 <- st_as_sf(clc2018) %>%
  filter(code_18 %in% c("231", "243", "321", "322", "323"))
cellPasture2018 <- st_intersection(gridFr_sf , pasture2018) %>%
  mutate(areaInter = st_area(.)) %>%
  group_by(ID) %>%
  summarise(areaPastureCell = sum(areaInter)) %>%
  mutate(propPastureCell2018 = areaPastureCell / 1e+08) %>%
  as_tibble() %>%
  select(ID, propPastureCell2018)

pastureGrid <- data.frame(ID = 1:5753)
pastureGrid <- merge(pastureGrid, as.data.frame(cellPasture2000), all = TRUE)
pastureGrid <- merge(pastureGrid, as.data.frame(cellPasture2006), all = TRUE)
pastureGrid <- merge(pastureGrid, as.data.frame(cellPasture2012), all = TRUE)
pastureGrid <- merge(pastureGrid, as.data.frame(cellPasture2018), all = TRUE)
pastureGrid[is.na(pastureGrid)] <- 0
pastureCov <- cbind(rep.col(pastureGrid[,2], 2000 - yearStart), rep.col(pastureGrid[,3], 2006 - 2000),
                    rep.col(pastureGrid[,4], 2012 - 2006), rep.col(pastureGrid[,5], yearEnd + 1 - 2012))
save(pastureCov, file = "data/pastureCov.RData")


# Proportion or rock cover
# 2000
rock2000 <- st_as_sf(clc2000) %>%
  filter(CODE_00 == "332")
cellRock2000 <- st_intersection(gridFr_sf , rock2000) %>%
  mutate(areaInter = st_area(.)) %>%
  group_by(ID) %>%
  summarise(arearockCell = sum(areaInter)) %>%
  mutate(propRockCell2000 = arearockCell / 1e+08) %>%
  as_tibble() %>%
  select(ID, propRockCell2000)
# 2006
rock2006 <- st_as_sf(clc2006) %>%
  filter(CODE_06 == "332")
cellRock2006 <- st_intersection(gridFr_sf , rock2006) %>%
  mutate(areaInter = st_area(.)) %>%
  group_by(ID) %>%
  summarise(arearockCell = sum(areaInter)) %>%
  mutate(propRockCell2006 = arearockCell / 1e+08) %>%
  as_tibble() %>%
  select(ID, propRockCell2006)
# 2012
rock2012 <- st_as_sf(clc2012) %>%
  filter(CODE_12 == "332")
cellRock2012 <- st_intersection(gridFr_sf , rock2012) %>%
  mutate(areaInter = st_area(.)) %>%
  group_by(ID) %>%
  summarise(arearockCell = sum(areaInter)) %>%
  mutate(propRockCell2012 = arearockCell / 1e+08) %>%
  as_tibble() %>%
  select(ID, propRockCell2012)
# 2018
rock2018 <- st_as_sf(clc2018) %>%
  filter(code_18 == "332")
cellRock2018 <- st_intersection(gridFr_sf , rock2018) %>%
  mutate(areaInter = st_area(.)) %>%
  group_by(ID) %>%
  summarise(arearockCell = sum(areaInter)) %>%
  mutate(propRockCell2018 = arearockCell / 1e+08) %>%
  as_tibble() %>%
  select(ID, propRockCell2018)

rockGrid <- data.frame(ID = 1:5753)
rockGrid <- merge(rockGrid, as.data.frame(cellRock2000), all = TRUE)
rockGrid <- merge(rockGrid, as.data.frame(cellRock2006), all = TRUE)
rockGrid <- merge(rockGrid, as.data.frame(cellRock2012), all = TRUE)
rockGrid <- merge(rockGrid, as.data.frame(cellRock2018), all = TRUE)
rockGrid[is.na(rockGrid)] <- 0
rockCov <- cbind(rep.col(rockGrid[,2], 2000 - yearStart), rep.col(rockGrid[,3], 2006 - 2000),
                 rep.col(rockGrid[,4], 2012 - 2006), rep.col(rockGrid[,5], yearEnd + 1 - 2012))
save(rockCov, file = "data/rockCov.RData")


# Elevation data
elev1 <- raster("data/elevation/032ab314564b9cb72c98fbeb093aeaf69720fbfd/eu_dem_v11_E30N20.TIF")
elev2 <- raster("data/elevation/97824c12f357f50638d665b5a58707cd82857d57/eu_dem_v11_E40N20.TIF")
elev3 <- raster("data/elevation/fdac92901442ed49bf899a997791c5a7faac0b87/eu_dem_v11_E40N30.TIF")
elev4 <- raster("data/elevation/6eb56e1310dfacaef22943d329e3e23f4665de04/eu_dem_v11_E30N30.TIF")
elevation <- merge(elev1, elev2, elev3, elev4)

# Proportion of high altitude (altitude higher than 2500 m)
# and mean altitude
gridFr_sf <- gridFr_sf %>%
  st_transform(crs = elevation@crs)
# Cannot do all cells at once (memory issue)
propHighAlt <- numeric()
meanAlt <- numeric()
# Need to split grid cells into chunks of cells
cutSeq <- seq(1, 5753, 100)
cutCells <- list()
for(i in 1:(length(cutSeq) - 1)){
  cutCells[[i]] <- seq(cutSeq[i], cutSeq[i + 1] - 1, by = 1)
}
cutCells[[58]] <- 5701:5753
for(j in 1:length(cutCells)){
  rasCell <- extract(elevation, gridFr_sf[cutCells[[j]], ])
  nCells <- unlist(lapply(rasCell, FUN = function(x) length(x)))
  nHighAlt <- unlist(lapply(rasCell, FUN = function(x) length(x[!is.na(x) & x > 2500])))
  propHighAlt_ <- nHighAlt / nCells
  meanAlt_ <- unlist(lapply(rasCell, FUN = function(x) mean(x, na.rm = TRUE)))
  propHighAlt <- c(propHighAlt, propHighAlt_)
  meanAlt <- c(meanAlt, meanAlt_)
  
  save(propHighAlt, meanAlt, file = "data/altiCov.RData")
  print(j)
}

altiCov <- cbind.data.frame(propHighAlt = propHighAlt, meanAlt = meanAlt)
save(altiCov, file = "data/altiCov.RData")


# Road data
allRoads2012 <- shapefile("data/roads/2012/TRONCON_ROUTE.shp")
allRoads2015 <- shapefile("data/roads/2015/TRONCON_ROUTE.shp")
allRoads2018 <- shapefile("data/roads/2018/TRONCON_ROUTE.shp")

# Road density (without highways)
roads2012 <- allRoads2012[!allRoads2012$VOCATION == "Type autoroutier", ]
roads2015 <- allRoads2015[!allRoads2015$VOCATION == "Type autoroutier", ]
roads2018 <- allRoads2018[!allRoads2018$VOCATION == "Type autoroutier", ]

gridFr_sf <- gridFr_sf %>%
  st_transform(crs = roads2012@proj4string)
# 2012
cellRoad2012 <- st_intersection(gridFr_sf, st_as_sf(roads2012)) %>%
  mutate(length = st_length(.)) %>%
  group_by(ID) %>%
  summarise(lengthRoadCell = sum(length, na.rm = TRUE)) %>%
  mutate(densRoadCell2012 = (lengthRoadCell / 1000) / 100) %>%
  as_tibble() %>%
  select(ID, densRoadCell2012)
# 2015
cellRoad2015 <- st_intersection(gridFr_sf, st_as_sf(roads2015)) %>%
  mutate(length = st_length(.)) %>%
  group_by(ID) %>%
  summarise(lengthRoadCell = sum(length, na.rm = TRUE)) %>%
  mutate(densRoadCell2015 = (lengthRoadCell / 1000) / 100) %>%
  as_tibble() %>%
  select(ID, densRoadCell2015)
# 2018
cellRoad2018 <- st_intersection(gridFr_sf, st_as_sf(roads2018)) %>%
  mutate(length = st_length(.)) %>%
  group_by(ID) %>%
  summarise(lengthRoadCell = sum(length, na.rm = TRUE)) %>%
  mutate(densRoadCell2018 = (lengthRoadCell / 1000) / 100) %>%
  as_tibble() %>%
  select(ID, densRoadCell2018)

cellRoadGrid <- data.frame(ID = 1:5753)
cellRoadGrid <- merge(cellRoadGrid, as.data.frame(cellRoad2012), all = TRUE)
cellRoadGrid <- merge(cellRoadGrid, as.data.frame(cellRoad2015), all = TRUE)
cellRoadGrid <- merge(cellRoadGrid, as.data.frame(cellRoad2018), all = TRUE)
cellRoadGrid[is.na(cellRoadGrid)] <- 0
cellRoadCov <- cbind(rep.col(cellRoadGrid[,2], 2012 - yearStart), rep.col(cellRoadGrid[,3], 2015-2012),
                     rep.col(cellRoadGrid[,4], yearEnd + 1 - 2015))
save(cellRoadCov, file = "data/cellRoadCov.RData")


# Distance to barriers
# Waterways
waterways <- shapefile("data/waterways/COURS_D_EAU.SHP")
fleuves <- waterways[waterways$TOPONYME %in% c("La Seine", "La Garonne", "La Loire", "Le Rh\xf4ne", "Le Rhin"), ]
fleuves$length <- gLength(fleuves, byid = TRUE)
fleuves5 <- fleuves[fleuves$length > 15000, ]
# Highways
highways2012 <- allRoads2012[allRoads2012$VOCATION == "Type autoroutier", ]
highways2015 <- allRoads2015[allRoads2015$VOCATION == "Type autoroutier", ]
highways2018 <- allRoads2018[allRoads2018$VOCATION == "Type autoroutier", ]

# Distance to the closest barrier (highways or waterways)
centroidGrid_sf <- st_as_sf(gridFr) %>%
  st_transform(crs = highways2012@proj4string) %>%
  st_centroid()
# 2012
linear2012 <- st_union(st_as_sf(highways2012), st_as_sf(fleuves5))
minDist2012 <- numeric()
for(j in 1:length(cutCells)){
  dist2012 <- centroidGrid_sf[cutCells[[j]], ] %>%
    st_distance(linear2012)
  
  dist2012Mat <- matrix(as.numeric(dist2012),
                        nrow = nrow(dist2012),
                        ncol = ncol(dist2012))
  
  minDist2012Cells <- apply(dist2012Mat, 1, min) / 1000
  minDist2012 <- c(minDist2012, minDist2012Cells)
  
  save(minDist2012, file = "data/barrierDistCov.RData")
  print(j)
}

# 2015
linear2015 <- st_union(st_as_sf(highways2015), st_as_sf(fleuves5))
minDist2015 <- numeric()
for(j in 1:length(cutCells)){
  dist2015 <- centroidGrid_sf[cutCells[[j]], ] %>%
    st_distance(linear2015)
  
  dist2015Mat <- matrix(as.numeric(dist2015),
                        nrow = nrow(dist2015),
                        ncol = ncol(dist2015))
  
  minDist2015Cells <- apply(dist2015Mat, 1, min) / 1000
  minDist2015 <- c(minDist2015, minDist2015Cells)
  
  save(minDist2012, minDist2015, file = "data/barrierDistCov.RData")
  print(j)
}

# 2018
linear2018 <- st_union(st_as_sf(highways2018), st_as_sf(fleuves5))
minDist2018 <- numeric()
for(j in 1:length(cutCells)){
  dist2018 <- centroidGrid_sf[cutCells[[j]], ] %>%
    st_distance(linear2018)
  
  dist2018Mat <- matrix(as.numeric(dist2018),
                        nrow = nrow(dist2018),
                        ncol = ncol(dist2018))
  
  minDist2018Cells <- apply(dist2018Mat, 1, min) / 1000
  minDist2018 <- c(minDist2018, minDist2018Cells)
  
  save(minDist2012, minDist2015, minDist2018, file = "data/barrierDistCov.RData")
  print(j)
}

distBarrierGrid <- cbind.data.frame(minDist2012, minDist2015, minDist2018)
distBarrierCov <- cbind(rep.col(distBarrierGrid[,1], 2012 - yearStart), rep.col(distBarrierGrid[,2], 2015-2012),
                        rep.col(distBarrierGrid[,3], yearEnd + 1 - 2015))
save(distBarrierCov, file = "data/distBarrierCov.RData")


# Short dispersal
# Number of neighbouring cells with observations
yCombined <- matrix(nrow = dim(y)[1], ncol = dim(y)[3])
for(k in 1:dim(y)[3]){
  rowSumTemp <- rowSums(y[, , k], na.rm = TRUE)
  rowSumTemp[rowSumTemp >= 1] <- 1
  yCombined[, k] <- rowSumTemp
}

yShortDisp <- yCombined
yShortDisp[,] <- 0

for(j in 1:dim(yCombined)[2]){
  
  gridFr_sfj <- gridFr_sf %>%
    mutate(obsOcc = yCombined[,j])
  yjBuffer <- gridFr_sf %>%
    mutate(obsOcc = yCombined[,j]) %>%
    st_centroid() %>%
    st_buffer(dist = 10000)
  yj <- st_intersection(yjBuffer, gridFr_sfj) %>%
    filter(ID != ID.1) %>%
    group_by(ID) %>%
    summarise(nCellsAround = sum(obsOcc.1, na.rm = TRUE)) %>%
    as_tibble() %>%
    select(ID, nCellsAround)
  allCellsAround <- left_join(x = tibble(ID = 1:dim(yCombined)[1]), y = yj[,c("ID", "nCellsAround")])
  allCellsAround[is.na(allCellsAround$nCellsAround), "nCellsAround"] <- 0
  yShortDisp[, j] <- allCellsAround$nCellsAround
  print(j)
}
save(yShortDisp, file = "data/yShortDisp.RData")


# Long dispersal
# Number of cells within a 150 km radius
yLongDisp <- yCombined
yLongDisp[,] <- 0

for(j in 1:dim(yCombined)[2]){
  
  gridFr_sfj <- gridFr_sf %>%
    mutate(obsOcc = yCombined[,j])
  yjBuffer <- gridFr_sf %>%
    mutate(obsOcc = yCombined[,j]) %>%
    st_centroid() %>%
    st_buffer(dist = 150000)
  yj <- st_intersection(yjBuffer, gridFr_sfj) %>%
    filter(ID != ID.1) %>%
    group_by(ID) %>%
    summarise(nCellsAround = sum(obsOcc.1, na.rm = TRUE)) %>%
    as_tibble() %>%
    select(ID, nCellsAround)
  allCellsAround <- left_join(x = tibble(ID = 1:dim(yCombined)[1]), y = yj[,c("ID", "nCellsAround")])
  allCellsAround[is.na(allCellsAround$nCellsAround), "nCellsAround"] <- 0
  yLongDisp[, j] <- allCellsAround$nCellsAround
  print(j)
}
save(yLongDisp, file = "data/yLongDisp.RData")

