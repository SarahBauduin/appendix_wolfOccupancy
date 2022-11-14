#####################
## Occupancy model ##
#####################

setwd("C:/Users/sarah.bauduin/Documents/GitHub/appendix_wolfOccupancy")

library(nimble)
library(MCMCvis)
library(raster)


# Load data and covariates
load("data/y.RData")
load("data/effort.RData")
load("data/forestCov.RData")
load("data/farmCov.RData")
load("data/pastureCov.RData")
load("data/rockCov.RData")
load("data/altiCov.RData")
load("data/cellRoadCov.RData")
load("data/distBarrierCov.RData")
load("data/yShortDisp.RData")
load("data/yLongDisp.RData")


# Identify cells which were never sampled
effortNA <- effort
effortNA[effortNA == 0] <- NA
noSampl <- apply(effortNA, 1, function(x) all(is.na(x)))
sum(noSampl) # 1110 cells never sampled across all the years

# Remove these cells from y and the covariates
y <- y[-noSampl, , ]
effort <- effort[-noSampl, ]
forestCov <- forestCov[-noSampl, ]
farmCov <- farmCov[-noSampl, ]
pastureCov <- pastureCov[-noSampl, ]
rockCov <- rockCov[-noSampl, ]
altiCov <- altiCov[-noSampl, ]
cellRoadCov <- cellRoadCov[-noSampl, ]
distBarrierCov <- distBarrierCov[-noSampl, ]
yShortDisp <- yShortDisp[-noSampl, ]
yLongDisp <- yLongDisp[-noSampl, ]

# Function to determine whether a cell is sampled or not
is_sampled <- nimbleFunction(
  run = function(effort = double(0)) {
    if(effort == 0) return(0)
    return(1)
    returnType(double())
  })


##########
# Models #
##########

# Full model

model <- nimbleCode({
  
  # Intercept priors
  alpha.psi ~ dunif(0, 1) # Initial occupancy
  alpha.phi ~ dnorm(0, sd = 1.5) # Extinction 
  alpha.gamma ~ dnorm(0, sd = 1.5) # Colonization
  alpha.p ~ dnorm(0, sd = 1.5) # Detection probability
  
  # Covariate coefficient priors
  beta.forestGamma ~ dnorm(0, sd = 1.5) # Forest cover impact on colonization
  beta.farmlandGamma ~ dnorm(0, sd = 1.5) # Farmland cover impact on colonization
  beta.pastureGamma ~ dnorm(0, sd = 1.5) # Pasture cover impact on colonization
  beta.rockGamma ~ dnorm(0, sd = 1.5) # Rock cover impact on colonization
  beta.highElevGamma ~ dnorm(0, sd = 1.5) # Proportion of high elevation impact on colonization
  beta.meanElevGamma ~ dnorm(0, sd = 1.5) # Mean elevation impact on colonization
  beta.distBarrGamma ~ dnorm(0, sd = 1.5) # Distance to barriers impact on colonization
  beta.shortDispGamma ~ dnorm(0, sd = 1.5) # Short dispersal impact on colonization
  beta.longDispGamma ~ dnorm(0, sd = 1.5) # Long dispersal impact on colonization
  
  beta.roads ~ dnorm(0, sd = 1.5) # Road density (= access) impact on detection 
  beta.effort ~ dnorm(0, sd = 1.5) # Sampling effort impact detection
  beta.occ2 ~ dnorm(0, sd = 1.5) # Secondary occasion (months) 
  beta.occ3 ~ dnorm(0, sd = 1.5)
  beta.occ4 ~ dnorm(0, sd = 1.5)
  beta.occ5 ~ dnorm(0, sd = 1.5)
  
  # Ecological submodel: Define state conditional on parameters 
  
  # Initial occupancy constant on all sites (1st year = z[, 1])
  for (i in 1:nSite){
    z[i, 1] ~ dbern(psi1[i]) # State occupied or not
    psi1[i] <- alpha.psi # Occupancy probability
  }
  
  # Then, for the following years, we put the new occupancy estimate
  # from the logistic regression
  for(i in 1:nSite){
    for (k in 2:nYear){
      # Define the new occupancy probability for the year k from the year k-1
      muZ[i, k] <- z[i, k - 1] * (1 - epsilon[k - 1]) + (1 - z[i, k - 1]) * gamma[i, k - 1]
      # New state (occupied or not) for the site i for the year k
      z[i, k] ~ dbern(muZ[i, k])
      # Logistic regression for colonization specific to site i and year k-1
      gamma[i, k - 1] <- 1 / (1 + exp(- lgamma[i, k - 1]))
      lgamma[i, k - 1] <- alpha.gamma + 
        beta.forestGamma * forest[i, k - 1] +
        beta.farmlandGamma * farm[i, k - 1] +
        beta.pastureGamma * pasture[i, k - 1] +
        beta.rockGamma * rock[i, k - 1] +
        beta.highElevGamma * highElev[i] +
        beta.meanElevGamma * meanElev[i] +
        beta.distBarrGamma * distBarr[i, k - 1] +
        beta.shortDispGamma * shortDisp[i, k - 1] +
        beta.longDispGamma * longDisp[i, k - 1]
    } # k years
  } # i sites
  
  # Constant extinction
  for (k in 2:nYear){
    epsilon[k - 1] <- 1 / (1 + exp(- lepsilon[k - 1]))
    lepsilon[k - 1] <- alpha.phi
  }			
  
  # Observation model (detectability)
  for (i in 1:nSite){
    for(j in 1: nRep){
      for (k in 1:nYear){
        # y is the observation
        # If z = 0 then y = 0
        # If z = 1, then muy = p and y ~ dbern(muy * p)
        y[i,j,k] ~ dbern(muy[i,j,k])
        # Multiply the detection probability by the presence/absence on the site
        muy[i, j, k] <- z[i, k] * p[i, j, k]
        # Detection probability
        p[i, j, k] <- (1 / (1 + exp(- lp.lim[i, j, k]))) * is_sampled(effort[i, k])
        lp.lim[i,j,k] <- min(999, max(-999, lp[i, j, k])) # Center pl value between -999 and 999
        lp[i, j, k] <- alpha.p +
          beta.roads * roads[i,k] +
          beta.effort * effort[i,k] +
          beta.occ2 * equals(j, 2) +
          beta.occ3 * equals(j, 3) +
          beta.occ4 * equals(j, 4) +
          beta.occ5 * equals(j, 5) 
      } # k years
    } # j months
  } # i site
})

# Bundle and scale the data
altiCov$meanAlt[is.na(altiCov$meanAlt)] <- 0 # fixing NA elevation
mydata <- list(y = y,
               forest = scale(forestCov), 
               farm = scale(farmCov),
               pasture = scale(pastureCov),
               rock = scale(rockCov),
               highElev = as.numeric(scale(altiCov$propHighAlt)), 
               meanElev = as.numeric(scale(altiCov$meanAlt)),
               distBarr = scale(distBarrierCov),
               shortDisp = yShortDisp, # error when scaling this one so no scale
               longDisp = scale(yLongDisp), 
               roads = scale(cellRoadCov), 
               effort = scale(effort)
)
myconstants <- list(nSite = dim(y)[1], 
                    nRep = dim(y)[2], 
                    nYear = dim(y)[3])

# Initial values
y_aggregated <- apply(y, c(1, 3), sum, na.rm = TRUE)  # Sites by years, number of detections
z_inits <- (y_aggregated > 0) * 1
z_inits[z_inits == 0] <- NA
inits <- function() {list(z = z_inits, 
                          alpha.psi = runif(1, 0, 1), 
                          alpha.phi = rnorm(1, mean = 0, sd = 1.5), 
                          alpha.gamma = rnorm(1, mean = 0, sd = 1.5), 
                          alpha.p = rnorm(1, mean = 0, sd = 1.5))}

# Setup Nimble
occWolf <- nimbleModel(code = model,
                       data = mydata,
                       constants = myconstants,
                       inits = inits(),
                       calculate = FALSE)
CoccWolf <- compileNimble(occWolf)
occWolfConf <- configureMCMC(occWolf)
occWolfConf$addMonitors(c("z"))
occWolfMCMC <- buildMCMC(occWolfConf, useConjugacy = FALSE)
CoccWolfMCMC <- compileNimble(occWolfMCMC,
                              project = occWolf)

# Run Nimble
samples <- runMCMC(mcmc = CoccWolfMCMC, 
                   niter = 10000, 
                   nburnin = 2000,
                   nchains = 2)
save(samples, file = "modelOutputs/fullModel.RData")

# Numerical summary
MCMCsummary(object = samples, round = 2, params = c("alpha.gamma",
                                                    "alpha.p",
                                                    "alpha.phi",
                                                    "alpha.psi",
                                                    "beta.forestGamma",
                                                    "beta.farmlandGamma",
                                                    "beta.pastureGamma",
                                                    "beta.rockGamma",
                                                    "beta.highElevGamma",
                                                    "beta.meanElevGamma",
                                                    "beta.distBarrGamma",
                                                    "beta.shortDispGamma",
                                                    "beta.longDispGamma",
                                                    "beta.roads",
                                                    "beta.effort",
                                                    "beta.occ2",
                                                    "beta.occ3",
                                                    "beta.occ4",
                                                    "beta.occ5"))

# Visualize parameter posterior distributions
MCMCplot(object = samples, params = c("alpha.gamma",
                                      "alpha.p",
                                      "alpha.phi",
                                      "alpha.psi",
                                      "beta.forestGamma",
                                      "beta.farmlandGamma",
                                      "beta.pastureGamma",
                                      "beta.rockGamma",
                                      "beta.highElevGamma",
                                      "beta.meanElevGamma",
                                      "beta.distBarrGamma",
                                      "beta.shortDispGamma",
                                      "beta.longDispGamma",
                                      "beta.roads",
                                      "beta.effort",
                                      "beta.occ2",
                                      "beta.occ3",
                                      "beta.occ4",
                                      "beta.occ5"))

# Check convergence
MCMCtrace(object = samples,
          pdf = FALSE, # no export to PDF
          ind = TRUE, 
          params = c("alpha.gamma",
                     "alpha.p",
                     "alpha.phi",
                     "alpha.psi",
                     "beta.forestGamma",
                     "beta.farmlandGamma",
                     "beta.pastureGamma",
                     "beta.rockGamma",
                     "beta.highElevGamma",
                     "beta.meanElevGamma",
                     "beta.distBarrGamma",
                     "beta.shortDispGamma",
                     "beta.longDispGamma",
                     "beta.roads",
                     "beta.effort",
                     "beta.occ2",
                     "beta.occ3",
                     "beta.occ4",
                     "beta.occ5"))


# Model with only significant variables

model <- nimbleCode({
  
  # Intercept priors
  alpha.psi ~ dunif(0, 1) # Initial occupancy
  alpha.phi ~ dnorm(0, sd = 1.5) # Extinction 
  alpha.gamma ~ dnorm(0, sd = 1.5) # Colonization
  alpha.p ~ dnorm(0, sd = 1.5) # Detection probability
  
  # Covariate coefficient priors
  beta.forestGamma ~ dnorm(0, sd = 1.5) # Forest cover impact on colonization
  beta.pastureGamma ~ dnorm(0, sd = 1.5) # Pasture cover impact on colonization
  beta.highElevGamma ~ dnorm(0, sd = 1.5) # Proportion of high elevation impact on colonization
  beta.meanElevGamma ~ dnorm(0, sd = 1.5) # Mean elevation impact on colonization
  beta.distBarrGamma ~ dnorm(0, sd = 1.5) # Distance to barriers impact on colonization
  beta.shortDispGamma ~ dnorm(0, sd = 1.5) # Short dispersal impact on colonization
  beta.longDispGamma ~ dnorm(0, sd = 1.5) # Long dispersal impact on colonization
  
  beta.roads ~ dnorm(0, sd = 1.5) # Road density (= access) impact on detection 
  beta.effort ~ dnorm(0, sd = 1.5) # Sampling effort impact detection
  beta.occ2 ~ dnorm(0, sd = 1.5) # Secondary occasion (months) 
  beta.occ3 ~ dnorm(0, sd = 1.5)
  beta.occ4 ~ dnorm(0, sd = 1.5)
  beta.occ5 ~ dnorm(0, sd = 1.5)
  
  # Ecological submodel: Define state conditional on parameters 
  
  # Initial occupancy constant on all sites (1st year = z[, 1])
  for (i in 1:nSite){
    z[i, 1] ~ dbern(psi1[i]) # State occupied or not
    psi1[i] <- alpha.psi # Occupancy probability
  }
  
  # Then, for the following years, we put the new occupancy estimate
  # from the logistic regression
  for(i in 1:nSite){
    for (k in 2:nYear){
      # Define the new occupancy probability for the year k from the year k-1
      muZ[i, k] <- z[i, k - 1] * (1 - epsilon[k - 1]) + (1 - z[i, k - 1]) * gamma[i, k - 1]
      # New state (occupied or not) for the site i for the year k
      z[i, k] ~ dbern(muZ[i, k])
      # Logistic regression for colonization specific to site i and year k-1
      gamma[i, k - 1] <- 1 / (1 + exp(- lgamma[i, k - 1]))
      lgamma[i, k - 1] <- alpha.gamma + 
        beta.forestGamma * forest[i, k - 1] +
        beta.pastureGamma * pasture[i, k - 1] +
        beta.highElevGamma * highElev[i] +
        beta.meanElevGamma * meanElev[i] +
        beta.distBarrGamma * distBarr[i, k - 1] +
        beta.shortDispGamma * shortDisp[i, k - 1] +
        beta.longDispGamma * longDisp[i, k - 1]
    } # k years
  } # i sites
  
  # Constant extinction
  for (k in 2:nYear){
    epsilon[k - 1] <- 1 / (1 + exp(- lepsilon[k - 1]))
    lepsilon[k - 1] <- alpha.phi
  }			
  
  # Observation model (detectability)
  for (i in 1:nSite){
    for(j in 1: nRep){
      for (k in 1:nYear){
        # y is the observation
        # If z = 0 then y = 0
        # If z = 1, then muy = p and y ~ dbern(muy * p)
        y[i,j,k] ~ dbern(muy[i,j,k])
        # Multiply the detection probability by the presence/absence on the site
        muy[i, j, k] <- z[i, k] * p[i, j, k]
        # Detection probability
        p[i, j, k] <- (1 / (1 + exp(- lp.lim[i, j, k]))) * is_sampled(effort[i, k])
        lp.lim[i,j,k] <- min(999, max(-999, lp[i, j, k])) # Center pl value between -999 and 999
        lp[i, j, k] <- alpha.p +
          beta.roads * roads[i,k] +
          beta.effort * effort[i,k] +
          beta.occ2 * equals(j, 2) +
          beta.occ3 * equals(j, 3) +
          beta.occ4 * equals(j, 4) +
          beta.occ5 * equals(j, 5) 
      } # k years
    } # j months
  } # i site
})

# Bundle and scale the data
altiCov$meanAlt[is.na(altiCov$meanAlt)] <- 0 # fixing NA elevation
mydata <- list(y = y,
               forest = scale(forestCov), 
               pasture = scale(pastureCov),
               highElev = as.numeric(scale(altiCov$propHighAlt)), 
               meanElev = as.numeric(scale(altiCov$meanAlt)),
               distBarr = scale(distBarrierCov),
               shortDisp = yShortDisp, # error when scaling this one so no scale
               longDisp = scale(yLongDisp), 
               roads = scale(cellRoadCov), 
               effort = scale(effort)
)
myconstants <- list(nSite = dim(y)[1], 
                    nRep = dim(y)[2], 
                    nYear = dim(y)[3])

# Initial values
y_aggregated <- apply(y, c(1, 3), sum, na.rm = TRUE)  # Sites by years, number of detections
z_inits <- (y_aggregated > 0) * 1
z_inits[z_inits == 0] <- NA
inits <- function() {list(z = z_inits, 
                          alpha.psi = runif(1, 0, 1), 
                          alpha.phi = rnorm(1, mean = 0, sd = 1.5), 
                          alpha.gamma = rnorm(1, mean = 0, sd = 1.5), 
                          alpha.p = rnorm(1, mean = 0, sd = 1.5))}

# Setup Nimble
occWolf <- nimbleModel(code = model,
                       data = mydata,
                       constants = myconstants,
                       inits = inits(),
                       calculate = FALSE)
CoccWolf <- compileNimble(occWolf)
occWolfConf <- configureMCMC(occWolf)
occWolfConf$addMonitors(c("z"))
occWolfMCMC <- buildMCMC(occWolfConf, useConjugacy = FALSE)
CoccWolfMCMC <- compileNimble(occWolfMCMC,
                              project = occWolf)

# Run Nimble
samples <- runMCMC(mcmc = CoccWolfMCMC, 
                   niter = 10000, 
                   nburnin = 2000,
                   nchains = 2)
save(samples, file = "modelOutputs/modelSignificantVar.RData")

# Numerical summary
MCMCsummary(object = samples, round = 2, params = c("alpha.gamma",
                                                    "alpha.p",
                                                    "alpha.phi",
                                                    "alpha.psi",
                                                    "beta.forestGamma",
                                                    "beta.pastureGamma",
                                                    "beta.highElevGamma",
                                                    "beta.meanElevGamma",
                                                    "beta.distBarrGamma",
                                                    "beta.shortDispGamma",
                                                    "beta.longDispGamma",
                                                    "beta.roads",
                                                    "beta.effort",
                                                    "beta.occ2",
                                                    "beta.occ3",
                                                    "beta.occ4",
                                                    "beta.occ5"))

# Visualize parameter posterior distributions
MCMCplot(object = samples, params = c("alpha.gamma",
                                      "alpha.p",
                                      "alpha.phi",
                                      "alpha.psi",
                                      "beta.forestGamma",
                                      "beta.pastureGamma",
                                      "beta.highElevGamma",
                                      "beta.meanElevGamma",
                                      "beta.distBarrGamma",
                                      "beta.shortDispGamma",
                                      "beta.longDispGamma",
                                      "beta.roads",
                                      "beta.effort",
                                      "beta.occ2",
                                      "beta.occ3",
                                      "beta.occ4",
                                      "beta.occ5"))

# Check convergence
MCMCtrace(object = samples,
          pdf = FALSE, # no export to PDF
          ind = TRUE, 
          params = c("alpha.gamma",
                     "alpha.p",
                     "alpha.phi",
                     "alpha.psi",
                     "beta.forestGamma",
                     "beta.pastureGamma",
                     "beta.highElevGamma",
                     "beta.meanElevGamma",
                     "beta.distBarrGamma",
                     "beta.shortDispGamma",
                     "beta.longDispGamma",
                     "beta.roads",
                     "beta.effort",
                     "beta.occ2",
                     "beta.occ3",
                     "beta.occ4",
                     "beta.occ5"))


##################
# Occupancy maps #
##################

# Load Nimble results
load("modelOutputs/modelSignificantVar.RData")
# Remove the covariates and keep only the occupancy estimates
chain1 <- samples[[1]][, 18:ncol(samples[[1]])]
chain2 <- samples[[2]][, 18:ncol(samples[[2]])]
dim(chain1) # 8000 161056
# 8000 = 10000 (iterations) - 2000 (burn-in)
# 161056 = 5752 cells * 28 years

# Mean occupancy for each cell for each year from the two chains
# as well as the 2.5% and 97.5% quantiles
# Cannot rbind the two chains because of memory issue so need to use a loop 
# to compute the mean and quantiles
meanLoop <- c()
Q0025Loop <- c()
Q0975Loop <- c()
for(i in 1:ncol(chain1)){
  meanLoop <- c(meanLoop, mean(c(chain1[,i], chain2[,i])))
  Q0025Loop <- c(Q0025Loop, quantile(c(chain1[,i], chain2[,i]), probs = 0.025))
  Q0975Loop <- c(Q0975Loop, quantile(c(chain1[,i], chain2[,i]), probs = 0.975))
}
seqCell <- seq(1, 161056, by = 5752)
listMap <- list()
listMapQ0025 <- list()
listMapQ0975 <- list()
for(i in 1:28){
  listMap[[i]] <- meanLoop[seqCell[i]:(seqCell[i] + 5752 - 1)]
  listMapQ0025[[i]] <- Q0025Loop[seqCell[i]:(seqCell[i] + 5752 - 1)]
  listMapQ0975[[i]] <- Q0975Loop[seqCell[i]:(seqCell[i] + 5752 - 1)]
  print(i)
}

# Create the raster maps of occupancy for each year
load("data/gridFr.RData")
gridRaster <- raster(xmn = extent(gridFr)[1], xmx =  extent(gridFr)[2], 
                     ymn =  extent(gridFr)[3], ymx =  extent(gridFr)[4], resolution = c(10000, 10000))
load("data/franceShape.RData") # shapefile of France
yearStart <- 1993
rasterYears <- list()

# Mean occupancy
for(i in 1:28){
  # Merge the occupancy value with the cell numbers
  df_ <- cbind.data.frame(cellID = rownames(effort), zProb = listMap[[i]], cellSampled = effort[,i])
  df_Full <- merge(df_, data.frame(cellID = 1:length(gridFr)), all.y = TRUE)
  df_FullSort <- df_Full[order(as.numeric(as.character(df_Full$cellID))), ]
  # NA for the cells not sampled during the year
  df_FullSort[is.na(df_FullSort[,"cellSampled"]) | df_FullSort[,"cellSampled"] == 0, "zProb"] <- NA
  
  # Put the occupancy value on the map
  grid_zprob <- rasterize(gridFr, gridRaster, df_FullSort$zProb)
  rasterYears[[i]] <- grid_zprob
  
  plot(rasterYears[[i]], main = paste0("Winter ",yearStart + i - 1, "-", yearStart + i))
  plot(franceShape, add = TRUE)
}
save(rasterYears, file = "modelOutputs/rasterYears.RData")

# 2.5% quantile occupancy
rasterYearsQ0025 <- list()
for(i in 1:28){
  # Merge the occupancy value with the cell numbers
    df_ <- cbind.data.frame(cellID = rownames(effort), zProb = listMapQ0025[[i]], cellSampled = effort[,i])
  df_Full <- merge(df_, data.frame(cellID = 1:length(gridFr)), all.y = TRUE)
  df_FullSort <- df_Full[order(as.numeric(as.character(df_Full$cellID))), ]
  # NA for the cells not sampled during the year
  df_FullSort[is.na(df_FullSort[,"cellSampled"]) | df_FullSort[,"cellSampled"] == 0, "zProb"] <- NA
  
  # Put the occupancy value on the map
  grid_zprob <- rasterize(gridFr, gridRaster, df_FullSort$zProb)
  rasterYearsQ0025[[i]] <- grid_zprob
  
  plot(rasterYearsQ0025[[i]], main = paste0(yearStart + i - 1, "-", yearStart + i))
  plot(franceShape, add = TRUE)
}
save(rasterYearsQ0025, file = "modelOutputs/rasterYearsQ0025.RData")

# 97.5% quantile occupancy
rasterYearsQ0975 <- list()
for(i in 1:28){
  # Merge the occupancy value with the cell numbers
  df_ <- cbind.data.frame(cellID = rownames(effort), zProb = listMapQ0975[[i]], cellSampled = effort[,i])
  df_Full <- merge(df_, data.frame(cellID = 1:length(gridFr)), all.y = TRUE)
  df_FullSort <- df_Full[order(as.numeric(as.character(df_Full$cellID))), ]
  # NA for the cells not sampled during the year
  df_FullSort[is.na(df_FullSort[,"cellSampled"]) | df_FullSort[,"cellSampled"] == 0, "zProb"] <- NA
  
  # Put the occupancy value on the map
  grid_zprob <- rasterize(gridFr, gridRaster, df_FullSort$zProb)
  rasterYearsQ0975[[i]] <- grid_zprob

  plot(rasterYearsQ0975[[i]], main = paste0(yearStart + i - 1, "-", yearStart + i))
  plot(franceShape, add = TRUE)
  
}
save(rasterYearsQ0975, file = "modelOutputs/rasterYearsQ0975.RData")

