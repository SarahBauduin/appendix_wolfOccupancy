
##This repository contains all the data and R codes linked to the article **Calibrating an occupancy metric to monitor elusive territorial species at large scale: application to the gray wolf population in France**.

**1.PrepareData.R** is the script to produce all the covariates necessary to run the wolf dynamic occupancy model. The variables produced are stored in the **data** folder.

**2.OccupancyModel.R** uses the files in the **data** folder to run the wolf dynamic occupancy model and produce the maps (mean, 2.5th and 97.5th quantiles) of predicted wolf occupancy, stored in the **modelOutputs** folder.

**3.MetricCalibration.R** uses the maps from the **modelOutputs** folder to calibrate parameters of the different monitoring-focus scenarios (results stored in the **calibrationOutputs** folder) and choose the best one.

**4.UseIndicator.R** uses the chosen scenario to produce results and maps of the spatial metric.

 