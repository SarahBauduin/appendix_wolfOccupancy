
This repository contains all the data and R codes linked to the article **Calibrating an occupancy metric to monitor elusive territorial species at large scale: application to the gray wolf population in France**.

**1.PrepareData.R** is the scripts to produce all the covariates necessary to run the wolf dynamic occupancy model. The variables produced by this scripts are in the **data** folder.

**2.OccupancyModel** uses the files in **data** to run the wolf dynamic occupancy model and produces the maps (mean, 2.5 and 97.5 quantile) of predicted wolf occupancy that are in the **modelOutputs** folder.

**3.MetricCalibration** uses the maps from the **modelOutputs** folder to calibrate the parameters for the different monitoring-focus scenario (producing results stored in the **calibrationOutputs** folder) and choose the best one.

**4.UseIndicator** uses the chosen scenario to produce results and maps for the spatial metric.

 