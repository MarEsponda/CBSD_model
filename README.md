CBSD Model for the SEA

In this project, the area affected by CBSD was estimated for four countries: Thailand, Cambodia, Laos, and Vietnam. The analysis was based on occurrence data of Cassava Mosaic Disease (CMD) recorded in these four countries from 2019 to 2022. Using the model, the affected areas were projected for a 10-year period following the arrival of the disease, with the initial infection point located in southern Cambodia to simulate the introduction of CMD. This process was performed in the RStudio environment, version 4.3.1, utilizing commonly used libraries within the R environment.

The first script to execute is Step_1.R, which creates the necessary matrices for running the model. The input data corresponds to the raster file "Study_zone1.tif", included in the repository. All files in the repository must be downloaded into a single folder, which will serve as the working directory. The path to this directory should be defined by replacing line 15 (setwd("your_file_location")). The script generates multiple files of type bigmatrix, and the matrix creation process can be time-consuming; thus, it is recommended to use a machine with robust computational resources.
Once the first script has been executed, it is possible to run Step_2.R, which loads the matrices created in the first script. This script evaluates the propagation trends, calculates the RMSE for each iteration, and finally exports a map for each cultivation cycle.

Contact:
Maria del Mar Esponda
m.esponda@cgiar.org



