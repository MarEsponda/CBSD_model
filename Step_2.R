#Code 2
#Cassava disease dispersion model
#Inverse Power Version
#Maria del Mar Esponda 2024

library(sf)           # Spatial data manipulation
library(terra)        # Handling and analysis of raster data
library(raster)       # Library for raster data
library(geosphere)    # Geographic calculations
library(dplyr)        # Data manipulation
library(units)        # Handling units
library(bigmemory)    # Managing large datasets in memory
library(bench)        # Quantifies process runtime
library(Metrics)      # Calculates Error

# Load matrices----------------------------------------------------------------------------------------------------------------
setwd("your_file_location")                                                     # Set directory
distpi2.5 <- attach.big.matrix("distpi2.5.desc")                                # Load input matrix from Blocks
distkm <- attach.big.matrix("distanciaskm.desc")
within_distance <- attach.big.matrix("within_distance.desc")
prob_infection <- attach.big.matrix("prob_infection.desc")

ifelse(distpi2.5[2,1]==0.0853 & distpi2.5[1,1]==1,"Correct","Incorrect matrix") # Control point to ensure correct input usage

# Set a threshold for distances greater than 350 km (based on literature, Delaquis 2018) -> Correct
#"No trade events were observed at > 500 km distance, so probabilities of exchange between nodes > 500 km apart were set to zero" Andersen, 2024 -> Incorrect
pix <- 140810
for(i in 1:pix){
  condition <- distkm[i,] > 350                                                # All distances greater than 350 km
  indices <- which(condition, arr.ind = TRUE)                                  # List positions where the condition is met
  distpi3.5[i,indices] <- 0                                                    # Assign 0 to distances greater than 350 km
  print(paste0("cycle", i))
}

# Infection trendline---------------------------------------------------------
# The exponential trendline is not possible in Excel; coefficients are calculated with nls 
x <- c(2019, 2020, 2021, 2022)
x_normalized <- x - min(x)
y <- c(238372.3, 2040803.1, 2951148.5, 4025957.1)                              # Affected area according PestDisplace observations             

# Exponential fit (y = a * exp(b * x))
y_normalized <- y / max(y)  # Normalize between 0 and 1
model <- nls(y_normalized ~ a * exp(b * x_normalized), start = list(a = 0.1, b = 0.001))

# Restore original scale
coefficients <- coef(model)
coefficients["a"] <- coefficients["a"] * max(y)

# Model summary
summary(model)
rmse(model$pixels, model$acerc)

# Infected area per year-----------------------------------------------------------------
# Calculate the assumed infected area based on PestDisplace data
# This will serve as the model input to project cumulative infection over 10 years
ADM3 <- st_read("ADM_3union2.shp")
ADM3 <- st_zm(ADM3)  # Remove Z/M coordinates
geodatabase_path <- "your_file_location/LaoVillagesShp.gdb"
points <- st_read(dsn = geodatabase_path, layer = "Points")
p19 <- st_read("Priv_2019.shp")
p20 <- st_read("Priv_2020.shp")
p21 <- st_read("Priv_2021.shp")
p22 <- st_read("Priv_2022.shp")

points_list <- list("2019" = p19,
                    "2020" = p20,
                    "2021" = p21,
                    "2022" = p22)
points_list <- lapply(points_list, st_zm)

results <- data.frame(Year = character(), Total_Area = numeric())

for (year in names(points_list)) {
  points <- points_list[[year]]
  # Find indices of polygons containing points
  indices <- st_intersects(ADM3, points)
  # Identify polygons with at least one point
  polygons_with_points <- ADM3[lengths(indices) > 0, ]
  total_area <- sum(polygons_with_points$AREA_ha, na.rm = TRUE)  
  results <- rbind(results, data.frame(Year = year, Total_Area = total_area))
}

# Show results
print(results)

# Disease spread simulation----------------------------------
# Start infection in some pixels (pixels 1 and 2)          
set.seed(123)
pix <- 140810
initial_infected <- sample(1:pix, 2, replace = FALSE)                             # IDs of pixels that will start spreading the disease
infection_status <- rep(0, pix)                                                   # Vector of size equal to the number of pixels (initial disease status)
infection_status[initial_infected] <- 1                                           # Mark initial pixels as infected (equal to 1)

num_steps <- 10                                                                   # Number of iterations (each corresponds to a 10-month cycle)
monthly_steps <- 10                                                               # Number of sub-iterations per cycle
whitefly <- 3                                                                     # Distance threshold in km for new infections
prob_threshold <- 0.5                                                             # Probability threshold for infection randomization
results <- data.frame(Cycle = rep(NA, 10), Infected_Pixels = rep(NA, 10),
                      cycle1_fly = rep(NA, 10), cycle2_fly = rep(NA, 10),
                      cycle3_fly = rep(NA, 10), cycle4_fly = rep(NA, 10),
                      cycle5_fly = rep(NA, 10), cycle6_fly = rep(NA, 10),
                      cycle7_fly = rep(NA, 10), cycle8_fly = rep(NA, 10),
                      cycle9_fly = rep(NA, 10), cycle10_fly = rep(NA, 10))       # Number of infected pixels per major cycle
print(table(infection_status))                                                    # Check the number of infected pixels before starting iterations (1 = infected, 0 = healthy)

# List to store infection status for each monthly step
infection_status_list <- vector("list", num_steps)

for (step in 1:num_steps) {
  new_infections <- rep(0, pix)                                                   # Vector of NEW infections based on probability
  for (i in which(infection_status == 1)) {                                       # Loop through currently infected pixels (seed infection)
    prob_infection <- distpi2.5[i, ]                                              # Row vector with weights for the infected district i
    new_infections <- new_infections | (runif(pix) < prob_infection)              # New infections based on probability
  }
  infection_status <- infection_status | new_infections                           # Update infection status with new infections based on probability
  print(table(infection_status))                                                  # Check the number of infected pixels (1 = infected, 0 = healthy)
  
  monthly_infection_status <- vector("list", monthly_steps)                       # Create a list with 10 objects to store monthly data
  for (month in 1:monthly_steps) {                                                # Whitefly infection
    new_infections_dist <- rep(0, pix)                                            # Vector of NEW infections based on distance
    for (i in which(infection_status == 1)) {                                     # Loop through currently infected pixels
      within_distance <- distkm[i, ] < whitefly                                   # Logical vector for pixels within the distance threshold
      randomized <- runif(sum(within_distance)) < prob_threshold                  # Randomize pixels within the distance threshold
      new_infections_dist[within_distance] <- new_infections_dist[within_distance] | randomized
    }
    infection_status <- infection_status | new_infections_dist                    # Update infection status with new distance-based infections
    print(paste0("Whitefly spread, cycle ", month))
    print(table(infection_status))                                                # Check the number of infected pixels
    results[step, (month+2)] <- sum(infection_status)
    monthly_infection_status[[month]] <- infection_status                         # Store logical vector of infections month by month
  }
  infection_status_list[[step]] <- monthly_infection_status                       # Store monthly reports for each cycle
  results$Cycle[step] <- step
  results$Infected_Pixels[step] <- sum(infection_status)
  print(paste("Cycle:", step))
  print(table(infection_status))
}

# The rest of the script follows the same logic.

# Grid Search Optimization--------------------------
# Search for the best prob_threshold value to fit the model results
prob_values <- seq(0.75, 1, by = 0.05)                                             # Threshold probabilities from 0.01 to 0.2, in steps of 0.01 (later refined to 0.08, 0.12, by = 0.001)
grid_results <- data.frame(prob_threshold = prob_values, rmse = NA)               # Create a data frame to store the results
simulate_and_evaluate <- function(prob_threshold) {                               # Function to simulate disease spread and calculate RMSE
  set.seed(123)                                                                   # For reproducibility
  pix <- 140810                                                                   # Number of pixels
  initial_infected <- sample(1:pix, 2, replace = FALSE)                           # IDs of pixels that will start spreading the disease
  infection_status <- rep(0, pix)                                                 # Vector with 0s representing the initial disease status
  infection_status[initial_infected] <- 1                                         # Set initial pixels as infected (equal to 1)
  num_steps <- 10
  monthly_steps <- 10
  whitefly_distance <- 3
  results <- data.frame(Cycle = rep(NA, 10), Infected_Pixels = rep(NA, 10),
                        cycle1_fly = rep(NA, 10), cycle2_fly = rep(NA, 10),
                        cycle3_fly = rep(NA, 10), cycle4_fly = rep(NA, 10),
                        cycle5_fly = rep(NA, 10), cycle6_fly = rep(NA, 10),
                        cycle7_fly = rep(NA, 10), cycle8_fly = rep(NA, 10),
                        cycle9_fly = rep(NA, 10), cycle10_fly = rep(NA, 10)) 
  for (step in 1:num_steps) {
    new_infections <- rep(0, pix)                                                 # Empty vector for NEW infections based on probability
    for (i in which(infection_status == 1)) {                                     # Loop through each currently infected pixel (seed infection)
      prob_infection <- distpi2.5[i, ]                                            # Row vector of weights for the infected district i
      new_infections <- new_infections | (runif(pix) < prob_infection)            # Generate new infections based on probability (every 10 months)
    }
    infection_status <- infection_status | new_infections                         # Update infection status with new probability-based infections
    print(table(infection_status))                                                # Check the number of infected pixels (1 = infected, 0 = healthy)
    monthly_infection_status <- vector("list", monthly_steps)                     # Create a list with 10 elements, each storing the vector of infected pixels for that month
    for (month in 1:monthly_steps) {                                              # Whitefly infection
      new_infections_dist <- rep(0, pix)                                          # Empty vector for NEW infections based on distance
      for (i in which(infection_status == 1)) {                                   # Loop through each currently infected pixel
        within_distance <- distkm[i, ] < whitefly_distance                        # Pixels within the distance threshold
        randomized <- runif(sum(within_distance)) < prob_threshold                # Apply randomness to pixels within the distance
        new_infections_dist[within_distance] <- new_infections_dist[within_distance] | randomized
      }
      infection_status <- infection_status | new_infections_dist                  # Update infection status with new distance-based infections
      print(paste0("Whitefly spread, cycle ", month))
      print(table(infection_status))                                              # Check the number of infected pixels
    }
    results$Cycle[step] <- step
    results$Infected_Pixels[step] <- sum(infection_status)
    print(paste("Cycle:", step))
    print(table(infection_status))
  }
  # Complete simulation, now calculate RMSE with actual data
  model <- data.frame(month = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100), 
                      pixels = c(1109, 1793, 2901, 4692, 7590, 12277, 19859, 32122, 51958, 84043))
  model$approximation <- results$Infected_Pixels
  return(rmse(model$pixels, model$approximation))
}

# Iterate over all prob_threshold values
for (i in seq_along(prob_values)) {
  grid_results$rmse[i] <- simulate_and_evaluate(prob_values[i])                   # Populate the table with RMSE results for each iteration
}

# Obtained values
grid_results <- data.frame(prob_threshold = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1),
                           rmse = c(33635.91, 29044.17, 22708.11, 20187.84, 17554.18, 15327.67, 14694.19, 13272.29, 12947.45, 13391.78, 12838.93, 13467.55, 13391.09, 13600.88, 14323.64, 14300.39, 14689.67, 15002.07, 14737.2, 15306.92, 15640.14))

# Find the best prob_threshold
best_result <- grid_results[which.min(grid_results$rmse), ]
print(paste0("The best prob_threshold is: ", best_result$prob_threshold))
print(paste0("The corresponding RMSE is: ", best_result$rmse))

# Plot results for visualization
plot(grid_results$prob_threshold, grid_results$rmse, type = "b", 
     xlab = "Threshold Probability (prob_threshold)", 
     ylab = "RMSE", 
     main = "Grid Search Optimization of prob_threshold")

# Export graphics to rasters------------------------------------------------------------------
raster <- raster("Study_zone1.tif")                                             # New raster (correcting the provinces that produce)
raster_wgs84 <- projectRaster(raster,                                           
                              crs = CRS("+proj=longlat +datum=WGS84 +no_defs")) # Project raster to geographic coordinates
coords <- rasterToPoints(raster_wgs84)                                          # Generate a point for each pixel
# dim(coords)
# Assign infection values to each raster
for (i in 1:num_steps){
  raster_base <- raster_wgs84
  logical_vector <- infection_status_list[[i]][[10]]                            # Get the last logical vector of the current element (total infection at the end of each cycle)
  binary_values <- as.numeric(logical_vector)                                   # Convert the logical vector to 0 and 1 values
  puntos <- coords                                                              # Create a copy of the coordinates
  puntos[,3] <- binary_values                                                   # Assign values to the coordinates dataframe
  puntos_vect <- terra::vect(puntos[, c("x", "y")],                             # Convert to points
                             atts = puntos[, "Study_zone1", drop = FALSE], 
                             crs = as.character(crs(raster_base)))
  puntos_vect <- as(puntos_vect, "Spatial")
  new_raster <- rasterize(puntos_vect, raster_base, field = "Study_zone1")
  writeRaster(new_raster, filename = paste0("Infection_", i, ".tif"),           # Save raster
              overwrite = TRUE)
}





