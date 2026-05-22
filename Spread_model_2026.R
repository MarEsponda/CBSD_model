# ==============================================================================
# Project: Cassava Brown Streak Disease (CBSD) Epidemiological Spread Model
# Language: R
# Author: Maria del Mar Esponda
# Date: May 2026
# Organization: Alliance of Bioversity International and CIAT
# Description: This script processes spatial data from SPAM 2020 to delimit 
#              the physical cassava area, perform spatial downscaling, and 
#              generate a high-resolution distance matrix for CBSD 
#              epidemiological connectivity and spread simulation.
#              Includes the "Border Control" Policy Scenario at 75% effectiveness.
# ==============================================================================

# ---- LOAD DEPENDENCIES ----
library(sf)             # Spatial vector data manipulation
library(terra)          # Modern raster data analysis
library(raster)         # Prior raster data analysis (legacy support)
library(geosphere)      # Geospatial distance calculations
library(dplyr)          # Data manipulation and plumbing
library(units)          # Advanced measurement unit management
library(bigmemory)      # Out-of-core matrix management for large datasets
library(bench)          # High-precision process benchmarking
library(Metrics)        # Statistical error metrics
library(xlsx)           # Excel format handling
library(writexl)        # Light-weight Excel writing
library(readxl)         # Excel data ingestion
library(exactextractr)  # Zonal statistics extraction
library(drc)            # Dose-Response curve modeling

# ---- PORTABLE DIRECTORY SETUP ----
# Modify this path to point to your local root directory
base_dir <- "path/to/your/workspace"

input_dir     <- file.path(base_dir, "data", "input")
processed_dir <- file.path(base_dir, "data", "processed")
output_dir    <- file.path(base_dir, "output")

# Ensure environment sub-folders exist
if(!dir.exists(input_dir))     dir.create(input_dir, recursive = TRUE)
if(!dir.exists(processed_dir)) dir.create(processed_dir, recursive = TRUE)
if(!dir.exists(output_dir))    dir.create(output_dir, recursive = TRUE)


# Spatial Analysis & Downscaling==============================================================================
# Description: This script processes spatial data from SPAM 2020 to delimit 
#              the physical cassava area, perform spatial downscaling, and 
#              generate a high-resolution distance matrix for CBSD 
#              epidemiological connectivity and spread simulation.

setwd(input_dir)

# Load SPAM 2020 global physical area raster for Cassava (H_CASS_A)
cassava_area <- rast("spam2020_v1r0_global_H_CASS_A.tif") 

# Load study zone boundary polygon (Market Intelligence Unit)
study_zone <- st_read("Study_zonep.shp")

# Align Coordinate Reference Systems (CRS) to match the base raster
study_zone_wgs84 <- st_transform(study_zone, crs(cassava_area))

# Spatial Clipping and Masking
# Crop and mask the global dataset to the target study area boundaries
cassava_area <- crop(cassava_area, study_zone_wgs84) 
cassava_area <- mask(cassava_area, study_zone_wgs84) 

# Boundary Refinement by Country
# Ingest Administrative Boundaries (Level 0) to mask out non-production pixels
country_boundaries <- st_read("ADM0.shp")
country_rasterized <- rasterize(country_boundaries, cassava_area, field = "value")

# Apply custom logic threshold: Keep pixels within boundaries with >= 100 physical hectares
refined_area <- c(cassava_area, country_rasterized)
refined_area <- lapp(refined_area, fun = function(val, value) {
  ifelse(is.na(value) | val < 100, NA, val)
})

# Spatial Disaggregation (Downscaling)
# Subdivide each 9x9 km cell into a 3x3 grid of 3x3 km sub-cells (factor = 3)
raster_disagg <- disaggregate(refined_area, fact = 3) 

# Proportionally distribute the physical area values across the 9 new sub-cells
final_raster <- raster_disagg / 9 

# Generate a binary mask layout for the finalized study zone
study_zone_mask <- final_raster * 0 + 1

# Coordinate Extraction
# Extract geographic centroids for each production pixel
pixel_coords <- as.data.frame(rasterToPoints(final_raster))[, 1:2]
write.csv(pixel_coords, file = file.path(processed_dir, "coords1.csv"), row.names = FALSE)

# Big Matrix Initialization for Cohort Routing
# Initialize a file-backed big matrix to handle memory-intensive distance arrays
setwd(processed_dir)
total_pixels <- 40248
backing_file <- "dist.bin"
descriptor_file <- "dist.desc"

filebacked_block <- big.matrix(
  nrow = total_pixels, 
  ncol = total_pixels, 
  type = "double", 
  backingfile = backing_file, 
  descriptorfile = descriptor_file
)

distance_km_matrix <- attach.big.matrix("dist.desc")

# Geodesic Distance Computation (Haversine)
# Loop to compute the great-circle distance between all pixel centroids (for CBSD spread matrix)
for (i in 1:total_pixels) {
  distance_km_matrix[i, ] <- distm(pixel_coords[i, ], pixel_coords, fun = distHaversine)
  
  # Progress tracking output
  if (i %% 500 == 0) {
    print(paste("Processing infection/connectivity vectors for iteration:", i, "of", total_pixels))
  }
}


# Geodesic Distance Matrix Transformation==============================================================================
# Description: Converts the geodesic distance matrix from meters to kilometers 
#              and computes the log-transformed distances for the CBSD epidemiological 
#              spread kernel, optimizing memory for 1.6 billion data points.
setwd(processed_dir)

# Attach the file-backed distance matrix (originally computed in meters)
distance_meters_matrix <- attach.big.matrix("dist.desc")
total_pixels <- 40248

# Initialize Big Matrix for Kilometer Conversion
# Define file-backed structures to prevent memory overflow (RAM)
backing_file_km    <- "distanciaskm.bin"
descriptor_file_km <- "distanciaskm.desc"

filebacked_block_km <- big.matrix(
  nrow = total_pixels, 
  ncol = total_pixels, 
  type = "double", 
  backingfile = backing_file_km, 
  descriptorfile = descriptor_file_km
)

# Attach the kilometer matrix framework
distance_km_matrix <- attach.big.matrix("distanciaskm.desc")

# Convert Meters to Kilometers (Batch Processing)
# Scale down distances and round to 3 decimal places (meter precision)
for (i in 1:total_pixels) {
  distance_km_matrix[i, ] <- round((distance_meters_matrix[i, ] / 1000), 3)
  
  # Progress tracking output every 1000 iterations
  if (i %% 1000 == 0) {
    print(paste("Converted", i, "rows out of", total_pixels, "to kilometers."))
  }
}

# Logarithmic Transformation and Infinity Filtering
# Vectorize the matrix into a single array for epidemiological kernel fitting
distance_vector <- as.vector(distance_km_matrix[, ])

# Apply logarithmic transformation 
# Note: Pixels evaluated against themselves (0 km) will naturally generate -Inf
log_distance_vector <- log(distance_vector)

# Highly optimized vector filter to isolate finite parameters (removes -Inf, Inf, and NA)
# Replaces time-intensive loops with high-performance C++ level indexing
finite_log_distances <- log_distance_vector[is.finite(log_distance_vector)]

# Verification
# Summary printout for quality assurance of the spread parameters
print("Logarithmic distance vector successfully processed and filtered:")
print(summary(finite_log_distances))


# IPD Matrix Construction==============================================================================
# Description: Generates the contact/probability matrix based on the Inverse 
#              Power Law Dispersal Kernel. This block presents the baseline scenario 
#              where k = 2. In extended experimental runs, this parameter varies 
#              from k = 1.75 to k = 4 to calibrate the CBSD spread dynamics.
setwd(processed_dir)

# Re-attach the base kilometer distance matrix generated in Block 2
distance_km_matrix <- attach.big.matrix("distanciaskm.desc")
total_pixels <- 40248

# Initialize File-Backed Probability Matrix
# Define specific binary structures to manage the power-law kernel arrays without RAM exhaustion
backing_file_ipl    <- "distpi2.bin"
descriptor_file_ipl <- "distpi2.desc"

filebacked_block_ipl <- big.matrix(
  nrow = total_pixels, 
  ncol = total_pixels, 
  type = "double", 
  backingfile = backing_file_ipl, 
  descriptorfile = descriptor_file_ipl
)

# Attach the Inverse Power Law matrix framework
prob_matrix_k2 <- attach.big.matrix("distpi2.desc")

# Unified Scaling, Diagonal Fix, and Kernel Computation
# Optimization parameter setup (Baseline Power Law Exponent)
k_exponent <- 2

# Unified pipeline loop to minimize CPU context switching and loop overhead
for (i in 1:total_pixels) {
  # Step A: Ingest kilometer distances for row 'i'
  row_distances <- distance_km_matrix[i, ]
  
  # Step B: Set self-distance diagonal element to 1 to eliminate division-by-zero errors (Inf/NaN prevention)
  row_distances[i] <- 1
  
  # Step C: Compute Inverse Power Law equation: f(d) = 1 / (d^k) and round to 4 decimals
  prob_matrix_k2[i, ] <- round((1 / (row_distances ^ k_exponent)), 4)
  
  # Performance tracking console output every 2000 iterations
  if (i %% 2000 == 0) {
    print(paste("Kernel Matrix Pipeline: Processed row", i, "of", total_pixels, "for k =", k_exponent))
  }
}

# Quality Assurance Verification
print("Inverse Power Law kernel probability matrix initialization successful.")


# Epidemiological Metrics Curve Fitting==============================================================================
# Description: Ingests pre-calculated historical annual infection metrics 
#              and fits non-linear regression models (Log-Logistic & NLS) 
#              to define the CBSD infection trajectory limits.
# Note on Dispersal Thresholds: According to Delaquis (2018), trade events are 
# constrained by distance thresholds (<= 350 km), rejecting looser thresholds 
# (> 500 km) as suggested in alternative models (Andersen, 2024).
setwd(processed_dir)

# Ingest Pre-Calculated Epidemiological Data
# Load final consolidated metrics (Raw vector inputs are kept confidential)
# The file 'extraction_metrics.csv' contains columns: Year and Total_Pixels
extraction_metrics_df <- read.csv("extraction_metrics.csv")
print(extraction_metrics_df)

# Non-Linear Epidemiologic Model Trajectory Fitting
# Ingest base spatial framework to extract target carrying capacity bounds
study_zone_envelope <- rast(file.path(processed_dir, "Area_resampled100.tif"))
carrying_capacity   <- global(study_zone_envelope, fun = "sum", na.rm = TRUE)$sum

# Structure time-series vector and clean target infection parameters (13 iterations)
regression_data <- data.frame(
  time_steps = rep(1:13),
  infection  = c(0, 0, 0, 0, 0, 8939.733, 59453.762, 108232.438, 658075.875, 669253.688, 938931.375, 940197.438, 940907.188)
)

# A. Three-Parameter Log-Logistic Fitting (Fixed upper asymptote)
log_logistic_fit <- drm(infection ~ time_steps, data = regression_data, fct = LL.3(fixed = c(NA, carrying_capacity, NA)))
summary(log_logistic_fit)

# B. Non-Linear Least Squares (NLS) Sigmoidal Curve Optimization
nls_sigmoidal_fit <- nls(
  formula = infection ~ carrying_capacity / (1 + (a / time_steps)^b),
  data    = regression_data,
  start   = list(a = 11, b = 6),
  control = list(maxiter = 1000, warnOnly = TRUE)
)
summary(nls_sigmoidal_fit)
coef(nls_sigmoidal_fit)


# Initial Infection Seeding==============================================================================
# Description: Establishes the initial geographic outbreak nodes for the CBSD
#              simulation. Based on literature evidence from the May 2015 
#              outbreak in Ratanakiri, Cambodia, the script isolates the province 
#              boundaries and samples random baseline production pixels as 
#              epidemiologic seeds.
# Reference: http://dx.doi.org/10.1094/PDIS-10-15-1228-PDN
setwd(processed_dir)

# Ingest and Convert Coordinates Layers
# Load geographic pixel centroids generated in Block 1
study_zone_coords <- read.csv("coords1.csv")

# Project point data coordinate matrix to Spatial Features (sf) framework (WGS84)
coords_sf <- st_as_sf(study_zone_coords, coords = c("x", "y"), crs = 4326)

# Isolate Initial Outbreak Boundaries (Ratanakiri)
# Load historical infection ground-zero vector polygon
ratanak           <- st_read(file.path(input_dir, "Ratanakkiri.shp"))
ratanak           <- st_transform(ratanak, crs(coords_sf))

# Geo-Intersection and Pixel Identification
# Execute geometric intersection to identify production pixels inside Ratanakiri
inside_polygon         <- st_intersects(coords_sf, ratanak, sparse = FALSE)
inicial                <- which(inside_polygon == TRUE)

# Random Reproducible Seeding (Patient Zero Selection)
# Set seed configuration to guarantee exact script reproducibility
set.seed(123)

# Randomly sample 2 baseline production nodes to serve as infection sources
initial_infected_nodes <- sample(inicial, 2, replace = FALSE)


# Epidemiological Simulation Model for CBSD==============================================================================
# Description: Executes a dual-component epidemiological simulation model for CBSD.
#              Calculates long-distance trade/spore events via power-law kernels
#              and handles localized short-distance expansion driven by whitefly 
#              vectors (Bemisia tabaci) using monthly nested stochastic steps.
setwd(processed_dir)

# Configuration Parameters and Algorithm Initialization
# Simulation thresholds and execution parameters
probability_threshold <- 0.60    # Local infection success probability
total_pixels          <- 40248
macro_cycles          <- 15      # Number of macro time-steps (iterations)
monthly_substeps      <- 3       # Nested internal micro-steps per cycle
whitefly_range_km     <- 3       # Local vector active flight distance boundary (SDD threshold)

# Attach big matrices from previous memory-mapped blocks
kernel_matrix_ipl  <- attach.big.matrix("distpi1.98.desc") # Calibrated power-law matrix (k=1.98)
distance_km_matrix <- attach.big.matrix("distanciaskm.desc")

# Ingest uniform spatial framework raster using terra
baseline_raster    <- rast("Area_resampled100.tif")
study_zone_grid    <- rast("studyzone_100.tif")
coords             <- rasterToPoints(raster(study_zone_grid))

# Initialize structures to capture epidemiologic history arrays
infection_history_log  <- vector("list", macro_cycles)
aggregate_infected_sum <- rep(NA, macro_cycles)

# Seed Context and Baseline State Injection
set.seed(123) # Guarantee exact stochastic simulation replication

# Vector initialization: 0 = Susceptible, 1 = Infected/Infectious
infection_status <- rep(0, total_pixels)

# Inject patient zero nodes selected from the Ratanakiri baseline trigger (Block 5)
infection_status[initial_infected_nodes] <- 1

print("Epidemiological simulation baseline matrix primed. Executing spread engine...")

# Spatio-Temporal Dual Simulation Engine Loop
for (cycle in 1:macro_cycles) {
  ldd_new_infections <- rep(0, total_pixels)
  active_infected_indices <- which(infection_status == 1)
  
  # COMPONENT A: LONG-DISTANCE DISPERSAL (LDD) - Spatial Power-Law Kernel
  for (i in active_infected_indices) {
    transmission_probabilities <- kernel_matrix_ipl[i, ]
    stochastic_draw            <- runif(total_pixels) < transmission_probabilities
    ldd_new_infections         <- ldd_new_infections | stochastic_draw
  }
  infection_status <- infection_status | ldd_new_infections
  
  print(paste("--- Macro Cycle:", cycle, "[LDD Phase Completed] ---"))
  print(table(infection_status))
  
  # COMPONENT B: SHORT-DISTANCE DISPERSAL (SDD) - Local Whitefly Flight Vector
  monthly_timeline_cache <- vector("list", monthly_substeps)
  for (month in 1:monthly_substeps) {
    sdd_new_infections <- rep(0, total_pixels)
    current_infected_indices <- which(infection_status == 1)
    
    for (i in current_infected_indices) {
      proximal_neighbors <- distance_km_matrix[i, ] < whitefly_range_km
      total_proximal     <- sum(proximal_neighbors)
      
      if (total_proximal > 0) {
        local_stochastic_draw <- runif(total_proximal) < probability_threshold
        sdd_new_infections[proximal_neighbors] <- sdd_new_infections[proximal_neighbors] | local_stochastic_draw
      }
    }
    infection_status <- infection_status | sdd_new_infections
    monthly_timeline_cache[[month]] <- infection_status
    print(paste0("Propagaci?n mosca, ciclo ", cycle, ", mes ", month))
  }
  infection_history_log[[cycle]] <- monthly_timeline_cache
  print(paste("Fin del Ciclo:", cycle))
}


# Spatial Raster Export Pipeline and Validation==============================================================================
# Description: Maps the logical matrix states from the CBSD stochastic simulation 
#              into geographic surfaces, intersections them with SPAM 2020 
#              cassava physical areas, exports time-series rasters for animations, 
#              and evaluates model fit performance via Root Mean Square Error (RMSE).
setwd(processed_dir)

# Setup Environmental Layers and Cache Storage
raster_base <- raster("studyzone_100.tif")
coords      <- rasterToPoints(raster_base)
area_total_infectada <- numeric(macro_cycles)
suma_pix             <- rep(NA, macro_cycles)

# Temporal Rasterization and Surface Intersection Loop
for (i in 1:macro_cycles){
  logical_vector <- infection_history_log[[i]][[3]]
  suma_pix[i]    <- sum(logical_vector)
  binary_values  <- as.numeric(logical_vector)
  puntos         <- coords
  puntos[, 3]    <- binary_values
  
  puntos_vect <- terra::vect(puntos[, c("x", "y")], 
                             atts = puntos[, "studyzone_100", drop = FALSE], 
                             crs = as.character(crs(raster_base)))
  puntos_vect <- as(puntos_vect, "Spatial")
  new_raster  <- rasterize(puntos_vect, raster_base, field = "studyzone_100")
  new_raster  <- rast(new_raster)
  
  # Superimpose with cultivated area raster
  raster_cultivado_infectado <- new_raster * rast(baseline_raster)
  raster_cultivado_infectado[raster_cultivado_infectado == 0] <- NA
  area_total_infectada[i] <- sum(values(raster_cultivado_infectado), na.rm = TRUE)
  
  # Save resulting raster
  writeRaster(raster_cultivado_infectado, 
              filename = file.path(output_dir, paste0("Ainfecc_3ciclos", i, ".tif")), 
              overwrite = TRUE)
}

# Model Fit and Statistical Performance Metrics
validation_metrics_df <- data.frame(
  Month            = seq(from = 10, to = (macro_cycles * 10), by = 10),
  Projected_Area_Y = c(64, 1129, 6034, 19715, 48909, 101175, 183035, 297501, 442022, 608229, 784004, 956917, 1117183, 1258979, 1380192),
  Simulated_Area_X = round(area_total_infectada, 0),
  Infected_Pixels  = suma_pix
)
validation_metrics_df <- validation_metrics_df %>% mutate(Model_RMSE = rmse(Projected_Area_Y, Simulated_Area_X))
print("====== CBSD EPIDEMIOLOGICAL MODEL VALIDATION COMPLETE ======")
print(validation_metrics_df)


# Optimization via Grid Search==============================================================================
# Description: Searches for the optimal 'prob_threshold' value that yields the 
#              best fit against model targets by iterating through a grid of 
#              probabilities and evaluating performance using RMSE.
setwd(processed_dir)

# Grid Search Parameter Space Setup
prob_values  <- seq(0.1, 1.0, by = 0.01)
grid_results <- data.frame(prob_threshold = prob_values, rmse = NA)

# Global simulation parameters and big matrices load
distpi2     <- attach.big.matrix("distpi1.95.desc")
distkm      <- attach.big.matrix("distanciaskm.desc")
raster_area <- rast("Area_resampled100.tif")
raster_base <- raster("studyzone_100.tif")
coords      <- rasterToPoints(raster_base)

# Simulation and Evaluation Function Architecture
simulate_and_evaluate <- function(prob_threshold) {
  set.seed(123) 
  infection_status <- rep(0, pix)
  infection_status[initial_infected_nodes] <- 1
  
  for (step in 1:num_steps) {
    new_infections <- rep(0, pix)
    for (i in which(infection_status == 1)) {
      prob_infection <- distpi2[i, ]
      new_infections <- new_infections | (runif(pix) < prob_infection)
    }
    infection_status <- infection_status | new_infections
    
    monthly_infection_status <- vector("list", monthly_steps)
    for (month in 1:monthly_steps) {
      new_infections_dist <- rep(0, pix)
      for (i in which(infection_status == 1)) {
        within_distance <- distkm[i, ] < whithelfy
        randomized <- runif(sum(within_distance)) < prob_threshold
        new_infections_dist[within_distance] <- new_infections_dist[within_distance] | randomized
      }
      infection_status <- infection_status | new_infections_dist
      monthly_infection_status[[month]] <- infection_status
    }
    infection_status_list[[step]] <- monthly_infection_status
  }
  
  for (i in 1:num_steps) {
    logical_vector <- infection_status_list[[i]][[3]]
    suma_pix[i]    <- sum(infection_status_list[[i]][[3]])
    binary_values  <- as.numeric(logical_vector)
    puntos         <- coords
    puntos[, 3]    <- binary_values
    
    puntos_vect <- terra::vect(puntos[, c("x", "y")], 
                               atts = puntos[, "studyzone_100", drop = FALSE], 
                               crs = as.character(crs(raster_base)))
    puntos_vect <- as(puntos_vect, "Spatial")
    new_raster  <- rasterize(puntos_vect, raster_base, field = "studyzone_100")
    new_raster  <- rast(new_raster)
    
    raster_cultivado_infectado <- new_raster * raster_area
    area_total_infectada[i]    <- sum(values(raster_cultivado_infectado), na.rm = TRUE)
  }
  
  modelo <- data.frame(
    mes      = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150), 
    areaproy = c(64, 1129, 6034, 19715, 48909, 101175, 183035, 297501, 442022, 608229, 784004, 956917, 1117183, 1258979, 1380192),
    area     = round(area_total_infectada, 0),
    pix      = suma_pix
  )
  modelo <- mutate(modelo, rmse = rmse(modelo$areaproy, modelo$area))
  return(modelo)
}

# Run Grid Search execution
for (i in seq_along(prob_values)) {
  resultado[[i]]       <- simulate_and_evaluate(prob_values[i])
  grid_results$area[i] <- resultado[[i]][15, 3]
  grid_results$pix[i]  <- resultado[[i]][15, 4]
  grid_results$rmse[i] <- resultado[[i]][15, 5]
}

best_result <- grid_results[which.min(grid_results$rmse), ]
print(paste0("El mejor prob_threshold es: ", best_result$prob_threshold))

# Export Diagnostic Optimization Chart
plot(
  grid_results$prob_threshold, grid_results$rmse, 
  type = "b", 
  xlab = "Probabilidad umbral (prob_threshold)", 
  ylab = "RMSE", 
  main = "Optimizacion de prob_threshold con Grid Search"
)


# Regional Impact Assessment - Losses by Province==============================================================================
# Description: Quantifies the temporal biophysical footprint of CBSD per province 
#              (ADM1). Extracts zonal sums for both infected cassava production 
#              (metric tons) and infected cultivation area (hectares) across 
#              the 15 macro simulation cycles.

# PART 1: IMPACT LOSSES BY PROVINCE (TONS)
setwd(processed_dir)

# Load administrative boundaries and baseline SPAM 2020 matrices
prov            <- st_read(file.path(input_dir, "ADM1_prod.shp"))
Harea           <- raster("Area_resampled100.tif")
Pton            <- rast("Prod_resampled100.tif")
z_estudio_wgs84 <- st_transform(st_read(file.path(input_dir, "Study_zonep.shp")), crs(Harea))

# Mask Intersection: Simulation Maps x Production Matrix
raster_list <- list()
for (i in 1:15) {
  raster_list[[i]] <- rast(file.path(output_dir, paste0("Ainfecc_3ciclos", i, ".tif")))
  r <- raster_list[[i]]
  values(r)[values(r) > 0] <- 1
  raster_list[[i]] <- r * Pton
}
raster_list[[16]] <- Pton

# Extraction of Zonal Production Metrics
resultados           <- lapply(raster_list, function(capa) { exact_extract(capa, prov, "sum") })
resultados           <- as.data.frame(resultados)
colnames(resultados) <- c(paste0("year", 1:15), "totalp")
resultados           <- round(resultados, 1)

# Compilation and Data Export
compendio <- as.data.frame(prov[, c(3, 4)]) %>% st_drop_geometry()
compendio <- cbind(compendio, resultados)[, -3]
write.xlsx(compendio, file.path(output_dir, "resultados_3ciclos_p.xlsx"))

# PART 2: IMPACT LOSSES BY PROVINCE (HECTARES)
Harea <- rast("Area_resampled100.tif")

# Mask Intersection: Simulation Maps x Harvested Area Matrix
raster_list <- list()
for (i in 1:15) {
  raster_list[[i]] <- rast(file.path(output_dir, paste0("Ainfecc_3ciclos", i, ".tif")))
  r <- raster_list[[i]]
  values(r)[values(r) > 0] <- 1
  raster_list[[i]] <- r * Harea
}
raster_list[[16]] <- Harea

# Extraction of Zonal Physical Area Metrics
resultados           <- lapply(raster_list, function(capa) { exact_extract(capa, prov, "sum") })
resultados           <- as.data.frame(resultados)
colnames(resultados) <- c(paste0("year", 1:15), "totalha")
resultados           <- round(resultados, 1)

# Compilation and Data Export
compendio <- as.data.frame(prov[, c(3, 4)]) %>% st_drop_geometry()
compendio <- cbind(compendio, resultados)[, -3]
write.xlsx(compendio, file.path(output_dir, "res_3ciclos_a.xlsx"))


# Policy Experimentation - "Border Control" Scenario (75% Effectiveness)==============================================================================
# Description: Generates a modified distance weight matrix simulating a 75% 
#              effective quarantine barrier between Zone 1 and Zone 2, runs the 
#              stochastic spread simulation, and extracts regional impacts.
setwd(output_dir)

# Study Zone Definition and GADM Administrative Boundaries
{
  paises      <- c("Vietnam", "Laos", "Thailand", "Cambodia")
  lista_munis <- lapply(paises, function(p) { gadm(country = p, level = 0, path = file.path(input_dir, "datos_gadm")) })
  z_estudio   <- do.call(rbind, lista_munis)
  writeVector(z_estudio, file.path(input_dir, "Cuatro_paises.shp"))
}

# Coupling Coordinates with Administrative Boundaries
setwd(processed_dir)
paises    <- st_read(file.path(input_dir, "Cuatro_paises.shp"))
coords    <- read.csv("coords1.csv")  
coords_sf <- st_as_sf(coords, coords = c("x", "y"), crs = 4326)

# Verify Coordinate Reference Systems alignment
st_crs(coords_sf) == st_crs(paises)

# Spatial join to assign country/zone values to grid centroids
coords_zonas <- st_join(coords_sf, paises["Value"])

# Downward imputation loop to fill unassigned border cells (118 edge pixels)
for (i in 2:nrow(coords_zonas)) {
  if (is.na(coords_zonas$Value[i])) coords_zonas$Value[i] <- coords_zonas$Value[i - 1]
}

# Isolate structural indices corresponding to Zone 1 and Zone 2
zona_vector <- coords_zonas$Value 
zona1_ids   <- which(zona_vector == 1) 
zona2_ids   <- which(zona_vector == 2) 

# Generate Binary Border Control Policy Weight Matrix (75% Effectiveness)
control_mat <- big.matrix(nrow = pix, ncol = pix, type = "double",
                          backingfile = "bordercontrol_75.bin", descriptorfile = "bordercontrol_75.desc")
control_mat[,] <- 1 # Initialize background matrix with uniform values of 1

# Apply a strict 75% reduction (factor = 0.25) to routes connecting Zone 1 to Zone 2
for (i in zona1_ids) {
  control_mat[i, zona2_ids] <- 0.25
}

# Hadamard Matrix Product for k=1.98 Power-Law Kernel
distpi1.98        <- attach.big.matrix("distpi1.98.desc")
distpi1.98_border <- big.matrix(nrow = pix, ncol = pix, type = "double",
                                backingfile = "distpi1.98_BC75.bin", descriptorfile = "distpi1.98_BC75.desc")

for(i in 1:pix){
  distpi1.98_border[i, ] <- distpi1.98[i, ] * control_mat[i, ]
}

# Core Epidemiological Simulation Run (75% Border Control)
set.seed(123)
infection_status                     <- rep(0, pix)
infection_status[as.integer(inicial[1:2])] <- 1 

distpi2     <- attach.big.matrix("distpi1.98_BC75.desc")
distkm      <- attach.big.matrix("distanciaskm.desc")
raster_area <- rast("Area_resampled100.tif")
raster_base <- raster("studyzone_100.tif")
coords      <- rasterToPoints(raster_base)

for (step in 1:num_steps) {
  new_infections <- rep(0, pix)
  for (i in which(infection_status == 1)) { 
    prob_infection <- distpi2[i, ] 
    new_infections = new_infections | (runif(pix) < prob_infection)
  }
  infection_status <- infection_status | new_infections
  
  # Component B: Short-Distance Local Vector Expansion
  monthly_timeline_cache <- vector("list", monthly_substeps)
  for (month in 1:monthly_substeps) { 
    new_infections_dist <- rep(0, pix)
    for (i in which(infection_status == 1)) { 
      within_distance <- distkm[i, ] < whithelfy 
      n_within        <- sum(within_distance)
      if (n_within > 0) {
        randomized                             <- runif(n_within) < prob_threshold
        new_infections_dist[within_distance] <- new_infections_dist[within_distance] | randomized
      }
    }
    infection_status                <- infection_status | new_infections_dist
    monthly_timeline_cache[[month]] <- infection_status
  }
  infection_status_list[[step]] <- monthly_timeline_cache
}

# Post-Processing: Spatial Output Mapping
for (i in 1:num_steps){
  logical_vector <- infection_status_list[[i]][[3]]
  suma_pix[i]    <- sum(logical_vector)
  binary_values  <- as.numeric(logical_vector)
  puntos         <- coords
  puntos[, 3]    <- binary_values
  
  puntos_vect <- terra::vect(puntos[, c("x", "y")], atts = puntos[, "studyzone_100", drop = FALSE], crs = as.character(crs(raster_base)))
  puntos_vect <- as(puntos_vect, "Spatial")
  new_raster  <- rasterize(puntos_vect, raster_base, field = "studyzone_100")
  new_raster  <- rast(new_raster)
  
  raster_cultivado_infectado <- new_raster * raster_area
  raster_cultivado_infectado[raster_cultivado_infectado == 0] <- NA
  area_total_infectada[i] <- sum(values(raster_cultivado_infectado), na.rm = TRUE)
  
  writeRaster(raster_cultivado_infectado, filename = file.path(output_dir, paste0("Ainfecc75_3ciclos", i, ".tif")), overwrite = TRUE)
}

# Policy Metrics Extractions (Tons & Hectares)
prov  <- st_read(file.path(input_dir, "ADM1_prod.shp"))
Pton  <- rast("Prod_resampled100.tif")
Harea <- rast("Area_resampled100.tif")

# Regional Impact Extraction: Losses by Province (Tons)
raster_list <- list()
for (i in 1:15) {
  raster_list[[i]] <- rast(file.path(output_dir, paste0("Ainfecc75_3ciclos", i, ".tif")))
  r                <- raster_list[[i]]
  values(r)[values(r) > 0] <- 1
  raster_list[[i]] <- r * Pton
}
raster_list[[16]]    <- Pton
resultados           <- lapply(raster_list, function(capa) { exact_extract(capa, prov, "sum") })
compendio            <- cbind(as.data.frame(prov[, c(3, 4)]) %>% st_drop_geometry(), round(as.data.frame(resultados), 1))[, -3]
colnames(compendio)  <- c("ADM0_NAME", "ADM1_NAME", paste0("year", 1:15), "totalp")
write.xlsx(compendio, file.path(output_dir, "75_bc_p.xlsx"))

# Regional Impact Extraction: Losses by Province (Hectares)
raster_list <- list()
for (i in 1:15) {
  raster_list[[i]] <- rast(file.path(output_dir, paste0("Ainfecc75_3ciclos", i, ".tif")))
  r                <- raster_list[[i]]
  values(r)[values(r) > 0] <- 1
  raster_list[[i]] <- r * Harea
}
raster_list[[16]]    <- Harea
resultados           <- lapply(raster_list, function(capa) { exact_extract(capa, prov, "sum") })
compendio            <- cbind(as.data.frame(prov[, c(3, 4)]) %>% st_drop_geometry(), round(as.data.frame(resultados), 1))[, -3]
colnames(compendio)  <- c("ADM0_NAME", "ADM1_NAME", paste0("year", 1:15), "totalha")
write.xlsx(compendio, file.path(output_dir, "75_bc_ha.xlsx"))