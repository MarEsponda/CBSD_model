# ==============================================================================
# Project: CBSD Epidemiological Model - Comprehensive Sensitivity Analysis Suite
# Language: R
# Author: Maria del Mar Esponda
# Date: May 2026
# Organization: Alliance of Bioversity International and CIAT
# Description: Unified sensitivity suite for Cassava Brown Streak Disease (CBSD).
#              Part 1: Monte Carlo Spatial Seed Location Vulnerability.
#              Part 2: Local Whitefly Transmission Parameter Perturbation (+/- 10%).
#              Part 3: Pure Algorithm Stochastic Uncertainty (Random Seeds Variance).
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
library(readxl)         # Excel data ingestion
library(exactextractr)  # Zonal statistics extraction
library(drc)            # Dose-Response curve modeling
library(future.apply)   # High-performance parallelized processing loops
library(gstat)          # Geostatistical modeling and spatial interpolation
library(ggplot2)        # Data visualization and publication graphics
library(scales)         # Coordinate formatting scale tools

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


# PART 1: Spatial Seed Location Vulnerability Analysis==============================================================================
# Description: Evaluates structural variations in total impacted physical area 
#              by randomizing initial patient zero nodes across target production 
#              zones, running a 100-iteration stochastic cascade per country.
setwd(processed_dir)

# Ingest uniform spatial grid matrices
coords1   <- read.csv("coords1.csv")
coords_sf <- st_as_sf(coords1, coords = c("x", "y"), crs = 4326)

# Global simulation parameters
prob_threshold <- 0.6
pix            <- 40248
num_steps      <- 15    
monthly_steps  <- 3     
whithelfy      <- 3     

# Attach big memories and core frameworks
distpi2               <- attach.big.matrix("distpi1.98.desc")
distkm                <- attach.big.matrix("distanciaskm.desc")
raster_area           <- rast("Area_resampled100.tif")
raster_base           <- raster("studyzone_100.tif")
coords                <- rasterToPoints(raster_base)
infection_status_list <- vector("list", num_steps)
area_total_infectada  <- numeric(num_steps)
suma_pix              <- rep(NA, num_steps)

# Set country iteration scope for the current active session
n_iteraciones <- 100
pais_actual   <- "Vietnam" 

# Load and project local administrative vector envelope
pais           <- st_read(file.path(input_dir, paste0(pais_actual, ".shp")), quiet = TRUE)
pais           <- st_transform(pais, st_crs(coords_sf))
inside_polygon <- st_intersects(coords_sf, pais, sparse = FALSE)
inicial        <- which(inside_polygon == TRUE)

# Define Spatial Stochastic Sensitivity Wrapper Function
sensibilidad_espacial <- function(id_iter) {
  set.seed(id_iter + as.numeric(Sys.time())) 
  initial_infected <- sample(inicial, 1, replace = FALSE)
  coords_inicio    <- coords1[initial_infected, ]
  
  infection_status <- rep(0, pix)
  infection_status[initial_infected] <- 1
  resultados_anuales <- numeric(num_steps)
  
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
        n_within        <- sum(within_distance)
        if (n_within > 0) {
          randomized                             <- runif(n_within) < prob_threshold
          new_infections_dist[within_distance] <- new_infections_dist[within_distance] | randomized
        }
      }
      infection_status                <- infection_status | new_infections_dist
      monthly_infection_status[[month]] <- infection_status
    }
    infection_status_list[[step]] <- monthly_infection_status
    
    logical_vector <- infection_status_list[[step]][[3]]
    suma_pix[step] <- sum(logical_vector)
    binary_values  <- as.numeric(logical_vector)
    
    puntos      <- coords
    puntos[, 3] <- binary_values
    puntos_vect <- terra::vect(puntos[, c("x", "y")], atts = puntos[, "studyzone_100", drop = FALSE], crs = as.character(crs(raster_base)))
    puntos_vect <- as(puntos_vect, "Spatial")
    
    new_raster <- rasterize(puntos_vect, raster_base, field = "studyzone_100")
    new_raster <- rast(new_raster)
    
    raster_cultivado_infectado <- new_raster * raster_area
    raster_cultivado_infectado[raster_cultivado_infectado == 0] <- NA
    resultados_anuales[step] <- sum(values(raster_cultivado_infectado), na.rm = TRUE)
  }
  
  #To characterize the epidemic velocity, two distinct metrics were computed: (i) the period-specific relative growth rate (tasas), 
  #which captures the proportional marginal increase in infected area between consecutive cycles ($\Delta A_t / A_{t-1}$),
  #and (ii) the overall intrinsic growth rate (tasa_global_r), an exponential per-capita scaling coefficient ($r$) derived from the 
  #initial and terminal states of the simulated outbreak footprint over the total temporal horizon.
  
  tasas         <- c(0, diff(resultados_anuales) / resultados_anuales[-length(resultados_anuales)])
  tasa_global_r <- (log(resultados_anuales[num_steps]) - log(1)) / num_steps
  
  df_iter <- data.frame(
    id_iteracin    = id_iter,
    pais_inicio    = pais_actual,
    lat_inicio     = coords_inicio[3],
    long_inicio    = coords_inicio[2],
    ciclo          = 1:num_steps,
    area_infectada = round(resultados_anuales, 1),
    tasa_anual     = tasas,
    tasa_global_r  = tasa_global_r
  )
  return(df_iter)
}

print("Iniciando análisis de sensibilidad espacial...")
lista_resultados_esp      <- future_lapply(1:n_iteraciones, sensibilidad_espacial, future.seed = TRUE)
df_final_sensibilidad_esp <- do.call(rbind, lista_resultados_esp)
saveRDS(lista_resultados_esp, file.path(output_dir, paste0(pais_actual, "_500_a.rds")))


# Multi-Country Consolidation and Descriptive Statistics Summary==============================================================================
# Description: Ingests independent country RDS files to compile absolute loss statistics.
setwd(output_dir)

Thai <- readRDS("Thai_100.rds")
Camb <- readRDS("Cambodia_100.rds")
Lao  <- readRDS("Laos_100.rds")
Viet <- readRDS("Vietnam_100.rds")

sensib1 <- rbind(do.call(rbind, Thai), do.call(rbind, Camb), do.call(rbind, Lao), do.call(rbind, Viet))
sensib2 <- filter(sensib1, ciclo == 15)

sensib3 <- sensib1 %>%
  filter(ciclo == 15) %>%
  group_by(pais_inicio) %>%
  summarize(mean           = mean(area_infectada),
            desv_s         = sd(area_infectada),
            Mediana        = median(area_infectada),
            CV_porcentaje  = (sd(area_infectada) / mean(area_infectada)) * 100,
            Minimo         = min(area_infectada),
            Maximo         = max(area_infectada),
            P5             = quantile(area_infectada, 0.05),
            P95            = quantile(area_infectada, 0.95),
            Error_Estandar = sd(area_infectada) / sqrt(n()),
            rate           = mean(tasa_global_r))

write.xlsx(sensib2, "Sensi_outbreak.xlsx")


# Spatial Vector Generation & IDW Spatial Interpolation==============================================================================
# Description: Generates geostatistical surface continuum layers using Inverse Distance Weighting.
puntos_sf <- st_as_sf(sensib2, coords = c("long_inicio", "lat_inicio"), crs = 4326)
st_write(puntos_sf, "Sensibilidad_p.shp", delete_layer = TRUE)

shape_ref     <- st_read(file.path(input_dir, "ADM0.shp"))
grid_template <- terra::rast(raster_area)
modelo_idw    <- gstat(formula = area_infectada ~ 1, data = puntos_sf, nmax = 12, set = list(idp = 2.0))
raster_interp <- terra::interpolate(raster_area, modelo_idw)

r_area_mask <- raster_area
r_area_mask[r_area_mask <= 0] <- NA
raster_final <- mask(raster_interp, r_area_mask)


# Sample Size Asymptotic Convergence Evaluation==============================================================================
# Description: Evaluates spatial verification limits comparing 100 vs 500 nodes in Vietnam.
Viet_500_a <- readRDS("Viet_500_a.rds")
Viet_500_b <- readRDS("Viet_500_b.rds")
Viet_500_c <- readRDS("Viet_500_c.rds")
Viet_500_d <- readRDS("Viet_500_d.rds")
Viet_500_e <- readRDS("Viet_500_e.rds")

sensib1_conv <- rbind(do.call(rbind, Viet_500_a), do.call(rbind, Viet_500_b), 
                      do.call(rbind, Viet_500_c), do.call(rbind, Viet_500_d), do.call(rbind, Viet_500_e))

sensib3_conv <- sensib1_conv %>%
  filter(ciclo == 15) %>%
  summarize(mean = mean(area_infectada), desv_s = sd(area_infectada))

print("====== ASYMPTOTIC CONVERGENCE TESTS (VIET_500) ======")
print(sensib3_conv)


# PART 2: Local Transmission Parameter Perturbation Analysis==============================================================================
# Description: Runs a 100-iteration loop over a continuous uniform perturbation 
#              range (+/- 10%) over the whitefly transmission probability parameter.
setwd(processed_dir)

ratanak        <- st_read(file.path(input_dir, "Ratanakkiri.shp")) %>% st_transform(crs(coords_sf))
inside_polygon <- st_intersects(coords_sf, ratanak, sparse = FALSE)
inicial        <- which(inside_polygon == TRUE)

set.seed(123)
initial_infected <- sample(inicial, 2, replace = FALSE)
coords_inicio    <- coords1[initial_infected, ]

sensibilidad_param <- function(id_iter) {
  set.seed(id_iter + as.numeric(Sys.time())) 
  prob_threshold     <- runif(1, min = p_opt * 0.9, max = p_opt * 1.1)
  infection_status   <- rep(0, pix)
  infection_status[initial_infected] <- 1
  resultados_anuales <- numeric(num_steps)
  
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
        n_within        <- sum(within_distance)
        if (n_within > 0) {
          randomized                             <- runif(n_within) < prob_threshold
          new_infections_dist[within_distance] <- new_infections_dist[within_distance] | randomized
        }
      }
      infection_status                <- infection_status | new_infections_dist
      monthly_infection_status[[month]] <- infection_status
    }
    infection_status_list[[step]] <- monthly_infection_status
    
    logical_vector <- infection_status_list[[step]][[3]]
    suma_pix[step] <- sum(logical_vector)
    binary_values  <- as.numeric(logical_vector)
    puntos        <- coords
    puntos[, 3]   <- binary_values
    
    puntos_vect <- terra::vect(puntos[, c("x", "y")], atts = puntos[, "studyzone_100", drop = FALSE], crs = as.character(crs(raster_base)))
    puntos_vect <- as(puntos_vect, "Spatial")
    
    new_raster <- rasterize(puntos_vect, raster_base, field = "studyzone_100")
    new_raster <- rast(new_raster)
    
    raster_cultivado_infectado <- new_raster * raster_area
    raster_cultivado_infectado[raster_cultivado_infectado == 0] <- NA
    resultados_anuales[step] <- sum(values(raster_cultivado_infectado), na.rm = TRUE)
  }
  
  tasas         <- c(0, diff(resultados_anuales) / resultados_anuales[-length(resultados_anuales)])
  tasa_global_r <- (log(resultados_anuales[num_steps]) - log(1)) / num_steps
  
  df_iter <- data.frame(
    id_iteracin    = id_iter,
    w_factor       = prob_threshold, 
    ciclo          = 1:num_steps,
    area_infectada = round(resultados_anuales, 1),
    tasa_anual     = tasas,
    tasa_global_r  = tasa_global_r
  )
  return(df_iter)
}

print("Iniciando análisis de sensibilidad paramétrica...")
lista_resultados_param      <- future_lapply(1:n_iteraciones, sensibilidad_param, future.seed = TRUE)
df_final_sensibilidad_param <- do.call(rbind, lista_resultados_param)
saveRDS(lista_resultados_param, file.path(output_dir, "P_thres.rds"))

# Statistical Reporting & Plots (Part 2)
setwd(output_dir)
sensib2_param <- filter(df_final_sensibilidad_param, ciclo == 15) %>% mutate(A_inf_kha = round(area_infectada / 1000, 0))
write.xlsx(sensib2_param, "P_thres.xlsx")

# Run non-parametric validation correlation models
cor_spearman <- cor.test(sensib2_param$w_factor, sensib2_param$area_infectada, method = "spearman")

# Plot Distribution curves
grafico_param <- ggplot(sensib2_param, aes(x = A_inf_kha)) +
  geom_density(fill = "#4CB391", color = "#2E6652", alpha = 0.6) +
  geom_vline(aes(xintercept = mean(A_inf_kha)), color = "red", linetype = "dashed", size = 0.8) +
  scale_x_continuous(labels = comma, n.breaks = 10) +
  labs(title = "Probability Density of Final Infected Area", x = "Total Infected Area (1,000 Ha)", y = "Probability Density") +
  theme_bw()


# PART 3: Pure Algorithm Stochastic Uncertainty Analysis==============================================================================
# Description: Hold spatial seed anchors and parameters strictly fixed. Evaluates 
#              the purely stochastic variance (random seed trajectories) across 100 loops.
setwd(processed_dir)

set.seed(123)
initial_infected <- sample(inicial, 2, replace = FALSE)  
coords_inicio    <- coords1[initial_infected, ]

sensibilidad_stoch <- function(id_iter) {
  set.seed(id_iter + as.numeric(Sys.time())) 
  infection_status   <- rep(0, pix)
  infection_status[initial_infected] <- 1
  resultados_anuales                 <- numeric(num_steps)
  
  for (step in 1:num_steps) {
    new_infections <- rep(0, pix)
    for (i in which(infection_status == 1)) { 
      prob_infection <- distpi2[i, ] 
      new_infections = new_infections | (runif(pix) < prob_infection)
    }
    infection_status <- infection_status | new_infections
    
    monthly_timeline_cache <- vector("list", monthly_steps)
    for (month in 1:monthly_steps) { 
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
    
    logical_vector <- infection_status_list[[step]][[3]]
    suma_pix[step] <- sum(logical_vector)
    binary_values  <- as.numeric(logical_vector)
    puntos         <- coords
    puntos[, 3]    <- binary_values
    
    puntos_vect <- terra::vect(puntos[, c("x", "y")], atts = puntos[, "studyzone_100", drop = FALSE], crs = as.character(crs(raster_base)))
    puntos_vect <- as(puntos_vect, "Spatial")
    new_raster  <- rasterize(puntos_vect, raster_base, field = "studyzone_100")
    new_raster  <- rast(new_raster)
    
    raster_cultivado_infectado <- new_raster * raster_area
    raster_cultivado_infectado[raster_cultivado_infectado == 0] <- NA
    resultados_anuales[step] <- sum(values(raster_cultivado_infectado), na.rm = TRUE)
  }
  
  tasas         <- c(0, diff(resultados_anuales) / resultados_anuales[-length(resultados_anuales)])
  tasa_global_r <- (log(resultados_anuales[num_steps]) - log(1)) / num_steps
  
  df_iter <- data.frame(
    id_iteracin    = id_iter,
    lat_inicio     = coords_inicio[1, 3], 
    long_inicio    = coords_inicio[1, 2], 
    ciclo          = 1:num_steps,
    area_infectada = round(resultados_anuales, 1),
    tasa_anual     = tasas,
    tasa_global_r  = tasa_global_r
  )
  return(df_iter)
}

print("Iniciando análisis de sensibilidad estocástica pura...")
lista_resultados_stoch      <- future_lapply(1:n_iteraciones, sensibilidad_stoch, future.seed = TRUE)
df_final_sensibilidad_stoch <- do.call(rbind, lista_resultados_stoch)
saveRDS(lista_resultados_stoch, file.path(output_dir, "Random_sensi.rds"))

# Multi-Temporal Confidence Profile Generation (Fan Chart Data)
setwd(output_dir)
write.xlsx(filter(df_final_sensibilidad_stoch, ciclo == 15), "Sensi_seed.xlsx")

almacen <- list()
for (i in 1:15) {
  almacen[[i]] <- df_final_sensibilidad_stoch %>%
    filter(ciclo == i) %>%
    summarize(Media = round(mean(area_infectada), 1), P5 = round(quantile(area_infectada, 0.05)), P95 = round(quantile(area_infectada, 0.95)),
              Mediana = median(area_infectada), Desviacion_Estandar = sd(area_infectada), CV_porcentaje = (sd(area_infectada) / mean(area_infectada)) * 100,
              Minimo = min(area_infectada), Maximo = max(area_infectada), Error_Estandar = sd(area_infectada) / sqrt(n()))
}
write.xlsx(do.call(rbind, almacen), "Fan_chart_random.xlsx")
print("====== COMPREHENSIVE SENSITIVITY PIPELINE COMPLETE ======")

```