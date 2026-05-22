# Spatio-Temporal Dispersal Modeling of Cassava Brown Streak Disease (CBSD)

This repository contains the R-based implementation of a spatial, stochastic cellular automaton designed to simulate the epidemiological spread of Cassava Brown Streak Disease (CBSD) across Southeast Asia. Developed under the framework of the Alliance of Bioversity International and CIAT, the model integrates spatial allocation parameters from the Spatial Production Allocation Model (SPAM 2020) with a coupled dual-component dispersal mechanism: Long-Distance Dispersal (LDD) driven by trade networks and simulated via an Inverse Power Law kernel, and Short-Distance Dispersal (SDD) driven by the biological flight vectors of the whitefly (*Bemisia tabaci*).

## Repository Architecture

The analytical framework is structured into two independent, standalone R scripts:

### 1. Core Spread and Policy Simulation Framework (`Spread_Model.R`)
This script serves as the primary simulation engine, executing data processing, parameter optimization, baseline calibration, and ex-post policy evaluation.
* **Spatial Processing & Downscaling:** Ingests global SPAM 2020 cassava harvested area layers, clips them to the specific study zone, filters out non-production noise using a $\ge 100$ hectares threshold, and executes spatial disaggregation to downscale the original grid into a high-resolution $3 \times 3\text{ km}$ sub-cell framework.
* **Geodesic Distance Computing:** Utilizes out-of-core memory management (`bigmemory`) to process a pairwise matrix of $40,248 \times 40,248$ pixel centroids, computing exact geodesic distances via the Haversine great-circle formula.
* **Dispersal Kernel Calibration:** Fits an Inverse Power Law (IPL) Dispersal Kernel ($f(d) = 1/d^k$) to map commercial or trade-mediated vector pathways across regional distances.
* **Parametric Grid Search Optimization:** Fits non-linear logistic and sigmoidal trajectories (Dose-Response framework) against historical data. Runs a multi-step sequence grid search to optimize the local whitefly transmission probability parameter (`prob_threshold`), minimizing the Root Mean Square Error ($RMSE$).
* **Policy Experimentation ("Border Control"):** Restricts transboundary pathways by implementing a user-defined $75\%$ quarantine efficiency barrier between specified geopolitical zones.
* **Impact Metrics Derivation:** Intersects multi-temporal binary infection mask surfaces with active agricultural grids to output Microsoft Excel spreadsheets (`.xlsx`) detailing regional biophysical impacts in metric tons and harvested hectares per province.

### 2. Comprehensive Sensitivity Analysis Suite (`Sensitivity_Analysis.R`)
This script evaluates structural parameter robustness and assesses model stability against stochastic noise through three parallelized (`future.apply`) Monte Carlo verification layers:
* **Part 1 (Spatial Seed Location Vulnerability):** Randomizes the initial patient zero origin node across regional boundaries ($100$ independent simulations per country) to calculate loss dispersion parameters (Mean, SD, Median, CV, Quantiles P5/P95, and Standard Error) for regional risk profiling.
* **Part 2 (Parameter Perturbation Analysis):** Holds the outbreak origin static while introducing a continuous uniform noise distribution ($\pm 10\%$ perturbation bounds) over the calibrated whitefly transmission threshold to evaluate coefficient stability via non-parametric correlation models (Spearman & Kendall).
* **Part 3 (Pure Stochastic Uncertainty):** Freezes both breakout locations and structural parameters under identical conditions across $100$ runs to isolate variance generated exclusively by pseudo-random seed numbers, compiling timeline vectors required to build confidence profile Fan Charts.
* **Geostatistical Interpolation (IDW):** Integrates with `gstat` to compute Inverse Distance Weighting interpolations, translating discrete simulation matrix endpoints into smoothed, continuous regional risk surfaces masked to active agricultural cells.
* **Asymptotic Convergence Verification:** Includes diagnostic logs validating model stability by contrasting asymptotic variance behaviors across sample cohorts ($100\text{ vs. }500$ operational simulations for the highest-variance territory).

---

## Repository Metadata & Indexing Information

For proper indexation within open data repositories (e.g., Dataverse, Zenodo), the following metadata specifications should be registered:

* **Time Period Covered:** `2015 – 2022`  
  *(Captures the temporal horizon from the initial observed field outbreak in Ratanakiri, Cambodia in 2015, up to the terminal validated historical tracking curve boundary in 2022).*
* **Time Period of Collection:** `2020 – 2026`  
  *(Reflects the SPAM baseline vintage framework coupled with active computational simulation engineering completed in 2026).*

---

## Core Simulation Constants

The following structural thresholds are established within the codebase:
* **Total Active Matrix Pixels (`pix`):** 40,248 production nodes.
* **Temporal Horizon Optimization (`num_steps`):** 15 Macro cycles (Each cycle represents a nested 10-month block).
* **Internal Timestep Resolution (`monthly_steps`):** 3 internal monthly iterations per macro cycle.
* **Biological Flight Boundary (`whithelfy`):** $3\text{ km}$ spatial buffer limit.
* **Optimal Calibration Transmission Weight (`prob_threshold`):** 0.60 (Calibrated via RMSE grid minimization).

---

## Portability & Execution Guidelines

Both scripts have been refactored to remove explicit absolute system dependencies, utilizing platform-independent relative structures via `file.path()`. 

To execute the computational pipeline:
1. Clone this repository onto your local workstation.
2. Open the scripts and modify the global tracking string pointer to map your local workspace root folder:
   ```R
   base_dir <- "path/to/your/local/root/workspace"

Contact:
Maria del Mar Esponda
m.esponda@cgiar.org



