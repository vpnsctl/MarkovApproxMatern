rm(list=ls())
source("aux_functions/calibration_functions.R")
source("aux_functions/aux_functions_cov.R")
source("aux_functions/probability_computations.R")
library(rSPDE)
library(excursions)
library(mvtnorm)

# Define paths for saving calibration and results
calibration_file <- "calibration_results/calibrated_m_list_fem.RDS"
partial_results_folder <- "prob_tables/partial_results"

# Create directories for saving results
if (!dir.exists("calibration_results")) {
    dir.create("calibration_results")
}
if (!dir.exists("prob_tables")) {
    dir.create("prob_tables")
}
if (!dir.exists(partial_results_folder)) {
    dir.create(partial_results_folder)
}

# Set parameters
range <- 0.5
sigma <- 1
sigma.e <- 0.1
n <- c(0, 25, 50, 100, 150, 250, 300, 400, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000)
n <- n[length(n):1]
n.obs <- 250
n.rep <- 10
samples_calibration <- 50
max_it_per_m <- 20
nu <- 1
alpha <- nu + 1/2
m_vec <- 1:6
domain_upper_limit <- 10
use.excursions <- TRUE
coverage <- 0.9

calibrated_m <- list()
previous_calibration <- NULL
kappa <- sqrt(8*nu)/range

locs <- lapply(n, function(x) if(x > 0) seq(from = 1, to = domain_upper_limit, length.out = x) else NULL)

if (file.exists(calibration_file)) {
    calibrated_m <- readRDS(calibration_file)
    cat("Loaded calibrated m from", calibration_file, "\n")
} 

    for (i in 1:length(n)) {
                    if(n[i] > 999){
                print("Calibrating...")
                previous_calibration <- auto_calibration_fem_rat(n=n[i]+n.obs, n_obs=n.obs, nu=nu, range=range, sigma=sigma, sigma_e=sigma.e, 
                    samples=samples_calibration, m_rat=m_vec, max_it_per_m = max_it_per_m, print=FALSE)
                calibrated_m[[as.character(n[i])]] <-  previous_calibration                   
                cat("Calibration:", previous_calibration, "\n")
            } else{
                calibrated_m[[as.character(n[i])]] <- calibrated_m[[as.character(1000)]]
            }
    }

saveRDS(calibrated_m, calibration_file)


