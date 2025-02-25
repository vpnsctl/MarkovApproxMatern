rm(list=ls())
source("aux_functions/calibration_functions.R")
source("aux_functions/aux_functions_cov.R")
source("aux_functions/probability_computations.R")
library(rSPDE)
library(excursions)
library(mvtnorm)

# Define paths for saving calibration and results
calibration_file <- "calibration_results/calibrated_m_list_fem_regular_mesh.RDS"
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
n            <- c(0:10, seq(20,50,by=10), seq(75,250,by=25), seq(300,500, by = 50))
n            <- rev(n)        # Reverse the vector
n.obs        <- 1001
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

full_mesh <- seq(0, domain_upper_limit, length.out = n.obs+max(n))

locs <- lapply(n, function(i) {
  full_mesh[1:(n.obs + i)]
})

if (file.exists(calibration_file)) {
    calibrated_m <- readRDS(calibration_file)
    cat("Loaded calibrated m from", calibration_file, "\n")
} 

    for (i in 1:length(n)) {
        if(n[i] > 249){
        print(paste("Calibrating for n =", n[i]))
                previous_calibration <- auto_calibration_fem_rat(n=n[i]+n.obs, n_obs=n.obs, nu=nu, range=range, sigma=sigma, sigma_e=sigma.e, 
                    samples=samples_calibration, m_rat=m_vec, max_it_per_m = max_it_per_m, print=TRUE)
                calibrated_m[[as.character(n[i])]] <-  previous_calibration                   
                cat("Calibration:", previous_calibration, "\n")
            } else{
                calibrated_m[[as.character(n[i])]] <- calibrated_m[[as.character(250)]]
            }
    }

calibrated_m <- lapply(calibrated_m,adjust_m)

saveRDS(calibrated_m, calibration_file)


