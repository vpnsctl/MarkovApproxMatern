rm(list=ls())
source("aux_functions/calibration_functions.R")

# Define the vector of `n` values
n            <- c(0:10, seq(20,50,by=10), seq(75,250,by=25), seq(300,500, by = 50))
n            <- rev(n)        # Reverse the vector
n.obs        <- 1001

range <- 1
sigma <- 1
sigma.e <- 0.1
nu <- 1
alpha <- nu + 1/2
samples_calibration <- 50
m_vec <- 1:6
max_it_per_m <- 20
calibrated_m <- list()
previous_calibration <- NULL

# Set the directory to save calibration
calibration_file <- "calibration_results/calibrated_m_list_regular_mesh.RDS"
if (!dir.exists("calibration_results")) {
    dir.create("calibration_results")
}

# Perform calibration and save it
for (i in 1:length(n)) {
    if(n[i] > 249){    
        print(paste("Calibrating for n =", n[i]))
        previous_calibration <- auto_calibration_nngp_rat(
            n = n[i] + n.obs, n_obs = n.obs, nu = nu, range = range,
            sigma = sigma, sigma_e = sigma.e, samples = samples_calibration, 
            m_rat = m_vec, previous_calibration = previous_calibration, max_it_per_m = max_it_per_m, print = FALSE
        )
        calibrated_m[[as.character(n[i])]] <- previous_calibration
        cat("Calibration for n =", n[i], ":", previous_calibration, "\n")
            } else{
                calibrated_m[[as.character(n[i])]] <- calibrated_m[[as.character(250)]]
            }        
}

calibrated_m <- lapply(calibrated_m,adjust_m)

# Save the calibration results
saveRDS(calibrated_m, file = calibration_file)
cat("Calibrations saved to:", calibration_file, "\n")
