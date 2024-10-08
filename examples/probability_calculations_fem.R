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

err_fem <- matrix(0, length(m_vec), length(n))

for (j in 1:n.rep) {
    obs_loc <- sort(runif(n.obs, 0, domain_upper_limit))
    while (min(diff(obs_loc)) < 1e-3) {
        obs_loc <- sort(runif(n.obs, 0, domain_upper_limit))
    }
    
    z <- matrix(rnorm(n.obs), ncol = n.obs)
    L <- chol(get_cov_mat(obs_loc, NULL, "true", nu, sqrt(8*nu)/range, sigma, NULL, NULL))
    y <- t(L) %*%t(z) + sigma.e * rnorm(n.obs)
    
    for (i in 1:length(n)) {
        partial_file <- paste0(partial_results_folder, "/partial_result_fem_rep", j, "_n", n[i], "_range", range, "_nu", nu, ".RDS")
        if (file.exists(partial_file)) {
                cat("Skipping rep =", j, "n =", n[i], "- results already exist\n")
                err_tmp <- readRDS(partial_file)
                err_fem[,i] <- err_tmp$err_fem
                next
        }
        if(n[i]>0){ 
            loc <- c(obs_loc, locs[[i]])

            obs.ind <- 1:length(loc)
            tmp <- sort(loc, index.return = TRUE)
            obs.ind[tmp$ix] <- 1:length(loc)
            loc <- tmp$x
            obs.ind <- obs.ind[1:n.obs]                

            if(any(diff(loc)<1e-3)){
                ind_remove <- which(diff(loc)<1e-3)
                ind_remove <- c(ind_remove, ind_remove+1)
                ind_remove <- setdiff(ind_remove, obs.ind)
                loc <- loc[-ind_remove]
                tolerance <- 1e-6
                obs.ind <- sapply(obs_loc, function(x) which(abs(loc - x) < tolerance)[1])     
            } 



        } else {
            loc <- obs_loc
            obs.ind <- 1:n.obs
        }

        if(is.null(calibrated_m[[as.character(n[i])]])){
            if(n[i] > 999){
                print("Calibrating...")
                previous_calibration <- auto_calibration_fem_rat(n=n[i]+n.obs, n_obs=n.obs, nu=nu, range=range, sigma=sigma, sigma_e=sigma.e, 
                    samples=samples_calibration, m_rat=m_vec, previous_calibration = previous_calibration, max_it_per_m = max_it_per_m, print=FALSE)
                calibrated_m[[as.character(n[i])]] <-  previous_calibration                   
                cat("Calibration:", previous_calibration, "\n")
            } else{
                calibrated_m[[as.character(n[i])]] <- calibrated_m[[as.character(1000)]]
            }
        } else{
            cat("rep =", j, "n =", n[i], "calibration loaded\n")
            cat("Calibration: ", calibrated_m[[as.character(n[i])]], "\n")
        }
        cat('compute truth\n')
        true_cov <- get_cov_mat(loc = loc, m = NULL, method = "true", nu = nu, kappa = kappa, sigma = sigma, 
                                samples = NULL, L=NULL)
        post_true <- posterior_constructor_cov(true_cov, y, sigma.e, obs.ind, 1:(n[i]), type = "all")
        post_mean_true <- as.vector(post_true$post_mean)
        post_cov_true <- as.matrix(post_true$post_cov)
        
        Q.true <- solve(post_cov_true)
        conf <- simconf(alpha = 1-coverage, mu = post_mean_true, Q = Q.true, vars  = diag(post_cov_true), n.iter=1e5)
        lb_prob <- conf$a
        ub_prob <- conf$b
        if(use.excursions && n[i] > 999) {
            prob_true <- gaussint(a=lb_prob,b=ub_prob,mu=post_mean_true,Q = Q.true, n.iter = 1e5)$P        
        } else {
            prob_true <- pmvnorm(lower=lb_prob,upper=ub_prob,mean=post_mean_true,sigma = post_cov_true)
        }
        
        for(k in 1:length(m_vec)){
            m <- m_vec[k]
            cat("rep = ", j,"n =", n[i], 'rational, m = ', m, '\n')
            mn <- calibrated_m[[as.character(n[i])]][k]
            post_fem_mat <- fem_post_calculations(y = y, loc = loc, idx_obs = obs.ind, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma.e, m = m, mesh_fem = mn)
            if(use.excursions && n[i] > 999) {
                prob_fem <- gaussint(a=lb_prob,b=ub_prob,mu=post_fem_mat$mu.fem,Q = solve(post_fem_mat$Sigma.fem), n.iter = 1e5)$P
            } else {
                prob_fem <- pmvnorm(lower=lb_prob,upper=ub_prob,mean=mu.fem,sigma = post_fem_mat$Sigma.fem)
            }
            err_fem[k,i] <- err_fem[k,i] + prob_fem-prob_true
            
            cat("prob fem: ", prob_fem, '\n')
            cat("error fem: ", err_fem[k,i]/j, '\n')
        }
        
        # Save partial results
        partial_res_list <- list(
            err_fem = err_fem[, i],
            n = n[i],
            rep = j,
            m = m_vec
        )
        saveRDS(partial_res_list, file = partial_file)
    }
}

# Average the errors across replicates
err_fem <- err_fem / n.rep

# Save final results
res_list <- list(err_fem = err_fem, nu = nu, m_vec = m_vec, n_vec = n)
saveRDS(res_list, file = paste0("prob_tables/fem_range", range, "_nu", nu, ".RDS"))
