source("aux_functions/aux_functions_cov.R")
source("aux_functions/predict_methods.R")
source("aux_functions/distance_computations.R")
source("aux_functions/aux_dist.R")


## General setup:
## Ranges considered are: 0.5, 1 and 2
## Domain length = N/100
## m for the different methods were calibrated by the total prediction time.
## sigma = 1
## No calibration were carried out against the parsimonious method.

nu_vec <- seq(from = 0.01, to = 2.49, by = 0.01)
idx <- (nu_vec + 0.5)%%1 > 1e-10
nu_vec_rat <- nu_vec[idx]
sigma = 1
m <- 1:6
m_rat <- 0:6


## Config 1:
## N = n_obs = 5000

N <- n_obs <- 5000

m_nngp_fun <- function(m, alpha){
            if(alpha<1) {
                mn <- m - 1
                if(mn < 1){
                    mn <- 1
                }
            } else if (alpha < 2) {
                m_vec <- c(1, 2, 13, 21, 27, 31) 
                mn <- m_vec[m]
            } else {
                m_vec <- c(15, 30, 37, 45, 51, 54)
                mn <- m_vec[m]
            }
            return(mn)
} 

m_pca_fun <- function(m, alpha){
            if(alpha<1) {
                m_vec <- c(268, 308, 355, 406, 433, 478)
                mn <- m_vec[m]
            } else if (alpha < 2) {
                m_vec <- c(380, 473, 561, 651, 708, 776)
                mn <- m_vec[m]
            } else {
                m_vec <- c(611, 810, 945, 1082, 1205, 1325)
                mn <- m_vec[m]
            }
            return(mn)
}

## Same as PCA, but with 10 samples
m_fourier_fun <- function(m, alpha){
            if(alpha<1) {
                m_vec <- c(268, 308, 355, 406, 433, 478)
                mn <- m_vec[m]
            } else if (alpha < 2) {
                m_vec <- c(380, 473, 561, 651, 708, 776)
                mn <- m_vec[m]
            } else {
                m_vec <- c(611, 810, 945, 1082, 1205, 1325)
                mn <- m_vec[m]
            }
            return(mn)
}

m_statespace_fun <- function(m, alpha){
    return(max(c(1,m - floor(alpha))))
}


### First range

range <- 0.5

dist_rat <- compute_distances_rational(N=N, n_obs = n_obs, m.vec=m_rat, nu.vec=nu_vec_rat, range=range, sigma=sigma)
saveRDS(dist_rat, "distance_tables/raw_tables/dist_rational_5000_range_05_calibrated.RDS")

dist_nngp <- compute_distances_nngp(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_nngp_fun = m_nngp_fun)
saveRDS(dist_nngp, "distance_tables/raw_tables/dist_nngp_5000_range_05_calibrated.RDS")

dist_pca <- compute_distances_pca(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_pca_fun = m_pca_fun)
saveRDS(dist_pca, "distance_tables/raw_tables/dist_pca_5000_range_05_calibrated.RDS")

dist_ss <- compute_distances_statespace(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_statespace_fun = m_statespace_fun)
saveRDS(dist_ss, "distance_tables/raw_tables/dist_ss_5000_range_05_calibrated.RDS")


### second range

range <- 1

dist_rat <- compute_distances_rational(N=N, n_obs = n_obs, m.vec=m_rat, nu.vec=nu_vec_rat, range=range, sigma=sigma)
saveRDS(dist_rat, "distance_tables/raw_tables/dist_rational_5000_range_1_calibrated.RDS")

dist_nngp <- compute_distances_nngp(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_nngp_fun = m_nngp_fun)
saveRDS(dist_nngp, "distance_tables/raw_tables/dist_nngp_5000_range_1_calibrated.RDS")

dist_pca <- compute_distances_pca(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_pca_fun = m_pca_fun)
saveRDS(dist_pca, "distance_tables/raw_tables/dist_pca_5000_range_1_calibrated.RDS")

dist_ss <- compute_distances_statespace(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_statespace_fun = m_statespace_fun)
saveRDS(dist_ss, "distance_tables/raw_tables/dist_ss_5000_range_1_calibrated.RDS")

### Third range

range <- 2

dist_rat <- compute_distances_rational(N=N, n_obs = n_obs, m.vec=m_rat, nu.vec=nu_vec_rat, range=range, sigma=sigma)
saveRDS(dist_rat, "distance_tables/raw_tables/dist_rational_5000_range_2_calibrated.RDS")

dist_nngp <- compute_distances_nngp(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_nngp_fun = m_nngp_fun)
saveRDS(dist_nngp, "distance_tables/raw_tables/dist_nngp_5000_range_2_calibrated.RDS")

dist_pca <- compute_distances_pca(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_pca_fun = m_pca_fun)
saveRDS(dist_pca, "distance_tables/raw_tables/dist_pca_5000_range_2_calibrated.RDS")

dist_ss <- compute_distances_statespace(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_statespace_fun = m_statespace_fun)
saveRDS(dist_ss, "distance_tables/raw_tables/dist_ss_5000_range_2_calibrated.RDS")


###############

source("aux_functions/aux_functions_cov.R")
source("aux_functions/predict_methods.R")
source("aux_functions/distance_computations.R")
source("aux_functions/aux_dist.R")


## General setup:
## Ranges considered are: 0.5, 1 and 2
## Domain length = N/100
## m for the different methods were calibrated by the total prediction time.
## sigma = 1
## No calibration were carried out against the parsimonious method.

nu_vec <- seq(from = 0.01, to = 2.49, by = 0.01)
idx <- (nu_vec + 0.5)%%1 > 1e-10
nu_vec_rat <- nu_vec[idx]
sigma = 1
m <- 1:6
m_rat <- 0:6



## Config 2:

# Calibration for N = 10000 and n_obs = 5000

N <- 10000
n_obs <- 10000


m_nngp_fun <- function(m, alpha){
            if(alpha<1) {
                mn <- m - 1
                if(mn < 1){
                    mn <- 1
                }
            } else if (alpha < 2) {
                m_vec <- c(1, 2, 3, 4, 14, 20)
                mn <- m_vec[m]
            } else {
                m_vec <- c(1, 22, 31, 39, 47, 54)
                mn <- m_vec[m]
            }
            return(mn)
} 

m_pca_fun <- function(m, alpha){
            if(alpha<1) {
                m_vec <- c(381, 448, 533, 588, 654, 727)
                mn <- m_vec[m]
            } else if (alpha < 2) {
                m_vec <- c(532, 704, 844, 953, 1065, 1162)
                mn <- m_vec[m]
            } else {
                m_vec <- c(904, 1202, 1431, 1622, 1808, 1965)
                mn <- m_vec[m]
            }
            return(mn)
}

## Same as PCA, but with 100 samples
m_fourier_fun <- function(m, alpha){
            if(alpha<1) {
                m_vec <- c(381, 448, 533, 588, 654, 727)
                mn <- m_vec[m]
            } else if (alpha < 2) {
                m_vec <- c(532, 704, 844, 953, 1065, 1162)
                mn <- m_vec[m]
            } else {
                m_vec <- c(904, 1202, 1431, 1622, 1808, 1965)
                mn <- m_vec[m]
            }
            return(mn)
}

m_statespace_fun <- function(m, alpha){
    return(max(c(1,m - floor(alpha))))
}


### First range

range <- 0.5

dist_rat <- compute_distances_rational(N=N, n_obs = n_obs, m.vec=m_rat, nu.vec=nu_vec_rat, range=range, sigma=sigma)
saveRDS(dist_rat, "distance_tables/raw_tables/dist_rational_10000_5000_range_05_calibrated.RDS")

dist_nngp <- compute_distances_nngp(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_nngp_fun = m_nngp_fun)
saveRDS(dist_nngp, "distance_tables/raw_tables/dist_nngp_10000_5000_range_05_calibrated.RDS")

dist_pca <- compute_distances_pca(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_pca_fun = m_pca_fun)
saveRDS(dist_pca, "distance_tables/raw_tables/dist_pca_10000_5000_range_05_calibrated.RDS")

dist_ss <- compute_distances_statespace(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_statespace_fun = m_statespace_fun)
saveRDS(dist_ss, "distance_tables/raw_tables/dist_ss_10000_5000_range_05_calibrated.RDS")


### second range

range <- 1

dist_rat <- compute_distances_rational(N=N, n_obs = n_obs, m.vec=m_rat, nu.vec=nu_vec_rat, range=range, sigma=sigma)
saveRDS(dist_rat, "distance_tables/raw_tables/dist_rational_10000_5000_range_1_calibrated.RDS")

dist_nngp <- compute_distances_nngp(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_nngp_fun = m_nngp_fun)
saveRDS(dist_nngp, "distance_tables/raw_tables/dist_nngp_10000_5000_range_1_calibrated.RDS")

dist_pca <- compute_distances_pca(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_pca_fun = m_pca_fun)
saveRDS(dist_pca, "distance_tables/raw_tables/dist_pca_10000_5000_range_1_calibrated.RDS")

dist_ss <- compute_distances_statespace(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_statespace_fun = m_statespace_fun)
saveRDS(dist_ss, "distance_tables/raw_tables/dist_ss_10000_5000_range_1_calibrated.RDS")

### Third range

range <- 2

dist_rat <- compute_distances_rational(N=N, n_obs = n_obs, m.vec=m_rat, nu.vec=nu_vec_rat, range=range, sigma=sigma)
saveRDS(dist_rat, "distance_tables/raw_tables/dist_rational_10000_5000_range_2_calibrated.RDS")

dist_nngp <- compute_distances_nngp(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_nngp_fun = m_nngp_fun)
saveRDS(dist_nngp, "distance_tables/raw_tables/dist_nngp_10000_5000_range_2_calibrated.RDS")

dist_pca <- compute_distances_pca(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_pca_fun = m_pca_fun)
saveRDS(dist_pca, "distance_tables/raw_tables/dist_pca_10000_5000_range_2_calibrated.RDS")

dist_ss <- compute_distances_statespace(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_statespace_fun = m_statespace_fun)
saveRDS(dist_ss, "distance_tables/raw_tables/dist_ss_10000_5000_range_2_calibrated.RDS")

#####################




## Config 3:
## N = 10000, n_obs = 10000

source("aux_functions/aux_functions_cov.R")
source("aux_functions/predict_methods.R")
source("aux_functions/distance_computations.R")
source("aux_functions/aux_dist.R")


## General setup:
## Ranges considered are: 0.5, 1 and 2
## Domain length = N/100
## m for the different methods were calibrated by the total prediction time.
## sigma = 1
## No calibration were carried out against the parsimonious method.

nu_vec <- seq(from = 0.01, to = 2.49, by = 0.01)
idx <- (nu_vec + 0.5)%%1 > 1e-10
nu_vec_rat <- nu_vec[idx]
sigma = 1
m <- 1:6
m_rat <- 0:6



N <- 10000
n_obs <- 10000

m_nngp_fun <- function(m, alpha){
            if(alpha<1) {
                mn <- m - 1
                if(mn < 1){
                    mn <- 1
                }
            } else if (alpha < 2) {
                m_vec <- c(1, 2, 7, 18, 24, 29)
                mn <- m_vec[m]
            } else {
                m_vec <- c(14, 28, 37, 44, 51, 57)
                mn <- m_vec[m]
            }
            return(mn)
} 

m_pca_fun <- function(m, alpha){
            if(alpha<1) {
                m_vec <- c(271, 311, 367, 412, 452, 493)
                mn <- m_vec[m]
            } else if (alpha < 2) {
                m_vec <- c(372, 493, 591, 672, 751, 821)
                mn <- m_vec[m]
            } else {
                m_vec <- c(640, 848, 1016, 1168, 1299, 1420)
                mn <- m_vec[m]
            }
            return(mn)
}

## Same as PCA, but with 100 samples
m_fourier_fun <- function(m, alpha){
            if(alpha<1) {
                m_vec <- c(271, 311, 367, 412, 452, 493)
                mn <- m_vec[m]
            } else if (alpha < 2) {
                m_vec <- c(372, 493, 591, 672, 751, 821)
                mn <- m_vec[m]
            } else {
                m_vec <- c(640, 848, 1016, 1168, 1299, 1420)
                mn <- m_vec[m]
            }
            return(mn)
}

m_statespace_fun <- function(m, alpha){
    return(max(c(1,m - floor(alpha))))
}


### First range

range <- 0.5

dist_rat <- compute_distances_rational(N=N, n_obs = n_obs, m.vec=m_rat, nu.vec=nu_vec_rat, range=range, sigma=sigma)
saveRDS(dist_rat, "distance_tables/raw_tables/dist_rational_10000_10000_range_05_calibrated.RDS")

dist_nngp <- compute_distances_nngp(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_nngp_fun = m_nngp_fun)
saveRDS(dist_nngp, "distance_tables/raw_tables/dist_nngp_10000_10000_range_05_calibrated.RDS")

dist_pca <- compute_distances_pca(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_pca_fun = m_pca_fun)
saveRDS(dist_pca, "distance_tables/raw_tables/dist_pca_10000_10000_range_05_calibrated.RDS")

dist_ss <- compute_distances_statespace(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_statespace_fun = m_statespace_fun)
saveRDS(dist_ss, "distance_tables/raw_tables/dist_ss_10000_10000_range_05_calibrated.RDS")


### second range

range <- 1

dist_rat <- compute_distances_rational(N=N, n_obs = n_obs, m.vec=m_rat, nu.vec=nu_vec_rat, range=range, sigma=sigma)
saveRDS(dist_rat, "distance_tables/raw_tables/dist_rational_10000_10000_range_1_calibrated.RDS")

dist_nngp <- compute_distances_nngp(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_nngp_fun = m_nngp_fun)
saveRDS(dist_nngp, "distance_tables/raw_tables/dist_nngp_10000_10000_range_1_calibrated.RDS")

dist_pca <- compute_distances_pca(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_pca_fun = m_pca_fun)
saveRDS(dist_pca, "distance_tables/raw_tables/dist_pca_10000_10000_range_1_calibrated.RDS")

dist_ss <- compute_distances_statespace(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_statespace_fun = m_statespace_fun)
saveRDS(dist_ss, "distance_tables/raw_tables/dist_ss_10000_10000_range_1_calibrated.RDS")

### Third range

range <- 2

dist_rat <- compute_distances_rational(N=N, n_obs = n_obs, m.vec=m_rat, nu.vec=nu_vec_rat, range=range, sigma=sigma)
saveRDS(dist_rat, "distance_tables/raw_tables/dist_rational_10000_10000_range_2_calibrated.RDS")

dist_nngp <- compute_distances_nngp(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_nngp_fun = m_nngp_fun)
saveRDS(dist_nngp, "distance_tables/raw_tables/dist_nngp_10000_10000_range_2_calibrated.RDS")

dist_pca <- compute_distances_pca(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_pca_fun = m_pca_fun)
saveRDS(dist_pca, "distance_tables/raw_tables/dist_pca_10000_10000_range_2_calibrated.RDS")

dist_ss <- compute_distances_statespace(N=N, n_obs = n_obs, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, m_statespace_fun = m_statespace_fun)
saveRDS(dist_ss, "distance_tables/raw_tables/dist_ss_10000_10000_range_2_calibrated.RDS")


###################

