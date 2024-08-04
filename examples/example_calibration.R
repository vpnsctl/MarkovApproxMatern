## Calibrations:

# NNGP (calibrations for m from 2 to 6):
# N = n_obs = 5000 
# (prediction) 
# 1 < alpha < 2: 2, 17, 20, 27, 30
# 2 < alpha < 3: 28, 37, 47, 52, 56
# (sampling) - including build the precision
# 1 < alpha < 2: 2,  9,  19, 22, 28
# 2 < alpha < 3: 30, 38, ,46, 54, 59

# N = 10000, n_obs = 5000
# (prediction) 
# (sampling) 
# 1 < alpha < 2: 2 2 12 22 29
# 2 < alpha < 3: 33 44 49 52 60

# N = 10000, n_obs = 10000
# (prediction) 
# 1 < alpha < 2: 2  3  9  20 31
# 2 < alpha < 3: 32 42 48 53 58
# (sampling) 
# 1 < alpha < 2: 3 3 3 16 23
# 2 < alpha < 3: 26 38 49 55 62


# NNGP Calibrations

# N= 5000, nobs = 5000
# alpha01: 268 308 355 406 433 478
# alpha12: 380 473 561 651 708 776
# alpha23: 611  810  945 1082 1205 1325

# N=10000, nobs = 5000
# alpha01: 381 448 533 588 654 727
# alpha12: 532  704  844  953 1065 1162
# alpha23: 904 1202 1431 1622 1808 1965

# N=10000, nobs = 10000
# alpha01: 271 311 367 412 452 493
# alpha12: 372 493 591 672 751 821
# alpha23: 640  848 1016 1168 1299 1420

# PCA Calibrations

# N= 5000, nobs = 5000
# alpha01: 268 308 355 406 433 478
# alpha12: 380 473 561 651 708 776
# alpha23: 611  810  945 1082 1205 1325

# N=10000, nobs = 5000
# alpha01: 381 448 533 588 654 727
# alpha12: 532  704  844  953 1065 1162
# alpha23: 904 1202 1431 1622 1808 1965

# N=10000, nobs = 10000
# alpha01: 271 311 367 412 452 493
# alpha12: 372 493 591 672 751 821
# alpha23: 640  848 1016 1168 1299 1420


rm(list=ls())
source("aux_functions/calibration_functions.R")

N <- 5000
n_obs <- 5000
m_min <- 2
m_max <- 1000
m_step <- 1
nu <- 0.4
range <- 1
sigma <- 1
sigma_e <- 0.1
m_rat <- 1:6
samples = 10

times.rat01 <- timing_rat(N, n_obs, nu = 0.4, range, sigma, sigma_e, m_rat, samples, type = "prediction")
times.rat12 <- timing_rat(N, n_obs, nu = 1.4, range, sigma, sigma_e, m_rat, samples, type = "prediction")
times.rat23 <- timing_rat(N, n_obs, nu = 2.4, range, sigma, sigma_e, m_rat, samples, type = "prediction")

times.rat <- c(times.rat01, times.rat12, times.rat23)

times.nngp01 <- timing_nngp(N = N, n_obs = n_obs, m_min = 2, m_max = 10, m_step = 1, nu = 0.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = 2, type = "estimation", est_nu = FALSE, only_optim = TRUE)
times.nngp12 <- timing_nngp(N = N, n_obs = n_obs, m_min = 2, m_max = 40, m_step = 1, nu = 1.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = 2, type = "estimation", est_nu = FALSE, only_optim = TRUE)
times.nngp23 <- timing_nngp(N = N, n_obs = n_obs, m_min = 20, m_max = 65, m_step = 1, nu = 2.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = 2, type = "estimation", est_nu = FALSE, only_optim = TRUE)

pca_calibration <- search_PCA(N, n_obs, m_min, m_max, times.rat, samples, type = "prediction")


pca_calibration01 <- pca_calibration$m.cal[1:6]
pca_calibration12 <- pca_calibration$m.cal[7:12]
pca_calibration23 <- pca_calibration$m.cal[13:18]


plot(pca_calibration$m.pca, pca_calibration$times.pca)
for(i in 1:length(times.rat)){
  lines(c(1,pca_calibration$m.cal[i]), times.rat[i]*c(1,1))
  lines(c(pca_calibration$m.cal[i],pca_calibration$m.cal[i]), c(0,times.rat[i]))
}


rm(list=ls())
source("aux_functions/calibration_functions.R")

N <- 5000
n_obs <- 5000
m_min <- 2
m_max <- 1000
m_step <- 1
nu <- 0.4
range <- 1
sigma <- 1
sigma_e <- 0.1
m_rat <- 1:6
samples = 10

m.nngp_pred_12 <- c(2:5, 8:12,17:21, 27:31)
times.nngp12_pred <- timing_nngp_m_vec(N = N, n_obs = n_obs, m.nngp = m.nngp_pred_12, nu = 1.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = 100, type = "prediction", est_nu = FALSE, only_optim = TRUE)

m.nngp_pred_23 <- c(20:29, 34:38, 43:56)
times.nngp23_pred <- timing_nngp_m_vec(N = N, n_obs = n_obs, m.nngp = m.nngp_pred_23, nu = 2.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = 100, type = "prediction", est_nu = FALSE, only_optim = TRUE)



rm(list=ls())
source("aux_functions/calibration_functions.R")

N <- 10000
n_obs <- 5000
m_min <- 2
m_max <- 1000
m_step <- 1
nu <- 0.4
range <- 1
sigma <- 1
sigma_e <- 0.1
m_rat <- 1:6
samples = 10

m.nngp_pred_12 <- 2:22
times.nngp12_pred <- timing_nngp_m_vec(N = N, n_obs = n_obs, m.nngp = m.nngp_pred_12, nu = 1.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = 100, type = "prediction", est_nu = FALSE, only_optim = TRUE)

m.nngp_pred_23 <- c(18:26, 33:47, 49:53)
times.nngp23_pred <- timing_nngp_m_vec(N = N, n_obs = n_obs, m.nngp = m.nngp_pred_23, nu = 2.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = 100, type = "prediction", est_nu = FALSE, only_optim = TRUE)


rm(list=ls())
source("aux_functions/calibration_functions.R")

N <- 10000
n_obs <- 10000
m_min <- 2
m_max <- 1000
m_step <- 1
nu <- 0.4
range <- 1
sigma <- 1
sigma_e <- 0.1
m_rat <- 1:6
samples = 10

m.nngp_pred_12 <- c(2:9, 16:20, 24:28, 31:35)
times.nngp12_pred <- timing_nngp_m_vec(N = N, n_obs = n_obs, m.nngp = m.nngp_pred_12, nu = 1.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = 100, type = "prediction", est_nu = FALSE, only_optim = TRUE)

m.nngp_pred_23 <- c(20:24, 32:36, 41:49, 50:59)
times.nngp23_pred <- timing_nngp_m_vec(N = N, n_obs = n_obs, m.nngp = m.nngp_pred_23, nu = 2.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = 100, type = "prediction", est_nu = FALSE, only_optim = TRUE)


rm(list=ls())
source("aux_functions/calibration_functions.R")

N <- 5000
n_obs <- 5000
m_min <- 2
m_max <- 1000
m_step <- 1
nu <- 0.4
range <- 1
sigma <- 1
sigma_e <- 0.1
m_rat <- 1:6
samples = 10

m.nngp_samp_12 <- c(2:7, 9:13, 19:23, 28:32)
times.nngp12_samp <- timing_nngp_m_vec(N = N, n_obs = n_obs, m.nngp = m.nngp_samp_12, nu = 1.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = 100, type = "sampling", est_nu = FALSE, only_optim = TRUE)

m.nngp_samp_23 <- c(18:22, 30:34, 38:42, 46:50, 54:62)
times.nngp23_samp <- timing_nngp_m_vec(N = N, n_obs = n_obs, m.nngp = m.nngp_samp_23, nu = 2.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = 100, type = "sampling", est_nu = FALSE, only_optim = TRUE)


rm(list=ls())
source("aux_functions/calibration_functions.R")

N <- 10000
n_obs <- 5000
m_min <- 2
m_max <- 1000
m_step <- 1
nu <- 0.4
range <- 1
sigma <- 1
sigma_e <- 0.1
m_rat <- 1:6
samples = 10

m.nngp_samp_12 <- c(2:6, 12:16, 22:26, 29:33)
times.nngp12_samp <- timing_nngp_m_vec(N = N, n_obs = n_obs, m.nngp = m.nngp_samp_12, nu = 1.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = 100, type = "sampling", est_nu = FALSE, only_optim = TRUE)

m.nngp_samp_23 <- c(18:22, 31:35, 40:44, 47:54, 58:62)
times.nngp23_samp <- timing_nngp_m_vec(N = N, n_obs = n_obs, m.nngp = m.nngp_samp_23, nu = 2.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = 100, type = "sampling", est_nu = FALSE, only_optim = TRUE)

rm(list=ls())
source("aux_functions/calibration_functions.R")

N <- 10000
n_obs <- 10000
m_min <- 2
m_max <- 1000
m_step <- 1
nu <- 0.4
range <- 1
sigma <- 1
sigma_e <- 0.1
m_rat <- 1:6
samples = 10

m.nngp_samp_12 <- c(3:12, 16:20, 23:27)
times.nngp12_samp <- timing_nngp_m_vec(N = N, n_obs = n_obs, m.nngp = m.nngp_samp_12, nu = 1.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = 100, type = "sampling", est_nu = FALSE, only_optim = TRUE)

m.nngp_samp_23 <- c(18:22, 26:30, 38:42, 45:49, 51:55, 58:62)
times.nngp23_samp <- timing_nngp_m_vec(N = N, n_obs = n_obs, m.nngp = m.nngp_samp_23, nu = 2.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = 100, type = "sampling", est_nu = FALSE, only_optim = TRUE)