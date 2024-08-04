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

# N= 5000, nobs = 5000 - with 200 replicates
# alpha01: 1 1 1 1 1 2
# alpha12: 1 2 13 22 28 32
# alpha12: 1 2 17 20 27 30
# alpha12: 1 2 15 23 28 33
# alpha23: 17 31 39 45 50 55
# alpha23: 20 28 37 47 52 56
# alpha23: 17 31 39 46 51 56


# N=10000, nobs = 5000 - with 200 replicates
# alpha01: 1 1 1 1 1 1
# alpha12: 1  1  1  2 15 21
# alpha12: 1  1  1  2 16 22
# alpha12: 1  1  1  6 16 22
# alpha12: 1  2  3  6 16 22
# alpha12: 1  1  1 12 18 23
# alpha23: 15 23 33 43 51 56 (descarta primeiro)
# alpha23: 15 24 32 39 48 53
# alpha23: 15 24 35 43 51 58

# N=10000, nobs = 10000 - with 200 replicates
# alpha01: 1 1 1 1 1 1
# alpha12: 1 2 6 18 25 30
# alpha12: 1 2 3 9  20 31
# alpha12: 1 2 8 18 26 31
# alpha12: 1 1 10 18 26 31
# alpha23: 15 29 40 46 52 57
# alpha23: 15 31 40 46 52 56
# alpha23: 16 30 39 47 54 60

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



## Tabulating rational timings - prediction


rm(list=ls())
source("aux_functions/calibration_functions.R")

N <- 5000
n_obs <- 5000
nu <- 0.4
range <- 1
sigma <- 1
sigma_e <- 0.1
m_rat <- 1:6
samples = 200

times.rat01 <- timing_rat(N, n_obs, nu = 0.4, range, sigma, sigma_e, m_rat, samples, type = "prediction")
saveRDS(times.rat01, "time_tables//times_tmp//times_rat_5000_01_v2.RDS")
times.rat12 <- timing_rat(N, n_obs, nu = 1.4, range, sigma, sigma_e, m_rat, samples, type = "prediction")
saveRDS(times.rat12, "time_tables//times_tmp//times_rat_5000_12_v2.RDS")
times.rat23 <- timing_rat(N, n_obs, nu = 2.4, range, sigma, sigma_e, m_rat, samples, type = "prediction")
saveRDS(times.rat23, "time_tables//times_tmp//times_rat_5000_v2_23.RDS")





rm(list=ls())
source("aux_functions/calibration_functions.R")

N <- 10000
n_obs <- 5000
nu <- 0.4
range <- 1
sigma <- 1
sigma_e <- 0.1
m_rat <- 1:6
samples = 200

times.rat01 <- timing_rat(N, n_obs, nu = 0.4, range, sigma, sigma_e, m_rat, samples, type = "prediction")
saveRDS(times.rat01, "time_tables//times_calibration//times_rat_10000_5000_01.RDS")
times.rat12 <- timing_rat(N, n_obs, nu = 1.4, range, sigma, sigma_e, m_rat, samples, type = "prediction")
saveRDS(times.rat12, "time_tables//times_calibration//times_rat_10000_5000_12_v2_200.RDS")
times.rat23 <- timing_rat(N, n_obs, nu = 2.4, range, sigma, sigma_e, m_rat, samples, type = "prediction")
saveRDS(times.rat23, "time_tables//times_calibration//times_rat_10000_5000_23_v2_200.RDS")




rm(list=ls())
source("aux_functions/calibration_functions.R")

N <- 10000
n_obs <- 10000
nu <- 0.4
range <- 1
sigma <- 1
sigma_e <- 0.1
m_rat <- 1:6
samples = 200

times.rat01 <- timing_rat(N, n_obs, nu = 0.4, range, sigma, sigma_e, m_rat, samples, type = "prediction")
saveRDS(times.rat01, "time_tables//times_calibration//times_rat_10000_10000_01.RDS")
times.rat12 <- timing_rat(N, n_obs, nu = 1.4, range, sigma, sigma_e, m_rat, samples, type = "prediction")
saveRDS(times.rat12, "time_tables//times_calibration//times_rat_10000_10000_12_v2.RDS")
times.rat23 <- timing_rat(N, n_obs, nu = 2.4, range, sigma, sigma_e, m_rat, samples, type = "prediction")
saveRDS(times.rat23, "time_tables//times_calibration//times_rat_10000_10000_23_v2.RDS")


#### NNGP

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
samples = 200

m.nngp_pred_01 <- 1:15
times.nngp01_pred <- timing_nngp_m_vec(N = N, n_obs = n_obs, m.nngp = m.nngp_pred_01, nu = 1.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = samples, type = "prediction")
saveRDS(times.nngp01_pred, "time_tables//times_calibration//times_nngp_5000_01_v2.RDS")

m.nngp_pred_12 <- 1:40
times.nngp12_pred <- timing_nngp_m_vec(N = N, n_obs = n_obs, m.nngp = m.nngp_pred_12, nu = 1.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = samples, type = "prediction")
saveRDS(times.nngp12_pred, "time_tables//times_calibration//times_nngp_5000_12.RDS")

m.nngp_pred_23 <- 15:65
times.nngp23_pred <- timing_nngp_m_vec(N = N, n_obs = n_obs, m.nngp = m.nngp_pred_23, nu = 2.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = samples, type = "prediction")
saveRDS(times.nngp23_pred, "time_tables//times_calibration//times_nngp_5000_23.RDS")



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
samples = 200

m.nngp_pred_01 <- 1:15
times.nngp01_pred <- timing_nngp_m_vec(N = N, n_obs = n_obs, m.nngp = m.nngp_pred_01, nu = 1.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = samples, type = "prediction")
saveRDS(times.nngp01_pred, "time_tables//times_calibration//times_nngp_10000_5000_01.RDS")

m.nngp_pred_12 <- 1:40
times.nngp12_pred <- timing_nngp_m_vec(N = N, n_obs = n_obs, m.nngp = m.nngp_pred_12, nu = 1.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = samples, type = "prediction")
saveRDS(times.nngp12_pred, "time_tables//times_calibration//times_nngp_10000_5000_12_v3.RDS")

m.nngp_pred_23 <- 15:65
times.nngp23_pred <- timing_nngp_m_vec(N = N, n_obs = n_obs, m.nngp = m.nngp_pred_23, nu = 2.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = samples, type = "prediction")
saveRDS(times.nngp23_pred, "time_tables//times_calibration//times_nngp_10000_5000_23_v3.RDS")


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
samples = 200

m.nngp_pred_01 <- 1:15
times.nngp01_pred <- timing_nngp_m_vec(N = N, n_obs = n_obs, m.nngp = m.nngp_pred_01, nu = 1.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = samples, type = "prediction")
saveRDS(times.nngp01_pred, "time_tables//times_calibration//times_nngp_10000_10000_01.RDS")

m.nngp_pred_12 <- 1:40
times.nngp12_pred <- timing_nngp_m_vec(N = N, n_obs = n_obs, m.nngp = m.nngp_pred_12, nu = 1.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = samples, type = "prediction")
saveRDS(times.nngp12_pred, "time_tables//times_calibration//times_nngp_10000_10000_12_v3.RDS")

m.nngp_pred_23 <- 15:65
times.nngp23_pred <- timing_nngp_m_vec(N = N, n_obs = n_obs, m.nngp = m.nngp_pred_23, nu = 2.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = samples, type = "prediction")
saveRDS(times.nngp23_pred, "time_tables//times_calibration//times_nngp_10000_10000_23_v2.RDS")

