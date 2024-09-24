# NNGP Calibrations

# N= 5000, nobs = 5000 - with 200 replicates
# alpha12: 1 2 13 21 27 31
# alpha23: 15 30 37 45 51 54
1
# N=10000, nobs = 5000 - with 200 replicates
# alpha12:  1  1  1  2 14 20 ->  1  2  3  4 14 20
# alpha23:  1 22 31 39 47 54

# N=10000, nobs = 10000 - with 200 replicates
# alpha12: 1  1  7 18 24 29 ->  1  2  7 18 24 29
# alpha23: 14 28 37 44 51 57

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



timings_rat_5000_12 <- readRDS("time_tables//times_calibration//times_rat_5000_12_cumsum_FALSE.RDS")
timings_rat_5000_23 <- readRDS("time_tables//times_calibration//times_rat_5000_23_cumsum_FALSE.RDS")
timings_rat_10000_5000_12 <- readRDS("time_tables//times_calibration//times_rat_10000_5000_12_cumsum_FALSE.RDS")
timings_rat_10000_5000_23 <- readRDS("time_tables//times_calibration//times_rat_10000_5000_23_cumsum_FALSE.RDS")
timings_rat_10000_10000_12 <- readRDS("time_tables//times_calibration//times_rat_10000_10000_12_cumsum_FALSE.RDS")
timings_rat_10000_10000_23 <- readRDS("time_tables//times_calibration//times_rat_10000_10000_23_cumsum_FALSE.RDS")

timings_nngp_5000_12 <- readRDS("time_tables//times_calibration//times_nngp_5000_12.RDS")
timings_nngp_5000_23 <- readRDS("time_tables//times_calibration//times_nngp_5000_23.RDS")
timings_nngp_10000_5000_12 <- readRDS("time_tables//times_calibration//times_nngp_10000_5000_12.RDS")
timings_nngp_10000_5000_23 <- readRDS("time_tables//times_calibration//times_nngp_10000_5000_23.RDS")
timings_nngp_10000_10000_12 <- readRDS("time_tables//times_calibration//times_nngp_10000_10000_12.RDS")
timings_nngp_10000_10000_23 <- readRDS("time_tables//times_calibration//times_nngp_10000_10000_23.RDS")

m.rat <- 1:6

m.nngp12 <- 1:40
m.nngp23 <- 5:65
m.nngp23_10000_5000 <- 1:65

calibration(m.rat = m.rat, times.rat = timings_rat_5000_12, m.other = m.nngp12, times.other = timings_nngp_5000_12)
calibration(m.rat = m.rat, times.rat = timings_rat_5000_23, m.other = m.nngp23, times.other = timings_nngp_5000_23)
calibration(m.rat = m.rat, times.rat = timings_rat_10000_5000_12, m.other = m.nngp12, times.other = timings_nngp_10000_5000_12)
calibration(m.rat = m.rat, times.rat = timings_rat_10000_5000_23, m.other = m.nngp23_10000_5000, times.other = timings_nngp_10000_5000_23)
calibration(m.rat = m.rat, times.rat = timings_rat_10000_10000_12, m.other = m.nngp12, times.other = timings_nngp_10000_10000_12)
calibration(m.rat = m.rat, times.rat = timings_rat_10000_10000_23, m.other = m.nngp23, times.other = timings_nngp_10000_10000_23)

N <- 250
n_obs <- 250
nu <- 0.4
range <- 1
sigma <- 1
sigma_e <- 0.1
m_rat <- 1:6
samples = 200

times.rat01 <- timing_rat(N, n_obs, nu = 0.4, range, sigma, sigma_e, m_rat, samples, type = "prediction")
times.rat12 <- timing_rat(N, n_obs, nu = 1.4, range, sigma, sigma_e, m_rat, samples, type = "prediction")
times.rat23 <- timing_rat(N, n_obs, nu = 2.4, range, sigma, sigma_e, m_rat, samples, type = "prediction")
saveRDS(times.rat01, "times_rat_250_01.RDS")
saveRDS(times.rat12, "times_rat_250_12.RDS")
saveRDS(times.rat23, "times_rat_250_23.RDS")

m.nngp_pred_01 <- 1:15
times.nngp01_pred <- timing_nngp_m_vec(N = N, n_obs = n_obs, m.nngp = m.nngp_pred_01, nu = 1.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = samples, type = "prediction")
saveRDS(times.nngp01_pred, "times_nngp_250_01.RDS")

m.nngp_pred_12 <- 1:20
times.nngp12_pred <- timing_nngp_m_vec(N = N, n_obs = n_obs, m.nngp = m.nngp_pred_12, nu = 1.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = samples, type = "prediction")
saveRDS(times.nngp12_pred, "times_nngp_250_12.RDS")

m.nngp_pred_23 <- 1:30
times.nngp23_pred <- timing_nngp_m_vec(N = N, n_obs = n_obs, m.nngp = m.nngp_pred_23, nu = 2.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = samples, type = "prediction")
saveRDS(times.nngp23_pred, "times_nngp_250_23.RDS")


m.rat <- 1:6

calibration(m.rat = m.rat, times.rat = times.rat12, m.other = m.nngp_pred_12, times.other = times.nngp12_pred)
calibration(m.rat = m.rat, times.rat = times.rat23, m.other = m.nngp_pred_23, times.other = times.nngp23_pred)




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
saveRDS(times.rat01, "time_tables//times_tmp//times_rat_5000_01_cumsum_FALSE.RDS")
times.rat12 <- timing_rat(N, n_obs, nu = 1.4, range, sigma, sigma_e, m_rat, samples, type = "prediction")
saveRDS(times.rat12, "time_tables//times_tmp//times_rat_5000_12_cumsum_FALSE.RDS")
times.rat23 <- timing_rat(N, n_obs, nu = 2.4, range, sigma, sigma_e, m_rat, samples, type = "prediction")
saveRDS(times.rat23, "time_tables//times_tmp//times_rat_5000_23_cumsum_FALSE.RDS")





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
saveRDS(times.rat01, "time_tables//times_calibration//times_rat_10000_5000_01_cumsum_FALSE.RDS")
times.rat12 <- timing_rat(N, n_obs, nu = 1.4, range, sigma, sigma_e, m_rat, samples, type = "prediction")
saveRDS(times.rat12, "time_tables//times_calibration//times_rat_10000_5000_12_cumsum_FALSE.RDS")
times.rat23 <- timing_rat(N, n_obs, nu = 2.4, range, sigma, sigma_e, m_rat, samples, type = "prediction")
saveRDS(times.rat23, "time_tables//times_calibration//times_rat_10000_5000_23_cumsum_FALSE.RDS")




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
saveRDS(times.rat01, "time_tables//times_calibration//times_rat_10000_10000_01_cumsum_FALSE.RDS")
times.rat12 <- timing_rat(N, n_obs, nu = 1.4, range, sigma, sigma_e, m_rat, samples, type = "prediction")
saveRDS(times.rat12, "time_tables//times_calibration//times_rat_10000_10000_12_cumsum_FALSE.RDS")
times.rat23 <- timing_rat(N, n_obs, nu = 2.4, range, sigma, sigma_e, m_rat, samples, type = "prediction")
saveRDS(times.rat23, "time_tables//times_calibration//times_rat_10000_10000_23_cumsum_FALSE.RDS")


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

m.nngp_pred_23 <- 5:14
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

m.nngp_pred_23 <- 1:4
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

m.nngp_pred_23 <- 5:14
times.nngp23_pred <- timing_nngp_m_vec(N = N, n_obs = n_obs, m.nngp = m.nngp_pred_23, nu = 2.4, range = range, sigma = sigma, sigma_e = sigma_e, samples = samples, type = "prediction")
saveRDS(times.nngp23_pred, "time_tables//times_calibration//times_nngp_10000_10000_23_v2.RDS")

