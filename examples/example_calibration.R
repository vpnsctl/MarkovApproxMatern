rm(list=ls())
source("aux_functions/calibration_functions.R")

N <- 1000
n_obs <- 1000
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
pca_calibration <- search_PCA(N, n_obs, m_min, m_max, times.rat, samples, type = "prediction")

pca_calibration01 <- pca_calibration$m.cal[1:6]
pca_calibration12 <- pca_calibration$m.cal[7:12]
pca_calibration23 <- pca_calibration$m.cal[13:18]


plot(pca_calibration$m.pca, pca_calibration$times.pca)
for(i in 1:length(times.rat)){
  lines(c(1,pca_calibration$m.cal[i]), times.rat[i]*c(1,1))
  lines(c(pca_calibration$m.cal[i],pca_calibration$m.cal[i]), c(0,times.rat[i]))
}

