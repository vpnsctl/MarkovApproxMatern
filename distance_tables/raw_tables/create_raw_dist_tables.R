source("aux_functions/aux_functions_cov.R")
source("aux_functions/predict_methods.R")
source("aux_functions/distance_computations.R")
source("aux_functions/aux_dist.R")

N <- c(500, 1000)

nu_vec=seq(0.01,2.99, by = 0.01)

idx <- (nu_vec + 0.5)%%1 > 1e-10
nu_vec <- nu_vec[idx]
range = 0.5
sigma = 1
m <- 1:6
samples = 10

nu <- 0.6
m <- 1


dist_rat <- compute_distances_rational(N=N, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma)
saveRDS(dist_rat, "distance_tables/dist_rational_range1.RDS")

dist_nngp <- compute_distances_nngp(N=N, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma)
saveRDS(dist_nngp, "distance_tables/dist_nngp_range1v2.RDS")

dist_pca <- compute_distances_pca(N=N, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma)
saveRDS(dist_pca, "distance_tables/dist_pca_range1.RDS")

dist_fourier <- compute_distances_fourier(N=N, m.vec=m, nu.vec=nu_vec, range=range, samples=samples)
saveRDS(dist_fourier, "distance_tables/dist_fourier_range1.RDS")

dist_ss <- compute_distances_statespace(N=N, m.vec=m, nu.vec=nu_vec, range=range)
saveRDS(dist_ss, "distance_tables/dist_ss_range1.RDS")







dist_nngp <- compute_distances_nngp(N=N, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, type = "simulation")
saveRDS(dist_nngp, "distance_tables/dist_nngp_samp_range05.RDS")

dist_pca <- compute_distances_pca(N=N, m.vec=m, nu.vec=nu_vec, range=range, sigma=sigma, type = "simulation")
saveRDS(dist_pca, "distance_tables/dist_pca_samp_range1.RDS")

dist_fourier <- compute_distances_fourier(N=N, m.vec=m, nu.vec=nu_vec, range=range, samples=samples, type = "simulation")
saveRDS(dist_fourier, "distance_tables/dist_fourier_samp_range1.RDS")

dist_ss <- compute_distances_statespace(N=N, m.vec=m, nu.vec=nu_vec, range=range, type = "simulation")
saveRDS(dist_ss, "distance_tables/dist_ss_samp_range1.RDS")