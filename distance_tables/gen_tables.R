source("aux_functions/aux_functions_cov.R")
source("aux_functions/distance_computations.R")
source("aux_functions/aux_dist.R")


# Range 0.2

dist_ss_02 <- readRDS("distance_tables/raw_tables/dist_ss_range02.RDS")
dist_pca_02 <- readRDS("distance_tables/raw_tables/dist_pca_range02.RDS")
dist_fourier_02 <- readRDS("distance_tables/raw_tables/dist_fourier_range02.RDS")
dist_rat_02 <- readRDS("distance_tables/raw_tables/dist_rational_range02.RDS")
dist_nngp_02 <- readRDS("distance_tables/raw_tables/dist_nngp_range02.RDS")

range_02 <- process_dist(dist_rat_02, dist_pca_02, dist_ss_02, dist_fourier_02)
saveRDS(range_02, "distance_tables/range_02.RDS")


# Range 0.2 sampling

dist_samp_ss_02 <- readRDS("distance_tables/raw_tables/dist_ss_samp_range02.RDS")
dist_samp_pca_02 <- readRDS("distance_tables/raw_tables/dist_pca_samp_range02.RDS")
dist_samp_fourier_02 <- readRDS("distance_tables/raw_tables/dist_fourier_samp_range02.RDS")
dist_samp_nngp_02 <- readRDS("distance_tables/raw_tables/dist_nngp_samp_range02.RDS")

range_02_samp <- process_dist(dist_rat_02, dist_samp_pca_02, dist_samp_ss_02, dist_samp_fourier_02, dist_samp_nngp_02)
saveRDS(range_02_samp, "distance_tables/range_02_samp.RDS")


# Range 0.5

dist_ss_05 <- readRDS("distance_tables/raw_tables/dist_ss_range05.RDS")
dist_pca_05 <- readRDS("distance_tables/raw_tables/dist_pca_range05.RDS")
dist_fourier_05 <- readRDS("distance_tables/raw_tables/dist_fourier_range05.RDS")
dist_rat_05 <- readRDS("distance_tables/raw_tables/dist_rational_range05.RDS")
dist_nngp_05 <- readRDS("distance_tables/raw_tables/dist_nngp_range05.RDS")

range_05 <- process_dist(dist_rat_05, dist_pca_05, dist_ss_05, dist_fourier_05)
saveRDS(range_05, "distance_tables/range_05.RDS")


# Range 0.5 sampling

dist_samp_ss_05 <- readRDS("distance_tables/raw_tables/dist_ss_samp_range05.RDS")
dist_samp_pca_05 <- readRDS("distance_tables/raw_tables/dist_pca_samp_range05.RDS")
dist_samp_fourier_05 <- readRDS("distance_tables/raw_tables/dist_fourier_samp_range05.RDS")
dist_samp_nngp_05 <- readRDS("distance_tables/raw_tables/dist_nngp_samp_range05.RDS")

range_05_samp <- process_dist(dist_rat_05, dist_samp_pca_05, dist_samp_ss_05, dist_samp_fourier_05, dist_samp_nngp_05)
saveRDS(range_05_samp, "distance_tables/range_05_samp.RDS")

# Range 1

dist_ss_1 <- readRDS("distance_tables/raw_tables/dist_ss_range1.RDS")
dist_pca_1 <- readRDS("distance_tables/raw_tables/dist_pca_range1.RDS")
dist_fourier_1 <- readRDS("distance_tables/raw_tables/dist_fourier_range1.RDS")
dist_rat_1 <- readRDS("distance_tables/raw_tables/dist_rational_range1.RDS")
dist_nngp_1 <- readRDS("distance_tables/raw_tables/dist_nngp_range1.RDS")

range_1 <- process_dist(dist_rat_1, dist_pca_1, dist_ss_1, dist_fourier_1)
saveRDS(range_1, "distance_tables/range_1.RDS")


# Range 1 sampling

dist_samp_ss_1 <- readRDS("distance_tables/raw_tables/dist_ss_samp_range1.RDS")
dist_samp_pca_1 <- readRDS("distance_tables/raw_tables/dist_pca_samp_range1.RDS")
dist_samp_fourier_1 <- readRDS("distance_tables/raw_tables/dist_fourier_samp_range1.RDS")
dist_samp_nngp_1 <- readRDS("distance_tables/raw_tables/dist_nngp_samp_range1.RDS")

range_1_samp <- process_dist(dist_rat_1, dist_samp_pca_1, dist_samp_ss_1, dist_samp_fourier_1, dist_samp_nngp_1)
saveRDS(range_1_samp, "distance_tables/range_1_samp.RDS")



attr(dist_ss_1, "m.vec") <- attr(dist_rat_05, "m.vec")
attr(dist_ss_1, "nu.vec") <- attr(dist_rat_05, "nu.vec")
attr(dist_pca_1, "m.vec") <- attr(dist_rat_05, "m.vec")
attr(dist_pca_1, "nu.vec") <- attr(dist_rat_05, "nu.vec")
attr(dist_fourier_1, "m.vec") <- attr(dist_rat_05, "m.vec")
attr(dist_fourier_1, "nu.vec") <- attr(dist_rat_05, "nu.vec")
attr(dist_samp_nngp_1, "m.vec") <- attr(dist_rat_1, "m.vec")
attr(dist_samp_nngp_1, "nu.vec") <- attr(dist_rat_1, "nu.vec")

saveRDS(dist_ss_1, "distance_tables/raw_tables/dist_ss_range1.RDS")
saveRDS(dist_pca_1, "distance_tables/raw_tables/dist_pca_range1.RDS")
saveRDS(dist_fourier_1, "distance_tables/raw_tables/dist_fourier_range1.RDS")
