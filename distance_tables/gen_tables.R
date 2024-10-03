source("aux_functions/aux_dist.R")


files <- list.files("distance_tables/raw_tables", full.names = TRUE)
files <- files[!grepl(".pkl$", files)]
df_dist <- process_all_files(files)

saveRDS(df_dist, "distance_tables/full_dists.RDS")

plot_dist(df_dist, N=5000, n_obs=5000, range_val = 0.5, distance="L2", methods = c("nnGP"))
