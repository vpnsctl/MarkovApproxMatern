source("aux_functions/aux_dist.R")


files <- list.files("distance_tables/raw_tables", full.names = TRUE)
files <- files[!grepl(".pkl$", files)]
df_dist <- process_all_files(files)

plot_dist(df_dist, n_loc=10000, n_obs=10000, range = 2, distance="L2", method = "Rational")
