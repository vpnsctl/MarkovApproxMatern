source("aux_functions/aux_dist.R")


files <- list.files("distance_tables/raw_tables", full.names = TRUE)
files <- files[!grepl(".pkl$", files)]
combined_results <- process_all_files(files)
