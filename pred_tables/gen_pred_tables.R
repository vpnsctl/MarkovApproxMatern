source("aux_functions/aux_dist.R")




files <- list.files("pred_tables/full_tables", full.names = TRUE)

df_dist <- process_all_files(files)