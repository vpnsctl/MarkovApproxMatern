source("aux_functions/aux_pred.R")

files <- list.files("pred_tables/full_tables", full.names = TRUE)

df_dist <- process_pred_files(files)