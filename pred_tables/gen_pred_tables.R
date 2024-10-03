source("aux_functions/aux_pred.R")

files <- list.files("pred_tables/full_tables", full.names = TRUE)

pred_error <- process_pred_files(files)
saveRDS(pred_error, "pred_tables/pred_error.RDS")
