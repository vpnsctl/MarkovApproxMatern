
process_pred_files <- function(file_paths) {
  all_dfs <- list()

  for (file in file_paths) {
      params <- extract_params_from_filename(file)
      N <- params$N
      n_obs <- params$n_obs
      range_val <- params$range

      if (grepl(".RDS$", file)) {
          pred_data <- readRDS(file)
      } else if (grepl(".csv$", file)) {
          pred_data <- read_csv(file)
      }
    
      pred_data <- tibble(data)
      colnames(pred_data) <- c(1:6, "nu")
      pred_data[["range"]] <-  range_val
      pred_data[["n_obs"]] <- n_obs
      pred_data[["N"]] <- N
      all_dfs[[file]] <- pred_data
  }
  
  all_dfs <- all_dfs

  combined_df <- process_dist(all_dfs)
  combined_df <- unique(combined_df)
  
  return(combined_df)
}
