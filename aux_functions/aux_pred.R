library(dplyr)
library(tidyr)
library(stringr)

extract_pred_params_from_filename <- function(filename) {
    file_base <- tools::file_path_sans_ext(basename(filename))
    
    match <- regmatches(file_base, regexec("res_(\\d+)_(\\d+)_range([0-9.]+)_([a-zA-Z]+)", file_base))
    
    if (length(match[[1]]) == 5) {
        N <- as.integer(match[[1]][2])
        n_obs <- as.integer(match[[1]][3])
        range <- as.numeric(match[[1]][4])
        method <- match[[1]][5]
        
        return(list(N = N, n_obs = n_obs, range = range, method = method))
    } else {
        warning(paste("Filename does not match expected pattern:", filename))
        return(NULL)
    }
}

process_pred_files <- function(file_paths) {
  all_dfs <- list()

  for (file in file_paths) {
      params <- extract_pred_params_from_filename(file)
    print(params)
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
