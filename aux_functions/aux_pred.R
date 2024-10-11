library(dplyr)
library(tidyr)
library(stringr)
library(readr)

extract_pred_params_from_filename <- function(filename) {
    file_base <- tools::file_path_sans_ext(basename(filename))
    
    match <- regmatches(file_base, regexec("res_(\\d+)_(\\d+)_range([0-9.]+)_([a-zA-Z]+)", file_base))
    
    if (length(match[[1]]) == 5) {
        N <- as.integer(match[[1]][2])
        n_obs <- as.integer(match[[1]][3])
        range <- as.numeric(match[[1]][4])
        method <- match[[1]][5]
        method <- switch(
          method,
          "fourier" = "Fourier",
          "statespace" = "State-Space",
          "rational" = "Rational",
          "nngp" = "nnGP",
          "pca" = "PCA",
          "fem" = "FEM",
          "taper" = "Taper",
          stop("Unknown method: ", method)
        )     
        
        return(list(N = N, n_obs = n_obs, range = range, method = method))
    } else {
        warning(paste("Filename does not match expected pattern:", filename))
        return(NULL)
    }
}

process_pred_files <- function(file_paths) {
  all_dfs <- list()

  for (file in file_paths) {
      cat("Processing file: ", file, "\n")
      params <- extract_pred_params_from_filename(file)
      N <- params$N
      n_obs <- params$n_obs
      range_val <- params$range
      method <- params$method

      if (grepl(".RDS$", file)) {
          pred_data <- readRDS(file)
          pred_data <- pred_data[[grep("^err\\.", names(pred_data))]]
          colnames(pred_data) <- c(1:6, "nu")
      } else if (grepl(".csv$", file)) {
          pred_data <- read.csv(file)
          colnames(pred_data) <- c("nu", 1:6)
      }
      cat("File was read successfully!\n")
      pred_data <- tibble(pred_data)
      pred_data[["range"]] <-  range_val
      pred_data[["n_obs"]] <- n_obs
      pred_data[["N"]] <- N
      pred_data[["Method"]] <- method
      all_dfs[[file]] <- pred_data
  }
  
  all_dfs <- do.call(bind_rows,all_dfs)
  all_dfs <- all_dfs |>
    pivot_longer(
      cols = `1`:`6`,  
      names_to = "m",  
      values_to = "pred_error"  
    )
    
  return(all_dfs)
}
