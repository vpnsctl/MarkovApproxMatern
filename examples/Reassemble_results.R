assemble_results <- function(input_dir, method) {
  error_element <- switch(method,
                          "rational" = "err.rat",
                          "nngp" = "err.nn",
                          "pca" = "err.pca",
                          "fourier" = "err.fourier",
                          "statespace" = "err.ss",
                          stop("Invalid method specified"))

  rds_files <- list.files(input_dir, pattern = "\\.RDS$", full.names = TRUE)

  res_pred <- list(nu = numeric(), err = NULL)

  for (file in rds_files) {
    res <- readRDS(file)
    
    nu_value <- res$nu
    
    res_pred$nu <- c(res_pred$nu, nu_value)

    err <- res[[error_element]]
    err$nu <- nu_value
    
    if (is.null(res_pred$err)) {
      res_pred$err <- err
    } else {
      res_pred$err <- rbind(res_pred$err, err)
    }
  }
  return(res_pred)
}

# Example
input_dir <- "pred_tables//10000_10000//range_2//rational"
method <- "rational"
res_pred <- assemble_results(input_dir, method)
