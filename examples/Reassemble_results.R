input_dir <- "pred_tables//10000_10000//range_2//rational"

rds_files <- list.files(input_dir, pattern = "\\.RDS$", full.names = TRUE)

res_pred <- list(nu = numeric(), err.rat = NULL)

for (file in rds_files) {
  
  res <- readRDS(file)

  nu_value <- res$nu

  res_pred$nu <- c(res_pred$nu, nu_value)

  err.rat <- res$err.rat
  err.rat$nu <- nu_value
  
  if (is.null(res_pred$err.rat)) {
    res_pred$err.rat <- err.rat
  } else {
    res_pred$err.rat <- rbind(res_pred$err.rat, err.rat)
  }
}

print(res_pred)
