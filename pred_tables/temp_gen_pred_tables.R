nu_values <- seq(0.01, 2.49, by = 0.01)

all_err_fem <- NULL

for (nu in nu_values) {
  file_path_with_two_decimals <- sprintf("pred_tables//5000_5000///range_2/fem/res_%.2f_5000_5000_range_2_fem.RDS", nu)
  file_path_with_one_decimal <- sprintf("pred_tables//5000_5000///range_2/fem/res_%.1f_5000_5000_range_2_fem.RDS", nu)
  file_path_whole_number <- sprintf("pred_tables//5000_5000///range_2/fem/res_%d_5000_5000_range_2_fem.RDS", as.integer(nu))

  if (file.exists(file_path_with_two_decimals)) {
    file_path <- file_path_with_two_decimals
  } else if (file.exists(file_path_with_one_decimal)) {
    file_path <- file_path_with_one_decimal
  } else if (file.exists(file_path_whole_number)) {
    file_path <- file_path_whole_number
  } else {
    message(sprintf("File not found for nu = %.2f", nu))
    next  
  }
  
  r <- readRDS(file_path)
  
  err_fem <- r$err.fem
  err_fem <- cbind(err_fem, nu = nu)
  
  all_err_fem <- rbind(all_err_fem, err_fem)
}

print(all_err_fem)
