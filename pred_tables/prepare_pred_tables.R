# Create empty vectors to store the data
nu_values <- seq(0.01, 2.49, by = 0.01)
pred_error <- data.frame()

# Loop through all nu values
for(nu in nu_values) {
  # Create the filename with the correct format
  filename <- sprintf("pred_tables//5000_5000//range_2//fem//res_%.2f_5000_5000_range_2_sigmae_0.10_fem.RDS", nu)
  
  # Read the RDS file
  data <- readRDS(filename)
  
  # Extract the values and create temporary data frame
  for(m in 1:6) {
    temp_df <- data.frame(
      nu = nu,
      method = "fem",
      mu_error = data$mu_errors[1,m],
      sigma_error = data$sigma_errors[1,m],
      m = m,
      sigma_e = 0.1
    )
    
    # Append to the main data frame
    pred_error <- rbind(pred_error, temp_df)
  }
}

# Loop through all nu values
for(nu in nu_values) {
  # Create the filename with the correct format
  filename <- sprintf("pred_tables//5000_5000//range_2//fem//res_%.2f_5000_5000_range_2_sigmae_0.32_fem.RDS", nu)
  
  # Read the RDS file
  data <- readRDS(filename)
  
  # Extract the values and create temporary data frame
  for(m in 1:6) {
    temp_df <- data.frame(
      nu = nu,
      method = "fem",
      mu_error = data$mu_errors[1,m],
      sigma_error = data$sigma_errors[1,m],
      m = m,
      sigma_e = sqrt(0.1)
    )
    
    # Append to the main data frame
    pred_error <- rbind(pred_error, temp_df)
  }
}

# Loop through all nu values
for(nu in nu_values) {
  # Create the filename with the correct format
  filename <- sprintf("pred_tables//5000_5000//range_2//nngp//res_%.2f_5000_5000_range_2_sigmae_0.10_nngp.RDS", nu)
  
  # Read the RDS file
  data <- readRDS(filename)
  
  # Extract the values and create temporary data frame
  for(m in 1:6) {
    temp_df <- data.frame(
      nu = nu,
      method = "nngp",
      mu_error = data$mu_errors[1,m],
      sigma_error = data$sigma_errors[1,m],
      m = m,
      sigma_e = 0.1
    )
    
    # Append to the main data frame
    pred_error <- rbind(pred_error, temp_df)
  }
}

# Loop through all nu values
for(nu in nu_values) {
  # Create the filename with the correct format
  filename <- sprintf("pred_tables//5000_5000//range_2//nngp//res_%.2f_5000_5000_range_2_sigmae_0.32_nngp.RDS", nu)
  
  # Read the RDS file
  data <- readRDS(filename)
  
  # Extract the values and create temporary data frame
  for(m in 1:6) {
    temp_df <- data.frame(
      nu = nu,
      method = "nngp",
      mu_error = data$mu_errors[1,m],
      sigma_error = data$sigma_errors[1,m],
      m = m,
      sigma_e = sqrt(0.1)
    )
    
    # Append to the main data frame
    pred_error <- rbind(pred_error, temp_df)
  }
}

# Loop through all nu values
for(nu in nu_values) {

    if(!(nu%in%c(0.5,1.5))){
        # Create the filename with the correct format
        filename <- sprintf("pred_tables//5000_5000//range_2//rational//res_%.2f_5000_5000_range_2_sigmae_0.10_rational.RDS", nu)
        
        # Read the RDS file
        data <- readRDS(filename)
        
        # Extract the values and create temporary data frame
        for(m in 1:6) {
          temp_df <- data.frame(
            nu = nu,
            method = "rational",
            mu_error = data$mu_errors[1,m],
            sigma_error = data$sigma_errors[1,m],
            m = m,
            sigma_e = 0.1
          )
          
          # Append to the main data frame
          pred_error <- rbind(pred_error, temp_df)
        }
    } else{
        for(m in 1:6) {
          temp_df <- data.frame(
            nu = nu,
            method = "rational",
            mu_error = 1e-16,
            sigma_error = 1e-16,
            m = m,
            sigma_e = 0.1
          )
          
          # Append to the main data frame
          pred_error <- rbind(pred_error, temp_df)
        }
    }

}

# Loop through all nu values
for(nu in nu_values) {

    if(!(nu%in%c(0.5,1.5))){
        # Create the filename with the correct format
        filename <- sprintf("pred_tables//5000_5000//range_2//rational//res_%.2f_5000_5000_range_2_sigmae_0.32_rational.RDS", nu)
        
        # Read the RDS file
        data <- readRDS(filename)
        
        # Extract the values and create temporary data frame
        for(m in 1:6) {
          temp_df <- data.frame(
            nu = nu,
            method = "rational",
            mu_error = data$mu_errors[1,m],
            sigma_error = data$sigma_errors[1,m],
            m = m,
            sigma_e = sqrt(0.1)
          )
          
          # Append to the main data frame
          pred_error <- rbind(pred_error, temp_df)
        }
    } else{
        for(m in 1:6) {
          temp_df <- data.frame(
            nu = nu,
            method = "rational",
            mu_error = 1e-16,
            sigma_error = 1e-16,
            m = m,
            sigma_e = sqrt(0.1)
          )
          
          # Append to the main data frame
          pred_error <- rbind(pred_error, temp_df)
        }
    }

}



# Loop through all nu values
for(nu in nu_values) {
  # Create the filename with the correct format
  filename <- sprintf("pred_tables//5000_5000//range_2//taper//res_%.2f_5000_5000_range_2_sigmae_0.10_taper.RDS", nu)
  
  # Read the RDS file
  data <- readRDS(filename)
  
  # Extract the values and create temporary data frame
  for(m in 1:6) {
    temp_df <- data.frame(
      nu = nu,
      method = "taper",
      mu_error = data$mu_errors[1,m],
      sigma_error = data$sigma_errors[1,m],
      m = m,
      sigma_e = 0.1
    )
    
    # Append to the main data frame
    pred_error <- rbind(pred_error, temp_df)
  }
}

# Loop through all nu values
for(nu in nu_values) {
  # Create the filename with the correct format
  filename <- sprintf("pred_tables//5000_5000//range_2//taper//res_%.2f_5000_5000_range_2_sigmae_0.32_taper.RDS", nu)
  
  # Read the RDS file
  data <- readRDS(filename)
  
  # Extract the values and create temporary data frame
  for(m in 1:6) {
    temp_df <- data.frame(
      nu = nu,
      method = "taper",
      mu_error = data$mu_errors[1,m],
      sigma_error = data$sigma_errors[1,m],
      m = m,
      sigma_e = sqrt(0.1)
    )
    
    # Append to the main data frame
    pred_error <- rbind(pred_error, temp_df)
  }
}

all_data <- list()

methods <- c("fourier", "statespace", "pca")
sigma_values <- c("0.10", "0.32")

for (method in methods) {
  for (sigma in sigma_values) {
    filename <- paste0("pred_tables//", method, "_pred_5000_5000_range2_sigmae_", sigma, ".csv")
    
    df <- read.csv(filename)
    
    df$sigma_e <- ifelse(sigma == "0.10", 0.1, sqrt(0.1))
    df$method <- method
    
    all_data[[paste(method, sigma)]] <- df
  }
}

combined_data <- do.call(rbind, all_data)
rownames(combined_data) <- NULL

final_data <- rbind(pred_error, combined_data)

saveRDS(final_data, "pred_tables//pred_error_range2_5000_5000_complete.RDS")
