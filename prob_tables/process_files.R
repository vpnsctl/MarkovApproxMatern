# Function to extract parameters from filename
get_params <- function(filename) {
  range <- as.numeric(gsub(".*range(\\d+(\\.\\d+)?).*", "\\1", filename))
  sigma_e <- as.numeric(gsub(".*sigmaE(\\d+(\\.\\d+)?).*", "\\1", filename))
  is_regular <- grepl("regular_mesh", filename)
  return(list(range = range, sigma_e = sigma_e, is_regular = is_regular))
}

# Function to process each file
process_file <- function(file_path) {
  data <- readRDS(file_path)
  params <- get_params(file_path)
  
  n_len <- length(data$n_vec)
  m_len <- length(data$m_vec)
  
  # Process rational errors
  df_rat <- data.frame(
    m = rep(data$m_vec, each = n_len),
    N = rep(data$n_vec, times = m_len),
    error = as.vector(t(data$err_rat)),
    Method = "Rational",
    range = params$range,
    sigma_e = params$sigma_e
  )
  
  # Process NNGP errors
  df_nngp <- data.frame(
    m = rep(data$m_vec, each = n_len),
    N = rep(data$n_vec, times = m_len),
    error = as.vector(t(data$err_nngp)),
    Method = "NNGP",
    range = params$range,
    sigma_e = params$sigma_e
  )
  
  # Process FEM errors
  df_fem <- data.frame(
    m = rep(data$m_vec, each = n_len),
    N = rep(data$n_vec, times = m_len),
    error = as.vector(t(data$err_fem)),
    Method = "FEM",
    range = params$range,
    sigma_e = params$sigma_e
  )
  
  combined_df <- rbind(df_rat, df_nngp, df_fem)
  return(combined_df)
}

# List all RDS files
files <- list.files(path = "prob_tables", pattern = "rat_nngp_fem.*\\.RDS$", full.names = TRUE)

# Separate regular and non-regular files
regular_files <- files[grepl("regular_mesh", files)]
nonregular_files <- files[!grepl("regular_mesh", files)]

# Process regular mesh files
regular_df <- do.call(rbind, lapply(regular_files, process_file))

# Process non-regular mesh files
nonregular_df <- do.call(rbind, lapply(nonregular_files, process_file))

# Read prob_errors.RDS
prob_errors <- readRDS("prob_tables/prob_errors.RDS")

# Add sigma_e column to prob_errors (set to 0.1) and rename prob_error to error
prob_errors_modified <- transform(prob_errors, 
                                sigma_e = 0.1,
                                error = prob_error)[, c("m", "N", "error", "Method", "range", "sigma_e")]

# Combine prob_errors with non-regular data
nonregular_df <- rbind(nonregular_df, prob_errors_modified)