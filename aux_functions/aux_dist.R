library(dplyr)
library(tidyr)
library(stringr)
library(rhdf5)

# Function to extract N, n_obs, and range from filenames
extract_params_from_filename <- function(filename) {
    # Extract numbers from the filename
    pattern <- "dist_\\w+_(\\d+)_(\\d+)_range_([\\d\\.]+)"
    matches <- str_match(filename, pattern)
    
    if (is.na(matches[1])) {
        stop("Filename format does not match the expected pattern.")
    }
    
    N <- as.integer(matches[2])
    n_obs <- as.integer(matches[3])
    range <- as.numeric(matches[4])
    
    return(list(N = N, n_obs = n_obs, range = range))
}

# Process H5 Fourier data
process_h5_fourier <- function(L2, Linf, N, n_obs) {
    dist_fourier <- list()
    dist_fourier$L2 <- t(L2)
    dist_fourier$Linf <- t(Linf)    
    attr(dist_fourier, "m.vec") <- 1:6
    attr(dist_fourier, "nu.vec") <- seq(2.49, 0.01, -0.01)
    attr(dist_fourier, "N") <- N
    attr(dist_fourier, "n_obs") <- n_obs
    attr(dist_fourier, "type") <- "Fourier"
    return(dist_fourier)
}

# Process different methods and distances
process_dist <- function(df_list) {
    args <- df_list
    
    # Initialize an empty data frame for accumulating results
    dist_df <- tibble()

    for (i in 1:length(args)) {
        possible_N <- as.character(attr(args[[i]], "N"))
        n_obs <- as.character(attr(args[[i]], "n_obs"))
        range_val <- as.character(attr(args[[i]], "range"))        
        m <- as.character(attr(args[[i]], "m.vec"))
        nu <- attr(args[[i]], "nu.vec")
        method <- attr(args[[i]], "type")
        
        # Process L2 distances
        tmp_df <- as.data.frame(args[[i]][["L2"]])
        colnames(tmp_df) <- m
        tmp_df[["Method"]] <- method
        tmp_df[["nu"]] <- nu
        tmp_df[["Dist"]] <- "L2"
        tmp_df[["N"]] <- possible_N
        tmp_df[["n_obs"]] <- n_obs
        tmp_df[["range"]] <- range_val
        tmp_df <- pivot_longer(tmp_df, cols = all_of(m), names_to = "m", values_to = "Error")        
        if(method == "Rational"){
            tmp_df_tmp <- data.frame(Dist = "L2", nu = c(0.5, 1.5), Error = 1e-16, n_obs = n_obs, N = possible_N, range = range_val, m = rep(m,2), Method = "Rational")
            tmp_df_tmp <- tibble(tmp_df_tmp)
            tmp_df = bind_rows(tmp_df,tmp_df_tmp)
        }

       
        # Append the processed data to the cumulative data frame
        dist_df <- bind_rows(dist_df, tmp_df)

        tmp_df <- as.data.frame(args[[i]][["Linf"]])
        colnames(tmp_df) <- m
        tmp_df[["Method"]] <- method
        tmp_df[["nu"]] <- nu
        tmp_df[["Dist"]] <- "Linf"
        tmp_df[["N"]] <- possible_N
        tmp_df[["n_obs"]] <- n_obs        
        tmp_df[["range"]] <- range_val
        tmp_df <- pivot_longer(tmp_df, all_of(m), names_to = "m", values_to = "Error")
        if(method == "Rational"){
            tmp_df_tmp <- data.frame(Dist = "Linf", nu = c(0.5, 1.5), Error = 1e-16, n_obs = n_obs, N = possible_N, range = range_val, m = rep(m,2), Method = "Rational")
            tmp_df_tmp <- tibble(tmp_df_tmp)
            tmp_df = bind_rows(tmp_df,tmp_df_tmp)
        }        

        
        # Append the Linf processed data to the cumulative data frame
        dist_df <- bind_rows(dist_df, tmp_df)
    }

    return(dist_df)
}

process_all_files <- function(file_paths) {
    all_dfs <- list()

    for (file in file_paths) {
        cat("Processing file: ", file, "\n")
        params <- extract_params_from_filename(file)
        N <- params$N
        n_obs <- params$n_obs
        range_val <- params$range

        if (grepl(".RDS$", file)) {
            dist_data <- readRDS(file)
        } else if (grepl(".h5$", file)) {
            L2 <- h5read(file, "L2")
            Linf <- h5read(file, "Linf")
            dist_data <- process_h5_fourier(L2, Linf, N, n_obs)
        }

        attr(dist_data,"range") <- range_val
        attr(dist_data,"n_obs") <- n_obs
        all_dfs[[file]] <- dist_data
    }
    
    all_dfs <- all_dfs

    combined_df <- process_dist(all_dfs)
    combined_df <- unique(combined_df)
    
    return(combined_df)
}


# Function to plot the L2 and Linfinity distances between the covariance functions
# If methods is null, all methods will be considered
# If n_loc is null, the first value of n_loc will be considered. n_loc must be a character
# n_loc is the number of locations

library(ggplot2)

plot_dist <- function(dist_df, N_plot, n_obs_plot, range_val, m_vec = NULL, distance = c("L2", "Linf"), methods = NULL, logscale = TRUE, style=1){
    N <- as.character(N_plot)
    n_obs <- as.character(n_obs_plot)
    range <- as.character(range_val)

    if(is.null(methods)){
        methods <- unique(dist_df[["Method"]])
    }

    if(is.null(m_vec)){
        m_vec <- unique(dist_df[["m"]])
    }

    idx_plot <- dist_df$Dist == distance & dist_df$N == N_plot & dist_df$n_obs == n_obs_plot & dist_df$range == range_val & dist_df$Method %in% methods & dist_df$m %in% m_vec
    df_dist_plot <- dist_df[idx_plot,]

    if(logscale){
        if(style==1){
            return(df_dist_plot |> ggplot2::ggplot() + 
                geom_line(aes(x = nu, y = Error, col = m,linetype=Method)) + scale_y_log10())
        } else{
                        return(df_dist_plot |> ggplot2::ggplot() + 
                geom_line(aes(x = nu, y = Error, col = Method,linetype=m)) + scale_y_log10())
        }
    } else{
        return(df_dist_plot |> ggplot2::ggplot() + 
            geom_line(aes(x = nu, y = Error, col = m,linetype=Method)))
    }
}
