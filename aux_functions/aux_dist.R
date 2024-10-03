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
    
    # Extract methods, N values, and other attributes for each argument
    methods <- unlist(lapply(args, function(i){attr(i, "type")}))
    
    # Initialize an empty data frame for accumulating results
    dist_df <- data.frame()

    for (i in 1:length(args)) {
        possible_N <- as.character(attr(args[[i]], "N"))
        n_obs <- as.character(attr(args[[i]], "n_obs"))
        range <- as.character(attr(args[[i]], "range"))        
        m <- as.character(attr(args[[i]], "m.vec"))
        nu <- attr(args[[i]], "nu.vec")
        
        # Process L2 distances
        tmp_df <- as.data.frame(args[[i]][["L2"]])
        colnames(tmp_df) <- m
        tmp_df[["Method"]] <- methods[i]
        tmp_df[["nu"]] <- nu
        tmp_df[["Dist"]] <- "L2"
        tmp_df[["N"]] <- possible_N
        tmp_df[["n_obs"]] <- n_obs
        tmp_df[["range"]] <- range
        tmp_df <- pivot_longer(tmp_df, cols = all_of(m), names_to = "m", values_to = "Error")
        
        # Append the processed data to the cumulative data frame
        dist_df <- rbind(dist_df, tmp_df)
    }

    # Process Linf distances
    for (i in 1:length(args)) {
        possible_N <- as.character(attr(args[[i]], "N"))
        n_obs <- as.character(attr(args[[i]], "n_obs"))
        range <- as.character(attr(args[[i]], "range"))           
        m <- as.character(attr(args[[i]], "m.vec"))
        nu <- attr(args[[i]], "nu.vec")
        
        tmp_df <- as.data.frame(args[[i]][["Linf"]])
        colnames(tmp_df) <- m
        tmp_df[["Method"]] <- methods[i]
        tmp_df[["nu"]] <- nu
        tmp_df[["Dist"]] <- "Linf"
        tmp_df[["N"]] <- possible_N
        tmp_df[["n_obs"]] <- n_obs        
        tmp_df[["range"]] <- range
        tmp_df <- pivot_longer(tmp_df, all_of(m), names_to = "m", values_to = "Error")
        
        # Append the Linf processed data to the cumulative data frame
        dist_df <- rbind(dist_df, tmp_df)
    }

    return(dist_df)
}

process_all_files <- function(file_paths) {
    all_dfs <- list()

    for (file in file_paths) {
        params <- extract_params_from_filename(file)
        N <- params$N
        n_obs <- params$n_obs
        range <- params$range
        print(range)

        if (grepl(".RDS$", file)) {
            dist_data <- readRDS(file)
        } else if (grepl(".h5$", file)) {
            L2 <- h5read(file, "L2")
            Linf <- h5read(file, "Linf")
            dist_data <- process_h5_fourier(L2, Linf, N, n_obs)
        }

        attr(dist_data,"range") <- range
        all_dfs[[file]] <- dist_data
    }
    
    all_dfs <<- all_dfs

    combined_df <- process_dist(all_dfs)
    
    return(combined_df)
}


# Function to plot the L2 and Linfinity distances between the covariance functions
# If methods is null, all methods will be considered
# If n_loc is null, the first value of n_loc will be considered. n_loc must be a character
# n_loc is the number of locations

library(ggplot2)

plot_dist <- function(df_dist, n_loc, n_obs, range, distance = c("L2", "Linf"), methods = NULL, logscale = TRUE, style=1){
    n_loc <- as.character(n_loc)
    n_obs <- as.character(n_obs)
    df_dist <- df_dist |> dplyr::filter(N == n_loc, n_obs == n_obs)
    if(is.null(methods)){
        methods <- unique(df_dist[["Method"]])
    }
    distance <- distance[[1]]
    if(logscale){
        if(style==1){
            return(df_dist |> dplyr::filter(Dist == distance, range == range, Method %in% methods) |> ggplot2::ggplot() + 
                geom_line(aes(x = nu, y = Error, col = m,linetype=Method)) + scale_y_log10())
        } else{
                        return(df_dist |> dplyr::filter(Dist == distance, range == range, Method %in% methods) |> ggplot2::ggplot() + 
                geom_line(aes(x = nu, y = Error, col = Method,linetype=m)) + scale_y_log10())
        }
    } else{
        return(df_dist |> dplyr::filter(Dist == distance,range == range, Method %in% methods) |> ggplot2::ggplot() + 
            geom_line(aes(x = nu, y = Error, col = m,linetype=Method)))
    }
}
