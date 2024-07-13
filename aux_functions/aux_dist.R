# Function to process the distance tables
library(dplyr)
library(tidyr)

process_dist <- function(...){
    args <- list(...)
    possible_N <- names(args[[1]][["L2"]])
    methods <- unlist(lapply(args, function(i){attr(i, "type")}))
    m <- attr(args[[1]], "m.vec")
    nu <- attr(args[[1]], "nu.vec")
    dist_df <- as.data.frame(args[[1]][["L2"]][[possible_N[1]]])
    colnames(dist_df) <- m
    dist_df[["Method"]] <- methods[1]
    dist_df[["nu"]] <- nu
    dist_df[["Dist"]] <- "L2"
    dist_df[["N"]] <- possible_N[1]
    dist_df <- pivot_longer(dist_df, all_of(m), names_to = "m", values_to = "Error")
    if(length(args)>1){
        for(i in 2:length(args)){
            tmp_df <- as.data.frame(args[[i]][["L2"]][[possible_N[1]]])
            m <- attr(args[[i]], "m.vec")
            nu <- attr(args[[i]], "nu.vec")            
            colnames(tmp_df) <- m
            tmp_df[["Method"]] <- methods[i]
            tmp_df[["nu"]] <- nu
            tmp_df[["Dist"]] <- "L2"
            tmp_df[["N"]] <- possible_N[1]            
            tmp_df <- pivot_longer(tmp_df, all_of(m), names_to = "m", values_to = "Error")            
            dist_df <- rbind(dist_df, tmp_df)
        }
    }
    if(length(possible_N)>1){
        for(i in 1:length(args)){
            for(j in 2:length(possible_N)){
                tmp_df <- as.data.frame(args[[i]][["L2"]][[possible_N[j]]])
                m <- attr(args[[i]], "m.vec")
                nu <- attr(args[[i]], "nu.vec")                    
                colnames(tmp_df) <- m
                tmp_df[["Method"]] <- methods[i]
                tmp_df[["nu"]] <- nu
                tmp_df[["Dist"]] <- "L2"
                tmp_df[["N"]] <- possible_N[j]                  
                tmp_df <- pivot_longer(tmp_df, all_of(m), names_to = "m", values_to = "Error")            
                dist_df <- rbind(dist_df, tmp_df)
            }
        }
    }

    for(i in 1:length(args)){
        for(j in 1:length(possible_N)){
            tmp_df <- as.data.frame(args[[i]][["Linf"]][[possible_N[j]]])
            m <- attr(args[[i]], "m.vec")
            nu <- attr(args[[i]], "nu.vec")                
            colnames(tmp_df) <- m
            tmp_df[["Method"]] <- methods[i]
            tmp_df[["nu"]] <- nu
            tmp_df[["Dist"]] <- "Linf"
            tmp_df[["N"]] <- possible_N[j]                  
            tmp_df <- pivot_longer(tmp_df, all_of(m), names_to = "m", values_to = "Error")            
            dist_df <- rbind(dist_df, tmp_df)
        }
    }
    return(dist_df)
}


# Function to plot the L2 and Linfinity distances between the covariance functions
# If methods is null, all methods will be considered
# If n_loc is null, the first value of n_loc will be considered. n_loc must be a character
# n_loc is the number of locations

library(ggplot2)

plot_dist <- function(df_dist, n_loc = NULL, distance = c("L2", "Linf"), methods = NULL, logscale = TRUE, style=1){
    if(is.null(n_loc)){
        n_loc <- df_dist[["N"]][1]
    }
    n_loc <- as.character(n_loc)
    df_dist <- df_dist |> dplyr::filter(N == n_loc)
    if(is.null(methods)){
        methods <- unique(df_dist[["Method"]])
    }
    distance <- distance[[1]]
    if(logscale){
        if(style==1){
            return(df_dist |> dplyr::filter(Dist == distance, Method %in% methods) |> ggplot2::ggplot() + 
                geom_line(aes(x = nu, y = Error, col = m,linetype=Method)) + scale_y_log10())
        } else{
                        return(df_dist |> dplyr::filter(Dist == distance, Method %in% methods) |> ggplot2::ggplot() + 
                geom_line(aes(x = nu, y = Error, col = Method,linetype=m)) + scale_y_log10())
        }
    } else{
        return(df_dist |> dplyr::filter(Dist == distance, Method %in% methods) |> ggplot2::ggplot() + 
            geom_line(aes(x = nu, y = Error, col = m,linetype=Method)))
    }
}
