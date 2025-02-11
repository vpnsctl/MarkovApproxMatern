
source("aux_functions/aux_functions_cov.R")


# # KL(Rational, TRUE)
# compute_kl_divergence_rational <- function(N, m.vec, nu.vec, range, sigma) {
#     print("Computing KL divergence for rational approximation")
#     N <- N[[1]]
#     loc <- seq(0, N/100, length.out = N)
#     D <- as.matrix(dist(loc))
#     kl.div <- matrix(0, length(nu.vec), length(m.vec))
    
#     print(paste("N = ", N))
#     print(paste("Domain length = ", max(loc)))
#     print(paste("range = ", range))
    
#     for(i in 1:length(nu.vec)) {
#         nu <- nu.vec[i]
#         alpha <- nu + 0.5
#         kappa <- sqrt(8*nu)/range
        
#         # Get true Matérn covariance
#         Sigma.t <- rSPDE::matern.covariance(h=D, kappa=kappa, nu=nu, sigma=sigma)
        
#         log_det_true <- determinant(Sigma.t, logarithm=TRUE)$modulus
        
#         for(j in 1:length(m.vec)) {
#             m <- m.vec[j]
            
#             # Compute rational approximation covariance
#             col_tmp <- rSPDE:::matern.rational.cov(h=loc, order=m, kappa=kappa, 
#                                                          nu=nu, sigma=sigma, 
#                                                          type_rational="brasil", 
#                                                          type_interp="spline")    
#             Sigma_rat <- toeplitz(x=drop(col_tmp), symmetric=TRUE)
            
#             tryCatch({
#                 # Compute log determinant of rational approximation
#                 Sigma_rat_inv <- solve(Sigma_rat)
#                 log_det_rat <- determinant(Sigma_rat, logarithm=TRUE)$modulus

#                 temp <- Sigma_rat_inv %*% Sigma.t
#                 trace_term <- sum(diag(temp))
#                 dim <- nrow(Sigma.t)
#                 kl.rat <- as.double(0.5 * (trace_term - dim  + log_det_rat - log_det_true))

#                 kl.div[i,j] <- kl.rat
#                 # Print partial KL components
#                 cat("\nKL components for nu =", nu, "m =", m, ":\n")
#                 cat("KL divergence:", kl.rat, "\n\n")                
#             }, error=function(e) {
#                 warning(paste("Error computing KL divergence for nu =", nu, "m =", m))
#                 kl.div[i,j] <- NA
#             })
#         }
#     }
    
#     ret <- list(KL=kl.div)
#     attr(ret, "type") <- "Rational"
#     attr(ret, "nu.vec") <- nu.vec
#     attr(ret, "m.vec") <- m.vec
#     attr(ret, "N") <- N
    
#     return(ret)
# }



# KL(Rational, TRUE)
compute_kl_divergence_rational <- function(N, m.vec, nu.vec, range, sigma) {
    print("Computing KL divergence for rational approximation")
    N <- N[[1]]
    loc <- seq(0, N/100, length.out = N)
    D <- as.matrix(dist(loc))
    kl.div <- matrix(0, length(nu.vec), length(m.vec))
    
    print(paste("N = ", N))
    print(paste("Domain length = ", max(loc)))
    print(paste("range = ", range))
    
    for(i in 1:length(nu.vec)) {
        nu <- nu.vec[i]
        alpha <- nu + 0.5
        kappa <- sqrt(8*nu)/range
        
        # Get true Matérn covariance
        Sigma.t <- rSPDE::matern.covariance(h=D, kappa=kappa, nu=nu, sigma=sigma)
        
        log_det_true <- determinant(Sigma.t, logarithm=TRUE)$modulus
        
        for(j in 1:length(m.vec)) {
            m <- m.vec[j]
            
            # Compute rational approximation covariance
            Qrat <-rSPDE:::matern.rational.ldl(loc = loc, order = m, nu = nu, kappa = kappa, 
                                               sigma = sigma, type_rational = "brasil", 
                                               type_interp =  "spline", equally_spaced = FALSE)    
            
            Q <- t(Qrat$L)%*%Qrat$D%*%Qrat$L
            Sigma_rat <- Qrat$A %*% solve(Q, t(Qrat$A))
            Sigma_rat <- as.matrix(Sigma_rat)
            
            tryCatch({
                # Compute log determinant of rational approximation
                Sigma_rat_inv <- solve(Sigma_rat)
                log_det_rat <- determinant(Sigma_rat, logarithm=TRUE)$modulus

                temp <- Sigma_rat_inv %*% Sigma.t
                trace_term <- sum(diag(temp))
                dim <- nrow(Sigma.t)
                kl.rat <- as.double(0.5 * (trace_term - dim  + log_det_rat - log_det_true))

                kl.div[i,j] <- kl.rat
                # Print partial KL components
                cat("\nKL components for nu =", nu, "m =", m, ":\n")
                cat("KL divergence:", kl.rat, "\n\n")                
            }, error=function(e) {
                warning(paste("Error computing KL divergence for nu =", nu, "m =", m))
                kl.div[i,j] <- NA
            })
        }
    }
    
    ret <- list(KL=kl.div)
    attr(ret, "type") <- "Rational"
    attr(ret, "nu.vec") <- nu.vec
    attr(ret, "m.vec") <- m.vec
    attr(ret, "N") <- N
    
    return(ret)
}


# ## KL(True, Rational)
# compute_kl_divergence_rational <- function(N, m.vec, nu.vec, range, sigma) {
#     print("Computing KL divergence for rational approximation")
#     N <- N[[1]]
#     loc <- seq(0, N/100, length.out = N)
#     D <- as.matrix(dist(loc))
#     kl.div <- matrix(0, length(nu.vec), length(m.vec))
    
#     print(paste("N = ", N))
#     print(paste("Domain length = ", max(loc)))
#     print(paste("range = ", range))
    
#     for(i in 1:length(nu.vec)) {
#         nu <- nu.vec[i]
#         kappa <- sqrt(8*nu)/range
        
#         # Get true Matérn covariance
#         Sigma.t <- rSPDE::matern.covariance(h=D, kappa=kappa, nu=nu, sigma=sigma)
        
#         log_det_true <- determinant(Sigma.t, logarithm=TRUE)$modulus

#         Sigma.t_inv <- solve(Sigma.t)
        
#         for(j in 1:length(m.vec)) {
#             m <- m.vec[j]
            

#             col_tmp <- rSPDE:::matern.rational.cov(h=loc, order=m, kappa=kappa, 
#                                                          nu=nu, sigma=sigma, 
#                                                          type_rational="chebfun", 
#                                                          type_interp="spline")    
#             Sigma_rat <- toeplitz(x=drop(col_tmp), symmetric=TRUE)
            
#             tryCatch({
#                 # Compute log determinant of rational approximation

#                 log_det_rat <- determinant(Sigma_rat, logarithm=TRUE)$modulus

#                 temp <- Sigma.t_inv %*% Sigma_rat
#                 trace_term <- sum(diag(temp))
                
#                 # Compute KL divergence
#                 dim <- nrow(Sigma.t)
#                 kl <- 0.5 * (trace_term - dim  + log_det_true - log_det_rat)
#                 kl.div[i,j] <- kl
#                 # Print partial KL components
#                 cat("\nKL components for nu =", nu, "m =", m, ":\n")
#                 cat("KL divergence:", kl, "\n\n")                
#             }, error=function(e) {
#                 warning(paste("Error computing KL divergence for nu =", nu, "m =", m))
#                 kl.div[i,j] <- NA
#             })
#         }
#     }
    
#     ret <- list(KL=kl.div)
#     attr(ret, "type") <- "Rational"
#     attr(ret, "nu.vec") <- nu.vec
#     attr(ret, "m.vec") <- m.vec
#     attr(ret, "N") <- N
    
#     return(ret)
# }

# Example
nu.vec <- seq(2.49, 0.01, by = -0.01)
m.vec <- 1:6
N <- 500
range_val <- 2
sigma <- 1

results <- compute_kl_divergence_rational(N=N, m.vec=m.vec, 
                                        nu.vec=nu.vec, range=range_val, sigma=sigma)


compute_hellinger_distance_rational <- function(N, m.vec, nu.vec, range, sigma) {
    print("Computing Hellinger distance for rational approximation")
    N <- N[[1]]
    loc <- seq(0, N/100, length.out = N)
    D <- as.matrix(dist(loc))
    hellinger.div <- matrix(0, length(nu.vec), length(m.vec))
    
    print(paste("N = ", N))
    print(paste("Domain length = ", max(loc)))
    print(paste("range = ", range))
    
    for(i in 1:length(nu.vec)) {
        nu <- nu.vec[i]
        alpha <- nu + 0.5
        kappa <- sqrt(8*nu)/range
        
        # Get true Matérn covariance
        Sigma.t <- rSPDE::matern.covariance(h=D, kappa=kappa, nu=nu, sigma=sigma)
        log_det_true <- determinant(Sigma.t, logarithm=TRUE)$modulus        
        
        for(j in 1:length(m.vec)) {
            m <- m.vec[j]
            
            col_tmp <- rSPDE:::matern.rational.cov(h=loc, order=m, kappa=kappa, 
                                                         nu=nu, sigma=sigma, 
                                                         type_rational="brasil", 
                                                         type_interp="spline")    
            Sigma_rat <- toeplitz(x=drop(col_tmp), symmetric=TRUE)
            
            tryCatch({               
                log_det_rat <- determinant(Sigma_rat, logarithm=TRUE)$modulus
                log_det_sum <- determinant((Sigma.t + Sigma_rat)/2, logarithm=TRUE)$modulus
                
                # Compute Hellinger distance using log determinants
                # H^2 = 1 - exp(0.25*log|Σ₁| + 0.25*log|Σ₂| - 0.5*log|(Σ₁+Σ₂)/2|)
                hellinger <- sqrt(1 - exp(0.25*log_det_true + 0.25*log_det_rat - 0.5*log_det_sum))
                hellinger.div[i,j] <- hellinger
                
                # Print partial components
                cat("\nHellinger components for nu =", nu, "m =", m, ":\n")
                cat("Hellinger distance:", hellinger, "\n\n")
                
            }, error=function(e) {
                warning(paste("Error computing Hellinger distance for nu =", nu, "m =", m))
                hellinger.div[i,j] <- NA
            })
        }
    }
    
    ret <- list(Hellinger=hellinger.div)
    attr(ret, "type") <- "Rational"
    attr(ret, "nu.vec") <- nu.vec
    attr(ret, "m.vec") <- m.vec
    attr(ret, "N") <- N
    
    return(ret)
}

# Example usage
nu.vec <- seq(2.49, 0.01, -0.01)
nu.vec <- 1.3
m.vec <- 1:6
N <- 1000
range_val <- 2
sigma <- 1

results <- compute_hellinger_distance_rational(N=N, m.vec=m.vec, 
                                            nu.vec=nu.vec, range=range_val, sigma=sigma)





###############################
######### NNGP ################
###############################


source("aux_functions/aux_functions_cov.R")

m_nngp_fun <- function(m, alpha, n, n.obs){
            if(alpha<1) {
                mn <- m - 1
                if(mn < 1){
                    mn <- 1
                }
            } else if (alpha < 2) {
                if(n == 5000){
                    m_vec <- c(1, 2, 13, 21, 27, 31) 
                } else if(n.obs == 5000){
                    m_vec <- c(1, 2, 3, 4, 14, 20)
                } else if(n.obs == 10000){
                    m_vec <- c(1, 2, 7, 18, 24, 29)
                } else if(n <= 1000){
                    m_vec <-  c(1, 2, 8, 14, 20, 24)
                } else{
                    stop("not implemented")
                }
                mn <- m_vec[m]
            } else {
                if(n == 5000){
                    m_vec <- c(15, 30, 37, 45, 51, 54)
                } else if(n.obs == 5000){
                    m_vec <- c(1, 22, 31, 39, 47, 54)
                } else if(n.obs == 10000){
                    m_vec <- c(14, 28, 37, 44, 51, 57)
                } else if(n <= 1000){
                    m_vec <- c(4, 17, 25, 29, 34, 40)
                } else{
                    stop("not implemented")
                }
                mn <- m_vec[m]
            }
            return(mn)
} 


# KL(NNGP, True)
compute_kl_divergence_nngp <- function(N, m.vec, nu.vec, range, sigma) {
    print("Computing KL divergence for NNGP approximation")
    N <- N[[1]]
    loc <- seq(0, N/100, length.out = N)
    D <- as.matrix(dist(loc))
    kl.div <- matrix(0, length(nu.vec), length(m.vec))
    
    print(paste("N = ", N))
    print(paste("Domain length = ", max(loc)))
    print(paste("range = ", range))
    
    for(i in 1:length(nu.vec)) {
        nu <- nu.vec[i]
        alpha <- nu + 0.5
        kappa <- sqrt(8*nu)/range
        
        # Get true Matérn covariance
        Sigma.t <- rSPDE::matern.covariance(h=D, kappa=kappa, nu=nu, sigma=sigma)
        
        log_det_true <- determinant(Sigma.t, logarithm=TRUE)$modulus
        
        for(j in 1:length(m.vec)) {
            m <- m.vec[j]
            i_m <- m_nngp_fun(m, nu + 0.5, N, N)
            Sigma_inv_nngp <- get.nnQ(loc=loc,kappa=kappa,nu=nu,sigma=sigma, n.nbr = i_m)
            
            tryCatch({
                # Compute log determinant of rational approximation

                log_det_inv_nngp <- determinant(Sigma_inv_nngp, logarithm=TRUE)$modulus

                temp <- Sigma_inv_nngp %*% Sigma.t
                trace_term <- sum(diag(temp))
                
                # Compute KL divergence
                dim <- nrow(Sigma.t)
                kl <- 0.5 * (trace_term - dim - log_det_inv_nngp - log_det_true)
                kl.div[i,j] <- kl
                # Print partial KL components
                cat("\nKL components for nu =", nu, "m =", m, ":\n")
                cat("KL divergence:", kl, "\n\n")                
            }, error=function(e) {
                warning(paste("Error computing KL divergence for nu =", nu, "m =", m))
                kl.div[i,j] <- NA
            })
        }
    }
    
    ret <- list(KL=kl.div)
    attr(ret, "type") <- "NNGP"
    attr(ret, "nu.vec") <- nu.vec
    attr(ret, "m.vec") <- m.vec
    attr(ret, "N") <- N
    
    return(ret)
}


# ## KL(True, NNGP)
# compute_kl_divergence_nngp <- function(N, m.vec, nu.vec, range, sigma, invert_prec=TRUE) {
#     print("Computing KL divergence for NNGP approximation")
#     N <- N[[1]]
#     loc <- seq(0, N/100, length.out = N)
#     D <- as.matrix(dist(loc))
#     kl.div <- matrix(0, length(nu.vec), length(m.vec))
    
#     print(paste("N = ", N))
#     print(paste("Domain length = ", max(loc)))
#     print(paste("range = ", range))
    
#     for(i in 1:length(nu.vec)) {
#         nu <- nu.vec[i]
#         alpha <- nu + 0.5
#         kappa <- sqrt(8*nu)/range
        
#         # Get true Matérn covariance
#         Sigma.t <- rSPDE::matern.covariance(h=D, kappa=kappa, nu=nu, sigma=sigma)
        
#         log_det_true <- determinant(Sigma.t, logarithm=TRUE)$modulus

#         Sigma.t.inv <- solve(Sigma.t)
        
#         for(j in 1:length(m.vec)) {
#             m <- m.vec[j]
#             i_m <- m_nngp_fun(m, nu + 0.5, N, N)
#             Sigma_nngp <- get.nnCov(loc=loc,kappa=kappa,nu=nu,sigma=sigma, n.nbr=m, invert_prec = invert_prec)
            
#             tryCatch({
#                 # Compute log determinant of rational approximation

#                 log_det_nngp <- determinant(Sigma_nngp, logarithm=TRUE)$modulus

#                 temp <- Sigma.t.inv %*% Sigma_nngp
#                 trace_term <- sum(diag(temp))
                
#                 # Compute KL divergence
#                 dim <- nrow(Sigma.t)
#                 kl <- 0.5 * (trace_term - dim + log_det_true- log_det_nngp)
#                 kl.div[i,j] <- kl
#                 # Print partial KL components
#                 cat("\nKL components for nu =", nu, "m =", m, ":\n")
#                 cat("KL divergence:", kl, "\n\n")                
#             }, error=function(e) {
#                 warning(paste("Error computing KL divergence for nu =", nu, "m =", m))
#                 kl.div[i,j] <- NA
#             })
#         }
#     }
    
#     ret <- list(KL=kl.div)
#     attr(ret, "type") <- "NNGP"
#     attr(ret, "nu.vec") <- nu.vec
#     attr(ret, "m.vec") <- m.vec
#     attr(ret, "N") <- N
    
#     return(ret)
# }

# Example
nu.vec <- seq(2.49, 0.01, -0.01)
m.vec <- 1:6
N <- 500
range_val <- 2
sigma <- 1

results <- compute_kl_divergence_nngp(N=N, m.vec=m.vec, 
                                        nu.vec=nu.vec, range=range_val, sigma=sigma)



compute_hellinger_distance_nngp <- function(N, m.vec, nu.vec, range, sigma, invert_prec = TRUE) {
    print("Computing Hellinger distance for NNGP approximation")
    N <- N[[1]]
    loc <- seq(0, N/100, length.out = N)
    D <- as.matrix(dist(loc))
    hellinger.div <- matrix(0, length(nu.vec), length(m.vec))
    
    print(paste("N = ", N))
    print(paste("Domain length = ", max(loc)))
    print(paste("range = ", range))
    
    for(i in 1:length(nu.vec)) {
        nu <- nu.vec[i]
        alpha <- nu + 0.5
        kappa <- sqrt(8*nu)/range
        
        # Get true Matérn covariance
        Sigma.t <- rSPDE::matern.covariance(h=D, kappa=kappa, nu=nu, sigma=sigma)
        log_det_true <- determinant(Sigma.t, logarithm=TRUE)$modulus        
        
        for(j in 1:length(m.vec)) {
            m <- m.vec[j]
            i_m <- m_nngp_fun(m, nu + 0.5, N, N)
            Sigma_nngp <- get.nnCov(loc=loc,kappa=kappa,nu=nu,sigma=sigma, n.nbr=m, invert_prec = invert_prec)
            
            tryCatch({                
                # Compute Hellinger distance using the formula:
                log_det_nngp <- determinant(Sigma_nngp, logarithm=TRUE)$modulus
                log_det_sum <- determinant((Sigma.t + Sigma_nngp)/2, logarithm=TRUE)$modulus
                
                # Compute Hellinger distance using log determinants
                # H^2 = 1 - exp(0.25*log|Σ₁| + 0.25*log|Σ₂| - 0.5*log|(Σ₁+Σ₂)/2|)
                hellinger <- sqrt(1 - exp(0.25*log_det_true + 0.25*log_det_nngp - 0.5*log_det_sum))
                hellinger.div[i,j] <- hellinger
                
                # Print partial components
                cat("\nHellinger components for nu =", nu, "m =", m, ":\n")
                cat("Hellinger distance:", hellinger, "\n\n")
                
            }, error=function(e) {
                warning(paste("Error computing Hellinger distance for nu =", nu, "m =", m))
                hellinger.div[i,j] <- NA
            })
        }
    }
    
    ret <- list(Hellinger=hellinger.div)
    attr(ret, "type") <- "NNGP"
    attr(ret, "nu.vec") <- nu.vec
    attr(ret, "m.vec") <- m.vec
    attr(ret, "N") <- N
    
    return(ret)
}

# Example usage
nu.vec <- 1.3
m.vec <- 1:6
N <- 5000
range_val <- 2
sigma <- 1

results <- compute_hellinger_distance_nngp(N=N, m.vec=m.vec, 
                                         nu.vec=nu.vec, range=range_val, sigma=sigma)




###############################
########## FEM ################
###############################

source("aux_functions/aux_functions_cov.R")

m_fem_fun <- function(m, alpha, n, n.obs){
            if(alpha<1) {
                if(n == 5000){
                m_vec <- c(1, 1, 2, 3, 3, 3)
                } else if(n == 1){
                m_vec <- 
                } else {
                    stop("not implemented")
                }
                mn <- m_vec[m]                
            } else if (alpha < 2) {
                if(n == 5000){
                    m_vec <- c(2, 4, 6, 7, 8, 9) 
                } else{
                    stop("not implemented")
                }
                mn <- m_vec[m]
            } else {
                if(n == 5000){
                    # m_vec <- c(7, 13, 17, 21, 21, 21)
                    m_vec <- c(5,4,3,2,2,1)  # Chose to stabilize the results, not based on cost.
                } else{
                    stop("not implemented")
                }
                mn <- m_vec[m]
            }
            return(mn)
} 


## KL(True, FEM)
compute_kl_divergence_fem <- function(N, m.vec, nu.vec, range, sigma) {
    print("Computing KL divergence for FEM approximation")
    N <- N[[1]]
    loc <- seq(0, N/100, length.out = N)
    D <- as.matrix(dist(loc))
    kl.div <- matrix(0, length(nu.vec), length(m.vec))
    
    print(paste("N = ", N))
    print(paste("Domain length = ", max(loc)))
    print(paste("range = ", range))
    
    for(i in 1:length(nu.vec)) {
        nu <- nu.vec[i]
        kappa <- sqrt(8*nu)/range
        
        # Get true Matérn covariance
        Sigma.t <- rSPDE::matern.covariance(h=D, kappa=kappa, nu=nu, sigma=sigma)
        
        log_det_true <- determinant(Sigma.t, logarithm=TRUE)$modulus
        
        for(j in 1:length(m.vec)) {
            m_rat <- m.vec[j]

            alpha <- nu + 0.5

            m <- m_fem_fun(m_rat, alpha, n = N, n.obs = N_obs)
            
            Sigma_fem <- fem_cov(m = m_rat, mesh_fem = m, loc = loc, kappa = kappa, nu = nu, sigma = sigma)

            Sigma_fem <- as.matrix(Sigma_fem)
            
            tryCatch({
                # Compute log determinant of rational approximation

                log_det_fem <- determinant(Sigma_fem, logarithm=TRUE)$modulus

                temp <- solve(Sigma_fem, Sigma.t)
                trace_term <- sum(diag(temp))
                
                # Compute KL divergence
                dim <- nrow(Sigma.t)
                kl <- 0.5 * (trace_term - dim + log_det_fem - log_det_true)
                kl.div[i,j] <- kl
                # Print partial KL components
                cat("\nKL components for nu =", nu, "m =", m, ":\n")
                cat("KL divergence:", kl, "\n\n")                
            }, error=function(e) {
                warning(paste("Error computing KL divergence for nu =", nu, "m =", m))
                kl.div[i,j] <- NA
            })
        }
    }
    
    ret <- list(KL=kl.div)
    attr(ret, "type") <- "FEM"
    attr(ret, "nu.vec") <- nu.vec
    attr(ret, "m.vec") <- m.vec
    attr(ret, "N") <- N
    
    return(ret)
}

# Example
nu.vec <- seq(2.49, 0.01, by = -0.01)
m.vec <- 1:6
N <- 1000
range_val <- 2
sigma <- 1

results <- compute_kl_divergence_fem(N=N, m.vec=m.vec, 
                                        nu.vec=nu.vec, range=range_val, sigma=sigma)




library(ggplot2)
library(tidyr)
library(dplyr)

plot_comparison <- function(rational_file, nngp_file) {
  # Read the saved results
  rational_results <- readRDS(rational_file)
  nngp_results <- readRDS(nngp_file)
  
  # Get attributes
  nu_vec <- attr(rational_results, "nu.vec")
  m_vec <- attr(rational_results, "m.vec")
  rational_type <- attr(rational_results, "type")
  nngp_type <- attr(nngp_results, "type")
  
  # Create data frames for each type
  rational_df <- as.data.frame(rational_results$KL)
  nngp_df <- as.data.frame(nngp_results$KL)
  
  # Add nu values and convert to long format
  rational_df$nu <- nu_vec
  nngp_df$nu <- nu_vec
  
  # Rename columns to match m values
  names(rational_df)[1:length(m_vec)] <- paste0("m", m_vec)
  names(nngp_df)[1:length(m_vec)] <- paste0("m", m_vec)
  
  # Convert to long format
  rational_long <- pivot_longer(rational_df, 
                              cols = starts_with("m"),
                              names_to = "m",
                              values_to = "KL") %>%
                  mutate(type = rational_type)
  
  nngp_long <- pivot_longer(nngp_df,
                           cols = starts_with("m"),
                           names_to = "m",
                           values_to = "KL") %>%
                  mutate(type = nngp_type)
  
  # Combine the data
  plot_data <- rbind(rational_long, nngp_long)
  
  # Create the plot
  ggplot(plot_data, aes(x = nu, y = KL, color = m, linetype = type)) +
    geom_line() +
    scale_y_log10() +
    scale_linetype_manual(values = c("solid", "dashed")) +
    labs(title = paste("Comparison of", rational_type, "vs", nngp_type),
         x = "nu",
         y = "KL Divergence (log scale)") +
    theme_minimal() +
    theme(legend.position = "right",
          legend.box = "vertical",
          legend.margin = margin()) +
    guides(color = guide_legend(title = "Parameters"),
           linetype = guide_legend(title = "Model Type"))
}

# Usage:
plot_comparison("rational_results.rds", "nngp_results.rds")