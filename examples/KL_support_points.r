library(rSPDE)

source("aux_functions/prediction_computations.R")

apply_matern_by_row <- function(h, kappa, sigma, nu, order) {
  if(is.matrix(h)){
    n <- nrow(h)
    result <- matrix(0, nrow=n, ncol=n)    
    for(i in 1:n) {
      row_dist <- h[i,]
      ord_idx <- order(row_dist)
      ordered_dist <- row_dist[ord_idx]
      ordered_result <- matern.rational.cov(ordered_dist, kappa=kappa, sigma=sigma, nu=nu, order=order)
      result[i,] <- ordered_result[order(ord_idx)]
    }
  } else{
    result <- matern.rational.cov(h, kappa=kappa, sigma=sigma, nu=nu, order=order)
  }  
  return(result)
}


build_matern_rational_cov <- function(loc, kappa, sigma, nu, order){
  # if(nu < 2 && order < 3){
    return(apply_matern_by_row(as.matrix(dist(loc)), kappa, sigma, nu, order))
  # } else{
  #     Qrat <- rSPDE:::matern.rational.precision(
  #       loc           = loc,
  #       order         = order,
  #       nu            = nu,
  #       kappa         = kappa,
  #       sigma         = sigma,
  #       type_rational = "brasil",
  #       type_interp   = "spline"
  #     )
  #     Sigma_rat <- as.matrix(Qrat$A %*% solve(Qrat$Q, t(Qrat$A)))
  #     return(Sigma_rat)
  # }
}

# KL(Rational, TRUE)
compute_kl_divergence_rational <- function(loc, m.vec, nu.vec, range, sigma) {
    print("Computing KL divergence for rational approximation")
    kl.div <- matrix(0, length(nu.vec), length(m.vec))
    
    print(paste("Domain length = ", max(loc)))
    print(paste("range = ", range))
    
    for(i in 1:length(nu.vec)) {
        nu <- nu.vec[i]
        alpha <- nu + 0.5
        kappa <- sqrt(8*nu)/range
        
        # Get true Matérn covariance
        Sigma.t <- rSPDE::matern.covariance(h=as.matrix(dist(loc)), kappa=kappa, nu=nu, sigma=sigma)
        
        log_det_true <- determinant(Sigma.t, logarithm=TRUE)$modulus
        
        for(j in 1:length(m.vec)) {
            m <- m.vec[j]
            
            Sigma_rat <- build_matern_rational_cov(loc = loc, kappa = kappa, sigma = sigma, nu = nu, order = m)
            
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
    
    return(ret)
}


# KL(NNGP, True)
compute_kl_divergence_nngp <- function(loc, obs.ind, m.vec, nu.vec, range, sigma) {
    print("Computing KL divergence for NNGP approximation")
    kl.div <- matrix(0, length(nu.vec), length(m.vec))
    
    print(paste("Domain length = ", max(loc)))
    print(paste("range = ", range))
    
    for(i in 1:length(nu.vec)) {
        nu <- nu.vec[i]
        alpha <- nu + 0.5
        kappa <- sqrt(8*nu)/range
        
        # Get true Matérn covariance
        Sigma.t <- rSPDE::matern.covariance(h=as.matrix(dist(loc)), kappa=kappa, nu=nu, sigma=sigma)
        
        log_det_true <- determinant(Sigma.t, logarithm=TRUE)$modulus
        
        for(j in 1:length(m.vec)) {
            m <- m.vec[j]
            N_m <- length(obs.ind)
            # i_m <- m_nngp_fun(m, nu + 0.5, N_m, N_m)
            i_m <- m_nngp_fun(m, nu + 0.5, 1000, 1000)
            Sigma_inv_nngp <- get.nnQ(loc=loc,kappa=kappa,nu=nu,sigma=sigma, n.nbr = i_m, S = obs.ind)
            
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
    
    return(ret)
}

# Example
n            <- c(0:10, seq(20,50,by=10), seq(75,250,by=25), seq(300,500, by = 50))
n            <- rev(n)        # Reverse the vector
n.obs        <- 1001
domain_upper_limit <- (n.obs + max(n)-1)/100
full_mesh <- seq(0, domain_upper_limit, length.out = n.obs+max(n))
obs_loc <- full_mesh[1:n.obs]

locs <- lapply(n, function(i) {
  full_mesh[1:(n.obs + i)]
})

nu.vec <- 1
m.vec <- 1:6
range_val <- 2
sigma <- 1

# results_rational <- compute_kl_divergence_rational(loc = locs[[1]], m.vec=m.vec, 
#                                         nu.vec=nu.vec, range=range_val, sigma=sigma)

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


# results_nngp <- compute_kl_divergence_nngp(loc = loc, obs.ind = 1:n.obs, m.vec=m.vec, 
#                                         nu.vec=nu.vec, range=range_val, sigma=sigma)
