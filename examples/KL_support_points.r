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

build_locations_with_support <- function(n.obs, n.pred, domain_upper_limit = max(n.pred)/100, tol = 1e-4){
  locs <- lapply(n.pred, function(x) {
    if(x > 0) {
      seq(from = 1, to = domain_upper_limit+5, length.out = x) 
    } else {
      NULL
    }
  })
  obs_loc <- sort(runif(n.obs, 0, domain_upper_limit))
  while (min(diff(obs_loc)) < tol) {
    obs_loc <- sort(runif(n.obs, 0, domain_upper_limit))
  }

  full_locs <- list()
  for(i in 1:length(n.pred)){
        obs.ind <- numeric(n.obs)
    if(n.pred[[i]] > 0){
      loc <- c(obs_loc, locs[[i]])
      tmp <- sort(loc, index.return = TRUE)
      loc <- tmp$x
      obs.ind[tmp$ix] <- seq_along(loc)
      obs.ind <- obs.ind[1:n.obs]
      if(any(diff(loc) < tol)) {
        ind_remove <- which(diff(loc) < tol)
        ind_remove <- unique(c(ind_remove, ind_remove + 1))
        # Don't remove actual observation points:
        ind_remove <- setdiff(ind_remove, obs.ind)
        loc <- loc[-ind_remove]
        # Re-match obs.ind by closeness
        tolerance <- 1e-6
        obs.ind <- sapply(obs_loc, function(x) {
          which(abs(loc - x) < tolerance)[1]
        })
      }      
    } else {
      # If n = 0, only observation locations
      loc     <- obs_loc
      obs.ind <- seq_len(n.obs)
    }
    full_locs[[i]] <- list(loc = loc, obs.ind = obs.ind, domain_upper_limit = domain_upper_limit)
  }
  return(full_locs)
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
locs <- build_locations_with_support(1000, c(0,1,2,4,8,16,32,64,128,512,1024), domain_upper_limit = 50, tol = 1e-4)

nu.vec <- seq(2.49, 0.01, by = -0.01)

nu.vec <- 1
m.vec <- 1:6
range_val <- 2
sigma <- 1

n <- 0

loc <- sort(runif(n.obs + n, 0, domain_upper_limit))
while (min(diff(loc)) < tol) {
  loc <- sort(runif(n.obs + n, 0, domain_upper_limit))
}

results_rational <- compute_kl_divergence_rational(loc = loc, m.vec=m.vec, 
                                        nu.vec=nu.vec, range=range_val, sigma=sigma)

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

n <- 16
n.obs <- 1000
domain_upper_limit <- 50
tol <- 1e-4

loc <- sort(runif(n.obs + n, 0, domain_upper_limit))
while (min(diff(loc)) < tol) {
  loc <- sort(runif(n.obs + n, 0, domain_upper_limit))
}

nu.vec <- 0.6

ind <- 3

results_nngp <- compute_kl_divergence_nngp(loc = loc, obs.ind = 1:1000, m.vec=m.vec, 
                                        nu.vec=nu.vec, range=range_val, sigma=sigma)
