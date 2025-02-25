# Set up basic parallel processing using base R
num_cores <- parallel::detectCores() - 1
cl <- parallel::makeCluster(num_cores)

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
range_val <- 0.5
sigma <- 1

library(rSPDE)
library(Matrix)

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

# Export necessary objects to all cores
parallel::clusterExport(cl, c("m.vec", "nu.vec", "range_val", "sigma", "n.obs", "compute_kl_divergence_rational", "build_matern_rational_cov", "apply_matern_by_row", "matern.covariance", "determinant", "dist", "matern.rational.cov"))

# Parallel computation
results_list <- parallel::parLapply(cl, locs, function(loc) {
  result <- compute_kl_divergence_rational(
    loc = loc,
    m.vec = m.vec,
    nu.vec = nu.vec,
    range = range_val,
    sigma = sigma
  )
  
  # Return with metadata
  list(
    KL = result$KL,
    type = attr(result, "type"),
    nu.vec = attr(result, "nu.vec"),
    m.vec = attr(result, "m.vec"),
    n = length(loc) - n.obs
  )
})

# Stop the cluster
parallel::stopCluster(cl)

# Organize results into a data frame
results_df <- do.call(rbind, lapply(seq_along(results_list), function(i) {
  data.frame(
    n = results_list[[i]]$n,
    KL_matrix = I(list(results_list[[i]]$KL)),
    type = results_list[[i]]$type,
    nu_vec = I(list(results_list[[i]]$nu.vec)),
    m_vec = I(list(results_list[[i]]$m.vec))
  )
}))



## NN GP - nearest neighbor

library(pracma)
library(rSPDE)

# Set up basic parallel processing using base R
num_cores <- parallel::detectCores() - 1
cl <- parallel::makeCluster(num_cores)

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
range_val <- 0.5
sigma <- 1

get.neighbor <- function(i, n.nbr, S = NULL) {
    if(is.null(S)) {
        S = 1:i
        S.i <- sort(S[S<i], decreasing = TRUE)
        if(length(S.i) < n.nbr) {
            return(sort(S.i))
        } else {
            return(sort(S.i[1:n.nbr]))
        }    
    } else {
        if(i %in% S) { 
            S.i <- sort(S[S<i], decreasing = TRUE)
            if(length(S.i) < n.nbr) {
                return(sort(S.i))
            } else {
                return(sort(S.i[1:n.nbr]))
            }    
        } else {
             dists <- abs(i - S)
             return(sort(sort(dists, index.return = TRUE)$ix[1:n.nbr]))
        }
    }
    
    
}

get.nn <- function(loc,kappa,nu,sigma, n.nbr, S = NULL) {
    k <- length(loc)
    ii <- numeric((2*k-n.nbr)*(n.nbr+1)/2)
    jj <- numeric((2*k-n.nbr)*(n.nbr+1)/2)
    val <- numeric((2*k-n.nbr)*(n.nbr+1)/2)
    
    Fs.d <- Fsi.d <- numeric(k)
    
    Fs.d[1] <- sigma^2
    Fsi.d[1] <- 1/sigma^2
    val[1] <- 1
    ii[1] <- 1
    jj[1] <- 1
    counter <- 1
    for(i in 2:k) {
        nbrs <- get.neighbor(i,n.nbr, S)
        if(isempty(nbrs)) {
            Fs.d[i] <- sigma^2
            Fsi.d[i] <- 1/sigma^2
            val[counter +1] <- 1
            ii[counter+1] <- i
            jj[counter+1] <- i
            counter = counter + 1
        } else {
            Sigma.in <- matern.covariance(h = abs(loc[nbrs]-loc[i]),
                                          kappa = kappa, nu = nu, sigma = sigma)
            Sigma.nn <- matern.covariance(h = as.matrix(dist(loc[nbrs])),
                                          kappa = kappa, nu = nu, sigma = sigma)
            N <- length(nbrs) + 1
            tmp <- solve(Sigma.nn, Sigma.in)
            val[counter + (1:N)] <- c(-t(tmp),1)
            ii[counter + (1:N)] <- rep(i,N)
            jj[counter + (1:N)] <- c(nbrs,i)
            Fs.d[i] <- sigma^2 - t(Sigma.in)%*%tmp
            Fsi.d[i] <- 1/Fs.d[i]    
            counter <- counter + N
        }
        
        
        
        
    }
    ii <- ii[1:counter]
    jj <- jj[1:counter]
    val <- val[1:counter]
    Bs <-  Matrix::sparseMatrix(i   = ii,
                                j    = jj,
                                x    = val,
                                dims = c(k, k))
    Fs <-  Matrix::Diagonal(k,Fs.d)
    Fsi <-  Matrix::Diagonal(k,Fsi.d)
    return(list(Bs = Bs, Fs = Fs, Fi = Fsi))
}

get.nnQ <- function(loc,kappa,nu,sigma, n.nbr,S=NULL) {
    tmp <- get.nn(loc = loc,kappa = kappa,nu = nu,sigma = sigma,n.nbr = n.nbr,S = S)
    
    return(t(tmp$Bs) %*% tmp$Fi %*%tmp$Bs)
}

get.nnCov <- function(loc,kappa,nu,sigma, n.nbr,S=NULL, invert_prec = TRUE) {
    if(!invert_prec){
        tmp <- get.nn(loc = loc,kappa = kappa,nu = nu,sigma = sigma,n.nbr = n.nbr,S = S)
        Bsi <- solve(tmp$Bs)
        return(Bsi %*% tmp$Fs %*%t(Bsi))
    } else{
        Q <- get.nnQ(loc = loc,kappa = kappa,nu = nu,sigma = sigma,n.nbr = n.nbr,S = S)
        return(solve(Q))
    }
}

get.nn.pred <- function(loc,kappa,nu,sigma, n.nbr, S = NULL) {
    
    n.S <- length(S)
    n.loc <- length(loc)
    k <- n.loc - n.S
    N <- k*n.nbr + n.S
    ii <- numeric(N)
    jj <- numeric(N)
    val <- numeric(N)
    Fs.d <- numeric(n.loc)
    counter <- 0
    ii[1:n.S] <- S
    jj[1:n.S] <- 1:n.S
    val[1:n.S] <- 1
    counter <- n.S
    not.observed <- setdiff(1:length(loc), S)
    reo <- c(S, not.observed)
    for(i in not.observed) {
            dists <- abs(loc[i] - loc[S])
            nbrs <- sort(sort(dists, index.return = TRUE)$ix[1:n.nbr])
            Sigma.in <- matern.covariance(h = abs(loc[S[nbrs]]-loc[i]),
                                          kappa = kappa, nu = nu, sigma = sigma)
            Sigma.nn <- matern.covariance(h = as.matrix(dist(loc[S[nbrs]])),
                                          kappa = kappa, nu = nu, sigma = sigma)
            
            tmp <- solve(Sigma.nn, Sigma.in)
            val[counter + (1:n.nbr)] <- t(tmp)
            ii[counter + (1:n.nbr)] <- rep(i,n.nbr)
            jj[counter + (1:n.nbr)] <- nbrs
            counter <- counter + n.nbr    
            Fs.d[i] <- sigma^2 - t(Sigma.in)%*%tmp
    }
    Bs <-  Matrix::sparseMatrix(i   = ii,
                                j    = jj,
                                x    = val,
                                dims = c(n.loc, n.S))
    Fs <-  Matrix::Diagonal(n.loc,Fs.d)
    return(list(B = Bs, F = Fs))
}

m_nngp <- readRDS("calibration_results//calibrated_m_list_regular_mesh.RDS")

m_nngp_fun <- function(m, N){
    return(m_nngp[[as.character(N)]][[m]])
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
            N_m <- length(loc) - n.obs
            i_m <- m_nngp_fun(m, N_m)
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

results_list <- lapply(locs, function(loc) {
  result <- compute_kl_divergence_nngp(
    loc = loc,
    obs.ind = 1:n.obs,
    m.vec = m.vec,
    nu.vec = nu.vec,
    range = range_val,
    sigma = sigma
  )
  
  # Return with metadata
  list(
    KL = result$KL,
    type = attr(result, "type"),
    nu.vec = attr(result, "nu.vec"),
    m.vec = attr(result, "m.vec"),
    n = length(loc) - n.obs
  )
})

# Organize results into a data frame
results_df <- do.call(rbind, lapply(seq_along(results_list), function(i) {
  data.frame(
    n = results_list[[i]]$n,
    KL_matrix = I(list(results_list[[i]]$KL)),
    type = results_list[[i]]$type,
    nu_vec = I(list(results_list[[i]]$nu.vec)),
    m_vec = I(list(results_list[[i]]$m.vec))
  )
}))


############
### FEM ####
############

source("aux_functions/aux_functions_cov.R")

# Set up basic parallel processing using base R
num_cores <- parallel::detectCores() - 1
cl <- parallel::makeCluster(num_cores)

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

m_fem <- readRDS("calibration_results//calibrated_m_list_fem_regular_mesh.RDS")

m_fem_fun <- function(m, N){
    return(m_fem[[as.character(N)]][[m]])
}

library(rSPDE)
library(Matrix)

compute_kl_divergence_fem <- function(loc, m.vec, nu.vec, range, sigma) {
    print("Computing KL divergence for FEM approximation")
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
            
            m_fem <- m_fem_fun(m, length(loc) - n.obs)

            Sigma_fem <- fem_cov(m=m, mesh_fem = m_fem, loc = loc, nu = nu, kappa = kappa, sigma = sigma)  
            
            tryCatch({
                # Compute log determinant of rational approximation
                Sigma_fem_inv <- solve(Sigma_fem)
                log_det_fem <- determinant(Sigma_fem, logarithm=TRUE)$modulus

                temp <- Sigma_fem_inv %*% Sigma.t
                trace_term <- sum(diag(temp))
                dim <- nrow(Sigma.t)
                kl.fem <- as.double(0.5 * (trace_term - dim  + log_det_fem - log_det_true))

                kl.div[i,j] <- kl.fem
                # Print partial KL components
                cat("\nKL components for nu =", nu, "m =", m, ":\n")
                cat("KL divergence:", kl.fem, "\n\n")                
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
    
    return(ret)
}

# Export necessary objects to all cores
parallel::clusterExport(cl, c("m.vec", "nu.vec", "range_val", "sigma", "n.obs", "compute_kl_divergence_fem", "fem_cov", "m_fem_fun", "m_fem", "matern.operators"))

# Parallel computation
results_list <- lapply(locs, function(loc) {
  result <- compute_kl_divergence_fem(
    loc = loc,
    m.vec = m.vec,
    nu.vec = nu.vec,
    range = range_val,
    sigma = sigma
  )
  
  # Return with metadata
  list(
    KL = result$KL,
    type = attr(result, "type"),
    nu.vec = attr(result, "nu.vec"),
    m.vec = attr(result, "m.vec"),
    n = length(loc) - n.obs
  )
})

# Organize results into a data frame
results_df <- do.call(rbind, lapply(seq_along(results_list), function(i) {
  data.frame(
    n = results_list[[i]]$n,
    KL_matrix = I(list(results_list[[i]]$KL)),
    type = results_list[[i]]$type,
    nu_vec = I(list(results_list[[i]]$nu.vec)),
    m_vec = I(list(results_list[[i]]$m.vec))
  )
}))
