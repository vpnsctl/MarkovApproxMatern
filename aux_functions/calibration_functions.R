library(rSPDE)
library(Matrix)
source("aux_functions/aux_functions_cov.R")

calibration_PCA <- function(N, n_obs, m_min, m_max, m_step, factor = 100, nu, range, sigma, sigma_e, m_rat, method_pca = "standard", samples, 
type = c("prediction", "sampling"), include_build_precision = TRUE, plot=FALSE){
    type <- type[[1]]
    if(!(type%in%c("prediction", "sampling"))){
        stop("type must be either 'prediction' or 'sampling'.")
    }
    if(m_min < 2 || m_max < 2){
        stop("m_min and m_max need to be greater of equal to 2")
    }
    if(m_min > N || m_max > N){
        stop("m_min and m_max need to be less than N")
    }
    if(n_obs > N){
        stop("n_obs needs to be less or equal to N")
    }
    if(m_min > m_max){
        stop("m_min needs to be less or equal to m_max")
    }
    kappa <- sqrt(8*nu)/range
    sigma.e <- sigma_e
    loc <- seq(0,N/100,length.out=N) 
    n <- N
    n.obs <- n_obs
    D_loc <- as.matrix(dist(loc)) 
    obs_ind <- sort(sample(1:n)[1:n.obs]) 
    cov_mat <- rSPDE::matern.covariance(h=D_loc,kappa=kappa,nu=nu,sigma=sigma) 
    R <- chol(cov_mat[obs_ind,obs_ind]) 
    X <- t(R)%*%rnorm(n.obs) 
    y <- X + sigma_e*rnorm(n.obs) 
    eigen_cov <- eigen(cov_mat) 
    times.rat <- rep(0,length(m_rat))
    for(i in 1:length(m_rat)){
        m.rat <- m_rat[i]
        for(jj in 1:samples){ 
            if(type == "prediction"){
                t1 <- Sys.time() 
                Qrat <- rSPDE:::matern.rational.precision(loc = loc, order = m.rat, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp = "spline")
                if(!include_build_precision){
                    t1 <- Sys.time()
                }
                Q <- Qrat$Q 
                Qhat.rat <- Q + t(Qrat$A[obs_ind,])%*%Qrat$A[obs_ind,]/sigma.e^2 
                mu.rat1 <- Qrat$A%*%solve(Qhat.rat, t(Qrat$A[obs_ind,])%*%y/sigma.e^2) 
                t2 <- Sys.time() 
                times.rat[i] <- times.rat[i] + as.numeric(t2-t1, units = "secs")
            } else if(type == "sampling"){
                t1 <- Sys.time()
                r <- rSPDE:::matern.rational.ldl(loc, order = m.rat, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp = "spline", ordering = "field")
                prec_mat <- t(r$L)%*%r$D%*%r$L
                if(!include_build_precision){
                    t1 <- Sys.time()
                }
                R <- chol(prec_mat)
                y <- matrix(rnorm(ncol(r$L)), ncol = ncol(r$L), nrow = 1)
                sim_full <- solve(R, t(y), system = "Lt")                
                sim <- (r$A)%*%sim_full
                t2 <- Sys.time()
                times.rat[i] <- times.rat[i] + as.numeric(t2-t1, units = "secs")
            } else{
                stop("invalid type")
            }
        }
        times.rat[i] <- times.rat[i]/samples
     } 
     m.pca <- seq(from=m_min,to=m_max, by=m_step) 
     times.pca <- rep(0, length(m.pca)) 

     for(i in 1:length(m.pca)){ 
        i_m <- m.pca[i] 
        for(jj in 1:samples){
            if(type == "prediction"){
                t1 <- Sys.time()                 
                B <- eigen_cov$vec[,1:i_m] 
                Bo <- B[obs_ind,] 
                if(!include_build_precision){
                    t1 <- Sys.time()
                }
                Sigma.w <- Diagonal(i_m,eigen_cov$val[1:i_m]) 
                Q.hat <- solve(Sigma.w) + t(Bo)%*%Bo/sigma_e^2 
                post_mean2 <- B%*%solve(Q.hat, t(Bo)%*%y/sigma_e^2) 
                t2 <- Sys.time() 
                times.pca[i] <- times.pca[i] + as.numeric(t2-t1, units="secs")
            } else if(type == "sampling"){
                t1 <- Sys.time()
                K <- eigen_cov$vec[,1:i_m]   
                if(!include_build_precision){
                    t1 <- Sys.time()
                } 
                D <- Diagonal(i_m, eigen_cov$val[1:i_m])  
                y <- matrix(rnorm(i_m), ncol = i_m, nrow = 1)
                sim <- K%*%sqrt(D)%*%t(y)
                t2 <- Sys.time()
                times.pca[i] <- times.pca[i] + as.numeric(t2 - t1, units="secs")
            } else{
                stop("invalid type")
            }
        }
        times.pca[i] <- times.pca[i]/samples
        } 
    m.cal <- rep(0,length(m_rat)) 
    for(i in 1:length(m_rat)){ 
        m.cal[i] = m.pca[which.min(abs(times.rat[i]-times.pca))] 
    } 

    if(plot){
        plot(m.pca,times.pca, xlab = "m", ylab = "time", ylim = c(0,max(times.pca)), type="l") 
        for(ii in 1:length(m_rat)){
            lines(c(m_min,m_max),c(times.rat[ii],times.rat[ii]),col=ii+1) 
        }
    }
    return(m.cal)
}

# ## Example:

# N <- 1000
# n_obs <- 800
# m_min <- 2
# m_max <- 700
# m_step <- 10
# nu <- 1.4
# range <- 1
# sigma <- 1
# sigma_e <- 0.1
# m_rat <- 2:6


# m_cal_pca_pred <- calibration_PCA(N=N, n_obs=n_obs, m_min=m_min, m_max=m_max, 
#     m_step=m_step, nu=nu, range=range, sigma=sigma, 
#     sigma_e=sigma_e, m_rat=m_rat, samples = 5, type = "prediction", plot=TRUE)

# m_cal_pca_samp <- calibration_PCA(N=N, n_obs=n_obs, m_min=m_min, m_max=m_max, m_step=m_step, nu=nu, range=range, sigma=sigma, 
#     sigma_e=sigma_e, m_rat=m_rat, samples = 10, type = "sampling", plot=TRUE, include_build_precision = FALSE)



#### Note that for nngp we always include the construction of the precision matrix in the timings.
### est_nu -> when type == estimation, should nu be estimated?
### only_optim -> should only the optimization be included in the timing, or should we include the prediction? If FALSE only optmization will be included.

calibration_nngp <- function(N, n_obs, m_min, m_max, m_step, factor = 100, nu, range, sigma, sigma_e, m_rat, samples, 
type = c("prediction", "sampling", "estimation"), est_nu, only_optim, plot=FALSE, print=FALSE){
    type <- type[[1]]
    if(!(type%in%c("prediction", "sampling", "estimation"))){
        stop("type must be either 'prediction', 'sampling' or 'estimation'.")
    }
    if(m_min < 2 || m_max < 2){
        stop("m_min and m_max need to be greater of equal to 2")
    }
    if(m_min > N || m_max > N){
        stop("m_min and m_max need to be less than N")
    }
    if(n_obs > N){
        stop("n_obs needs to be less or equal to N")
    }
    if(m_min > m_max){
        stop("m_min needs to be less or equal to m_max")
    }
    kappa <- sqrt(8*nu)/range
    sigma.e <- sigma_e
    loc <- seq(0,N/100,length.out=N) 
    n <- N
    n.obs <- n_obs
    D_loc <- as.matrix(dist(loc)) 
    obs_ind <- obs.ind <- sort(sample(1:n)[1:n.obs]) 
    cov_mat <- rSPDE::matern.covariance(h=D_loc,kappa=kappa,nu=nu,sigma=sigma) 
    R <- chol(cov_mat[obs_ind,obs_ind]) 
    X <- t(R)%*%rnorm(n.obs) 
    y <- X + sigma_e*rnorm(n.obs) 
    eigen_cov <- eigen(cov_mat) 
    times.rat <- rep(0,length(m_rat))
    if(print){
        print("starting rational")
    }
    for(i in 1:length(m_rat)){
        m.rat <- m_rat[i]
        if(print){
            print(paste("m =",m.rat))
        }
        for(jj in 1:samples){ 
            if(type == "prediction"){
                t1 <- Sys.time() 
                Qrat <- rSPDE:::matern.rational.precision(loc = loc, order = m.rat, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp = "spline")
                Q <- Qrat$Q 
                Qhat.rat <- Q + t(Qrat$A[obs_ind,])%*%Qrat$A[obs_ind,]/sigma.e^2 
                mu.rat1 <- Qrat$A%*%solve(Qhat.rat, t(Qrat$A[obs_ind,])%*%y/sigma.e^2) 
                t2 <- Sys.time() 
                times.rat[i] <- times.rat[i] + as.numeric(t2-t1, units = "secs")
            } else if(type == "sampling"){
                t1 <- Sys.time()
                r <- rSPDE:::matern.rational.ldl(loc, order = m.rat, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp = "spline", ordering = "field")
                prec_mat <- t(r$L)%*%r$D%*%r$L
                R <- chol(prec_mat)
                y <- matrix(rnorm(ncol(r$L)), ncol = ncol(r$L), nrow = 1)
                sim_full <- solve(R, t(y), system = "Lt")                
                sim <- (r$A)%*%sim_full
                t2 <- Sys.time()
                times.rat[i] <- times.rat[i] + as.numeric(t2-t1, units = "secs")
            } else if(type == "estimation"){
                t1 <- Sys.time()
                if(est_nu){
                    theta0 <- c(log(kappa), log(sigma), log(nu), log(sigma_e))
                    par <- optim(theta0, rat.like, method = "L-BFGS-B", loc = loc[obs_ind], Y = y, m = m.rat)
                    t2 <- Sys.time()
                    times.rat[i] <- times.rat[i] + as.numeric((t2-t1)/par$counts[1], units = "secs") 
                    if(!only_optim){
                        kappa_est = exp(par$par[1])
                        sigma_est = exp(par$par[2])
                        nu_est = exp(par$par[3])
                        sigma_e_est = exp(par$par[4])    
                        Qrat <- rSPDE:::matern.rational.precision(loc = loc, order = m.rat, nu = nu_est, kappa = kappa_est, cumsum = TRUE, ordering = "location", sigma = sigma_est, type_rational = "brasil", type_interp = "spline")
                        Q <- Qrat$Q
                        Qhat.rat <- Q + t(Qrat$A[obs_ind,])%*%Qrat$A[obs_ind,]/sigma_e_est^2        
                        R <- Matrix::Cholesky(Qhat.rat, perm = FALSE)         
                        mu.rat <- Qrat$A%*%solve(R, t(Qrat$A[obs_ind,])%*%y/sigma_e_est^2, system = "A")
                        t3 <- Sys.time()
                        times.rat[i] <- times.rat[i] + as.numeric((t3-t2), units = "secs")   
                    }
                } else{
                    t1 <- Sys.time()
                    theta0 <- c(log(kappa), log(sigma), log(sigma_e))
                    par <- optim(theta0, rat.like, method = "L-BFGS-B", loc = loc[obs_ind], Y = y, m = m.rat, nu = nu)
                    t2 <- Sys.time()    
                    times.rat[i] <- times.rat[i] + as.numeric((t2-t1)/par$counts[1], units = "secs") 
                    if(!only_optim){
                        kappa_est = exp(par$par[1])
                        sigma_est = exp(par$par[2])
                        sigma_e_est = exp(par$par[3]) 
                        Qrat <- rSPDE:::matern.rational.precision(loc = loc, order = m.rat, nu = nu, kappa = kappa_est, cumsum = TRUE, ordering = "location", sigma = sigma_est, type_rational = "brasil", type_interp = "spline")
                        Q <- Qrat$Q
                        Qhat.rat <- Q + t(Qrat$A[obs_ind,])%*%Qrat$A[obs_ind,]/sigma_e_est^2        
                        R <- Matrix::Cholesky(Qhat.rat, perm = FALSE)   
                        mu.rat <- Qrat$A%*%solve(R, t(Qrat$A[obs_ind,])%*%y/sigma_e_est^2, system = "A")
                        t3 <- Sys.time()
                        times.rat[i] <- times.rat[i] + as.numeric((t3-t2), units = "secs")         
                    }                   
                }
                
            } else{
                stop("invalid type")
            }
        }
        times.rat[i] <- times.rat[i]/samples
                if(print){
                    print(paste("time =", times.rat[i]))
                }       
     } 
     m.nngp <- seq(from=m_min,to=m_max, by=m_step) 
     times.nngp <- rep(0, length(m.nngp)) 

    if(print){
        print("starting nngp")
    }

     for(i in 1:length(m.nngp)){ 
        i_m <- m.nngp[i] 
        if(print){
            print(paste("m =", i_m))
        }
        for(jj in 1:samples){
            if(type == "prediction"){
                t1 <- Sys.time()                 
                Qnn <- get.nnQ(loc = loc[obs.ind],kappa = kappa,nu = nu,sigma = sigma, n.nbr = i_m)
                Qhat <- Qnn + Diagonal(n.obs)/sigma.e^2        
                mu.nn <- solve(Qhat, y/sigma.e^2)
                Bp <- get.nn.pred(loc = loc, kappa = kappa, nu = nu, sigma = sigma, n.nbr = i_m, S = obs.ind)$B
                mu.nn <- Bp%*%mu.nn
                t2 <- Sys.time() 
                times.nngp[i] <- times.nngp[i] + as.numeric(t2-t1, units="secs")
            } else if(type == "sampling"){
                t1 <- Sys.time()
                prec_mat <- get.nnQ(loc=loc,kappa=kappa,nu=nu,sigma=sigma, n.nbr = i_m)
                R <- chol(prec_mat)
                y <- matrix(rnorm(ncol(prec_mat)), ncol = ncol(prec_mat), nrow = 1)
                sim_full <- solve(R, t(y), system = "Lt")
                t2 <- Sys.time()
                times.nngp[i] <- times.nngp[i] + as.numeric(t2 - t1, units="secs")
            } else if(type == "estimation"){
                t1 <- Sys.time()
                if(est_nu){
                    theta0 <- c(log(kappa), log(sigma), log(nu), log(sigma_e))
                    par <- optim(theta0, nn.like, method = "L-BFGS-B", loc = loc[obs_ind], Y = y, n.nbr = i_m)
                    t2 <- Sys.time()
                    times.nngp[i] <- times.nngp[i] + as.numeric((t2-t1)/par$counts[1], units = "secs") 
                    if(!only_optim){
                        kappa_est = exp(par$par[1])
                        sigma_est = exp(par$par[2])
                        nu_est = exp(par$par[3])
                        sigma_e_est = exp(par$par[4])    
                        Qnn <- get.nnQ(loc = loc[obs_ind],kappa = kappa_est,nu = nu_est,sigma = sigma_est, n.nbr = i_m)
                        A = Diagonal(n.obs)
                        Qhat <- Qnn + t(A)%*%A/sigma_e^2        
                        mu.nn <- solve(Qhat, t(A)%*%y/sigma_e^2)
                        Bp <- get.nn.pred(loc = loc, kappa = kappa_est, nu = nu_est, sigma = sigma_est, n.nbr = i_m, S = obs_ind)$B
                        mu.nn <- Bp%*%mu.nn
                        t3 <- Sys.time()
                        times.nngp[i] <- times.nngp[i] + as.numeric(t3-t2, units="secs")
                    }
                } else{
                    t1 <- Sys.time()
                    theta0 <- c(log(kappa), log(sigma), log(sigma_e))
                    par <- optim(theta0, nn.like,method = "L-BFGS-B", loc = loc[obs_ind], Y = y, n.nbr = i_m, nu = nu)
                    t2 <- Sys.time()    
                    times.nngp[i] <- times.nngp[i] + as.numeric((t2-t1)/par$counts[1], units = "secs") 
                    if(!only_optim){
                        kappa_est = exp(par$par[1])
                        sigma_est = exp(par$par[2])
                        sigma_e_est = exp(par$par[3]) 
                        Qnn <- get.nnQ(loc = loc[obs_ind],kappa = kappa_est,nu = nu,sigma = sigma_est, n.nbr = i_m)
                        A = Diagonal(n.obs)
                        Qhat <- Qnn + t(A)%*%A/sigma_e^2        
                        mu.nn <- solve(Qhat, t(A)%*%y/sigma_e^2)
                        Bp <- get.nn.pred(loc = loc, kappa = kappa_est, nu = nu, sigma = sigma_est, n.nbr = i_m, S = obs_ind)$B
                        mu.nn <- Bp%*%mu.nn
                        t3 <- Sys.time()
                        times.nngp[i] <- times.nngp[i] + as.numeric(t3-t2, units="secs")
                    }                   
                }
            } else{
                stop("invalid type")
            }
        }
        times.nngp[i] <- times.nngp[i]/samples
        if(print){
            print(paste("time =", times.nngp[i]))
        }
        } 
    m.cal <- rep(0,length(m_rat)) 
    for(i in 1:length(m_rat)){ 
        m.cal[i] = m.nngp[which.min(abs(times.rat[i]-times.nngp))] 
    } 

    if(plot){
        plot(m.nngp,times.nngp, xlab = "m", ylab = "time", ylim = c(0,max(times.nngp)), type="l") 
        for(ii in 1:length(m_rat)){
            lines(c(m_min,m_max),c(times.rat[ii],times.rat[ii]),col=ii+1) 
        }
    }
    return(m.cal)
}

# ## Example:

# N <- 1000
# n_obs <- 800
# m_min <- 2
# m_max <- 30
# m_step <- 1
# nu <- 1.4
# range <- 1
# sigma <- 1
# sigma_e <- 0.1
# m_rat <- 2:6


# m_cal_nngp_pred <- calibration_nngp(N=N, n_obs=n_obs, m_min=m_min, m_max=m_max, 
#     m_step=m_step, nu=nu, range=range, sigma=sigma, 
#     sigma_e=sigma_e, m_rat=m_rat, samples = 5, type = "prediction", plot=TRUE)

# m_cal_nngp_samp <- calibration_nngp(N=N, n_obs=n_obs, m_min=m_min, m_max=m_max, m_step=m_step, nu=nu, range=range, sigma=sigma, sigma_e=sigma_e, m_rat=m_rat, samples = 5, type = "sampling", plot=TRUE)


# N <- 100
# n_obs <- 80
# m_min <- 15
# m_max <- 50
# m_step <- 3
# nu <- 1.4
# range <- 1
# sigma <- 1
# sigma_e <- 0.1
# m_rat <- 2:6

# # Estimating nu, only calibrating the optimization

# m_cal_nngp_est_nu_only_optim <- calibration_nngp(N=N, n_obs=n_obs, m_min=m_min, m_max=m_max, m_step=m_step, nu=nu, range=range, sigma=sigma, sigma_e=sigma_e, m_rat=m_rat, samples = 5, type = "estimation", est_nu = TRUE, only_optim = TRUE, plot=TRUE, print=TRUE)

# # Estimating nu, calibrating optimization and prediction

# m_cal_nngp_est_nu_only_optim <- calibration_nngp(N=N, n_obs=n_obs, m_min=m_min, m_max=m_max, m_step=m_step, nu=nu, range=range, sigma=sigma, sigma_e=sigma_e, m_rat=m_rat, samples = 5, type = "estimation", est_nu = TRUE, only_optim = FALSE, plot=TRUE, print=TRUE)

# # Not estimating nu, only calibrating the optimization

# m_cal_nngp_est_only_optim <- calibration_nngp(N=N, n_obs=n_obs, m_min=m_min, m_max=m_max, m_step=m_step, nu=nu, range=range, sigma=sigma, sigma_e=sigma_e, m_rat=m_rat, samples = 5, type = "estimation", est_nu = FALSE, only_optim = TRUE, plot=TRUE, print=TRUE)

# # Not estimating nu, calibrating optimization and prediction

# m_cal_nngp_est_only_optim <- calibration_nngp(N=N, n_obs=n_obs, m_min=m_min, m_max=m_max, m_step=m_step, nu=nu, range=range, sigma=sigma, sigma_e=sigma_e, m_rat=m_rat, samples = 5, type = "estimation", est_nu = FALSE, only_optim = FALSE, plot=TRUE, print=TRUE)
