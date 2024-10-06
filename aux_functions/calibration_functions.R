library(rSPDE)
library(Matrix)
source("aux_functions/aux_functions_cov.R")
source("aux_functions/prediction_computations.R")

# Calibrations obtained PCA:

# N= 5000, nobs = 5000
# alpha01: 268 308 355 406 433 478
# alpha12: 380 473 561 651 708 776
# alpha23: 611  810  945 1082 1205 1325
# N=10000, nobs = 5000
# alpha01: 381 448 533 588 654 727
# alpha12: 532  704  844  953 1065 1162
# alpha23: 904 1202 1431 1622 1808 1965
# N=10000, nobs = 10000
# alpha01: 271 311 367 412 452 493
# alpha12: 372 493 591 672 751 821
# alpha23: 640  848 1016 1168 1299 1420

calcTimePCA <- function(samples, n.obs, type, eigen_cov, include_build_precision, i_m, obs_ind) {
    sigma_e <- 1
    times.pca <- 0 
    for(jj in 1:samples){
        y <- rnorm(n.obs) 
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
            times.pca <- times.pca + as.numeric(t2-t1, units="secs")
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
            times.pca <- times.pca + as.numeric(t2 - t1, units="secs")
        } else{
            stop("invalid type")
        }
    }
    return(times.pca/samples)
}
search_PCA <- function(N, n_obs, m_min, m_max, times, samples, 
                       type = c("prediction", "sampling"), 
                       include_build_precision = TRUE){
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
    n <- N
    n.obs <- n_obs
    obs_ind <- sort(sample(1:n)[1:n.obs]) 
    eigen_cov <- list(vec = matrix(1, ncol=N, nrow=N), val = rep(1, N))
    m.pca <- seq(from=m_min,to=m_max, by=1) 
    times.pca <- rep(0, length(m.pca))
    m.cal <- rep(0, length(times))
    for(i in 1:length(times)){ 
        t.target <- times[i]
        ind.lo <- 1
        ind.up <- length(m.pca)
        m.lo <- m.pca[ind.lo]
        t.lo <- times.pca[ind.lo]
        m.up <- m.pca[ind.up]
        t.up <- times.pca[ind.up]
        found <- FALSE
        while(!found){
            ind.mid <- floor(ind.lo + (ind.up - ind.lo)/2)
            if(ind.mid == ind.lo){
                found <- TRUE
            }
            m.mid <- m.pca[ind.mid]
            t.mid <- times.pca[ind.mid]
            if(t.mid == 0){
                t.mid <- calcTimePCA(samples, n.obs, type, eigen_cov, include_build_precision, m.mid, obs_ind)    
                times.pca[ind.mid] = t.mid
            }
            if(print){
                cat("prog = ", i/length(times),", t.target = ", t.target, ", m.mid = ", m.mid, ", t.mid = ",t.mid, "\n")
            }
            if(t.mid > t.target){
                ind.up <- ind.mid
                m.up <- m.mid
                t.up <- t.mid
            }
            if(t.mid < t.target){
                ind.lo <- ind.mid
                m.lo <- ind.mid
                t.lo <- ind.mid
            }
            
        }
        m.cal[i] <- m.pca[ind.up]
    }
    return(list(m.cal = m.cal, m.pca = m.pca, times.pca = times.pca))
}


timing_PCA <- function(N, n_obs, m_min, m_max, m_step, max.time = Inf, samples, 
                       type = c("prediction", "sampling"), 
                       include_build_precision = TRUE){
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
    n <- N
    n.obs <- n_obs
    obs_ind <- sort(sample(1:n)[1:n.obs]) 
    eigen_cov <- list(vec = matrix(1, ncol=N, nrow=N), val = rep(1, N))
    m.pca <- seq(from=m_min,to=m_max, by=m_step) 
    times.pca <- rep(0, length(m.pca)) 
    cat("compute PCA times\n")
    for(i in 1:length(m.pca)){ 
        i_m <- m.pca[i] 
        for(jj in 1:samples){
            y <- rnorm(n.obs) 
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
        cat(i/length(m.pca), " ", times.pca[i]/max.time,"\n")
        if(times.pca[i]> max.time) {
            cat("Early stop\n")
            break
        }
    }
    return(list(m = m.pca,time = times.pca))
}


calibration <- function(m.rat, times.rat, m.other, times.other, 
                            multi_time = 1, 
                            plot=FALSE){
    
    m.cal <- rep(0,length(m.rat)) 
    for(i in 1:length(m.rat)){ 
        m.cal[i] = m.other[which.min(abs(multi_time*times.rat[i]-times.other))] 
    } 
    
    if(plot){
        plot(m.pca,times.pca, xlab = "m", ylab = "time", ylim = c(0,max(times.pca)), type="l") 
        for(ii in 1:length(m.rat)){
            lines(c(min(m.other),max(m.other)),c(multi_time*times.rat[ii],multi_time*times.rat[ii]),col=ii+1) 
        }
    }
    return(m.cal)
}


calibration_PCA <- function(N, n_obs, m_min, m_max, m_step, factor = 100, nu, range, sigma, 
                            sigma_e, m_rat, method_pca = "standard", samples, 
                            type = c("prediction", "sampling"), 
                            include_build_precision = TRUE, 
                            multi_time = 1, 
                            plot=FALSE, use_chol = TRUE){
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
    obs_ind <- sort(sample(1:n)[1:n.obs]) 
    y <- rnorm(n.obs) 
    times.rat <- rep(0,length(m_rat))
    eigen_cov <- list(vec = matrix(1, ncol=N, nrow=N), val = rep(1, N))
    for(i in 1:length(m_rat)){
        m.rat <- m_rat[i]
        cat("compute rational times, m = ", i, "\n")
        for(jj in 1:samples){ 
            cat(jj/samples, " ")
            if(type == "prediction"){
                t1 <- Sys.time() 
                Qrat <- rSPDE:::matern.rational.precision(loc = loc, order = m.rat, nu = nu, kappa = kappa, cumsum = TRUE, ordering = "location", sigma = sigma, type_rational = "brasil", type_interp = "spline")
                Q <- Qrat$Q 
                Qhat.rat <- Q + t(Qrat$A[obs_ind,])%*%Qrat$A[obs_ind,]/sigma_e^2        
                R <- Matrix::Cholesky(Qhat.rat, perm = FALSE)         
                mu.rat <- Qrat$A%*%solve(R, t(Qrat$A[obs_ind,])%*%y/sigma_e^2, system = "A")
                t2 <- Sys.time() 
                times.rat[i] <- times.rat[i] + as.numeric(t2-t1, units = "secs")
            } else if(type == "sampling"){
                t1 <- Sys.time()
                r <- rSPDE:::matern.rational.ldl(loc, order = m.rat, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp = "spline", ordering = "field")
                if(use_chol){
                    prec_mat <- t(r$L)%*%r$D%*%r$L
                } else{
                    R <- t(r$L)%*%sqrt(r$D)
                }

                if(!include_build_precision){
                    t1 <- Sys.time()
                }
                if(use_chol){
                    R <- chol(prec_mat)
                }
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
     cat("compute PCA times\n")
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
        cat(i/length(m.pca), " ", times.pca[i]/max(times.rat),"\n")
        if(times.pca[i]> 1.1*max(times.rat)) {
            cat("Early stop\n")
            break
        }
        } 
    m.cal <- rep(0,length(m_rat)) 
    for(i in 1:length(m_rat)){ 
        m.cal[i] = m.pca[which.min(abs(multi_time*times.rat[i]-times.pca))] 
    } 

    if(plot){
        plot(m.pca,times.pca, xlab = "m", ylab = "time", ylim = c(0,max(times.pca)), type="l") 
        for(ii in 1:length(m_rat)){
            lines(c(m_min,m_max),c(multi_time*times.rat[ii],multi_time*times.rat[ii]),col=ii+1) 
        }
    }
    return(m.cal)
}

# ## Example:

# N <- 1000
# n_obs <- 800
# m_min <- 2
# m_max <- 1000
# m_step <- 10
# nu <- 1.4
# range <- 1
# sigma <- 1
# sigma_e <- 0.1
# m_rat <- 1:6


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
    obs_ind <- obs.ind <- sort(sample(1:n)[1:n.obs]) 
    y <- rnorm(n.obs) 
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
                Qrat <- rSPDE:::matern.rational.precision(loc = loc, order = m.rat, nu = nu, kappa = kappa, cumsum = TRUE, ordering = "location", sigma = sigma, type_rational = "brasil", type_interp = "spline")
                Q <- Qrat$Q 
                Qhat.rat <- Q + t(Qrat$A[obs_ind,])%*%Qrat$A[obs_ind,]/sigma_e^2        
                R <- Matrix::Cholesky(Qhat.rat, perm = FALSE)         
                mu.rat <- Qrat$A%*%solve(R, t(Qrat$A[obs_ind,])%*%y/sigma_e^2, system = "A")
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


# Timing NNGP


timing_nngp <- function(N, n_obs, m_min, m_max, m_step, factor = 100, nu, range, sigma, sigma_e, samples, 
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
    loc <- seq(0,N/factor,length.out=N) 
    n <- N
    n.obs <- n_obs
    obs_ind <- obs.ind <- sort(sample(1:n)[1:n.obs]) 
    y <- rnorm(n.obs) 

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
            cat(jj/samples, " ")
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
    return(times.nngp)
}












timing_rat <- function(N, n_obs, nu, range, sigma, sigma_e, m_rat, samples, 
                            type = c("prediction", "sampling", "estimation"), 
                            include_build_precision = TRUE, 
                            multi_time = 1, 
                            plot=FALSE, use_chol = TRUE, est_nu = FALSE, only_optim = TRUE, factor = 100, print=FALSE){
    type <- type[[1]]
    if(!(type%in%c("prediction", "sampling", "estimation"))){
        stop("type must be either 'prediction', 'sampling' or 'estimation'.")
    }
    if(n_obs > N){
        stop("n_obs needs to be less or equal to N")
    }
    
    kappa <- sqrt(8*nu)/range
    sigma.e <- sigma_e
    loc <- seq(0,N/factor,length.out=N) 
    n <- N
    n.obs <- n_obs
    obs_ind <- sort(sample(1:n)[1:n.obs]) 
    y <- rnorm(n.obs) 
    times.rat <- rep(0,length(m_rat))
    eigen_cov <- list(vec = matrix(1, ncol=N, nrow=N), val = rep(1, N))
    for(i in 1:length(m_rat)){
        m.rat <- m_rat[i]
        if(print){
            cat("compute rational times, m = ", m.rat, "\n")
        }
        for(jj in 1:samples){ 
            if(print){
                cat(jj/samples, " ")
            }
            if(type == "prediction"){
                t1 <- Sys.time() 
                Qrat <- rSPDE:::matern.rational.precision(loc = loc, order = m.rat, nu = nu, kappa = kappa, cumsum = FALSE, ordering = "location", sigma = sigma, type_rational = "brasil", type_interp = "spline")
                Q <- Qrat$Q 
                Qhat.rat <- Q + t(Qrat$A[obs_ind,])%*%Qrat$A[obs_ind,]/sigma_e^2        
                R <- Matrix::Cholesky(Qhat.rat, perm = FALSE)         
                mu.rat <- Qrat$A%*%solve(R, t(Qrat$A[obs_ind,])%*%y/sigma_e^2, system = "A")
                t2 <- Sys.time() 
                times.rat[i] <- times.rat[i] + as.numeric(t2-t1, units = "secs")
            } else if(type == "sampling"){
                t1 <- Sys.time()
                r <- rSPDE:::matern.rational.ldl(loc, order = m.rat, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp = "spline", ordering = "field")
                if(use_chol){
                    prec_mat <- t(r$L)%*%r$D%*%r$L
                } else{
                    R <- t(r$L)%*%sqrt(r$D)
                }
                
                if(!include_build_precision){
                    t1 <- Sys.time()
                }
                if(use_chol){
                    R <- chol(prec_mat)
                }
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
        if(print){
            cat("\n")
        }
        times.rat[i] <- times.rat[i]/samples
        print(times.rat[i])
    } 
    return(times.rat)
}


timing_nngp_m_vec <- function(N, n_obs, m.nngp, factor = 100, nu, range, sigma, sigma_e, samples, 
type = c("prediction", "sampling", "estimation"), est_nu, only_optim, plot=FALSE, print=FALSE){
    type <- type[[1]]
    if(!(type%in%c("prediction", "sampling", "estimation"))){
        stop("type must be either 'prediction', 'sampling' or 'estimation'.")
    }
    if(n_obs > N){
        stop("n_obs needs to be less or equal to N")
    }

    kappa <- sqrt(8*nu)/range
    sigma.e <- sigma_e
    loc <- seq(0,N/factor,length.out=N) 
    n <- N
    n.obs <- n_obs
    obs_ind <- obs.ind <- sort(sample(1:n)[1:n.obs]) 
    y <- rnorm(n.obs) 

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
            if(print){
                cat(jj/samples, " ")
            }
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
    return(times.nngp)
}


timing_taper_m_vec <- function(N, n_obs, m.taper, factor = 100, nu, range, sigma, sigma_e, samples, print=FALSE){
        if(n_obs > N){
            stop("n_obs needs to be less or equal to N")
        }
    
        kappa <- sqrt(8*nu)/range
        sigma.e <- sigma_e
        loc <- seq(0,N/factor,length.out=N) 
        n <- N
        n.obs <- n_obs
        obs_ind <- obs.ind <- sort(sample(1:n)[1:n.obs]) 
        y <- rnorm(n.obs) 
    
        times.taper <- rep(0, length(m.taper))     
    
        if(print){
            print("starting taper")
        }
    
         for(i in 1:length(m.taper)){ 
            i_m <- m.taper[i] 
            if(print){
                print(paste("m =", i_m))
            }
            for(jj in 1:samples){
                if(print){
                    cat(jj/samples, " ")
                }
                    t1 <- Sys.time()                 
                    taper_pred(y = y, loc_full = loc, idx_obs = obs_ind, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e, m = i_m)
                    t2 <- Sys.time() 
                    times.taper[i] <- times.taper[i] + as.numeric(t2-t1, units="secs")
            }
            times.taper[i] <- times.taper[i]/samples
            if(print){
                print(paste("time =", times.taper[i]))
            }
            } 
        return(times.taper)
    }




auto_calibration_nngp_rat <- function(n, n_obs, nu, range, sigma, sigma_e, samples, m_rat, previous_calibration = NULL, max_it_per_m = 20, print=FALSE){
    ret_m <- numeric(length(m_rat))
    if(is.null(previous_calibration)){
        # get rational times
        if(print){
            print("Computing rational times")
        }
        times_rat <- timing_rat(N = n, n_obs = n_obs, nu = nu, range = range, sigma = sigma, sigma_e = sigma_e, m_rat = m_rat, samples = samples, type = "prediction")
        if(print){
            print("Calibrating nngp")
        }
        m_tmp <- 1
        times_tmp <- timing_nngp_m_vec(N = n, n_obs = n_obs, m.nngp = m_tmp, nu = nu, range = range, sigma = sigma, sigma_e = sigma_e, samples = samples, type = "prediction")
        for(i in 1:length(m_rat)){
            if(print){
                print(paste("Starting calibration for m =",m_rat[i]))
            }
            count <- 0
            if(times_tmp >= times_rat[i]){
                ret_m[i] <- m_tmp
            } else{
                while((times_tmp < times_rat[i]) && (count < max_it_per_m)){
                    m_tmp <- m_tmp + 1
                    times_tmp <- timing_nngp_m_vec(N = n, n_obs = n_obs, m.nngp = m_tmp, nu = nu, range = range, sigma = sigma, sigma_e = sigma_e, samples = samples, type = "prediction")
                    count <- count + 1
                }
                ret_m[i] <- m_tmp
            }
            if(print){
                print("Calibration found:")
                print(paste("m_nngp =", ret_m[i]))
            }
        }
    } else{
        if(print){
            print("Computing rational times")
        }
        times_rat <- timing_rat(N = n, n_obs = n_obs, nu = nu, range = range, sigma = sigma, sigma_e = sigma_e, m_rat = m_rat, samples = samples, type = "prediction")
        if(print){
            print("Calibrating nngp")
        }
        for(i in 1:length(m_rat)){
            if(print){
                print(paste("Starting calibration for m =",m_rat[i]))
            }
            count <- 0
            m_tmp <- previous_calibration[i]
            times_tmp <- timing_nngp_m_vec(N = n, n_obs = n_obs, m.nngp = m_tmp, nu = nu, range = range, sigma = sigma, sigma_e = sigma_e, samples = samples, type = "prediction")

            if(times_tmp > times_rat[i]){
                while((times_tmp > times_rat[i]) &&  (count < max_it_per_m) && (m_tmp > 1)){
                    m_tmp <- m_tmp - 1
                    times_tmp <- timing_nngp_m_vec(N = n, n_obs = n_obs, m.nngp = m_tmp, nu = nu, range = range, sigma = sigma, sigma_e = sigma_e, samples = samples, type = "prediction")
                    count <- count + 1
                }
                if(times_tmp < times_rat[i]){
                    m_tmp <- m_tmp + 1
                }
            } else{
                while((times_tmp < times_rat[i]) &&  (count < max_it_per_m)){
                    m_tmp <- m_tmp + 1
                    times_tmp <- timing_nngp_m_vec(N = n, n_obs = n_obs, m.nngp = m_tmp, nu = nu, range = range, sigma = sigma, sigma_e = sigma_e, samples = samples, type = "prediction")
                    count <- count + 1
                }
            }

            ret_m[i] <- m_tmp

            if(print){
                print("Calibration found:")
                print(paste("m_nngp =", ret_m[i]))
            }
        }
    }
    return(ret_m)
}

# time1 <- Sys.time()
# m_tmp <- auto_calibration_nngp_rat(n=250, n_obs=250, nu=1.4, range=0.5, sigma=1, sigma_e=0.1, samples=50, m_rat=1:6, previous_calibration = NULL, max_it_per_m = 20, print=TRUE) 
# time2 <- Sys.time()
# print(time2-time1)

# time1 <- Sys.time()
# m_tmp2 <- auto_calibration_nngp_rat(n=500, n_obs=250, nu=1.4, range=0.5, sigma=1, sigma_e=0.1, samples=50, m_rat=1:6, previous_calibration = m_tmp, max_it_per_m = 20, print=TRUE)
# time2 <- Sys.time()
# print(time2-time1)




auto_calibration_taper_rat <- function(n, n_obs, nu, range, sigma, sigma_e, samples, m_rat, step_size = 1, previous_calibration = NULL, max_it_per_m = 20, print=FALSE){
    ret_m <- numeric(length(m_rat))
    if(is.null(previous_calibration)){
        # get rational times
        if(print){
            print("Computing rational times")
        }
        times_rat <- timing_rat(N = n, n_obs = n_obs, nu = nu, range = range, sigma = sigma, sigma_e = sigma_e, m_rat = m_rat, samples = samples, type = "prediction", print = print)
        if(print){
            print("Calibrating taper")
        }
        m_tmp <- 1
        times_tmp <- timing_taper_m_vec(N = n, n_obs = n_obs, m.taper = m_tmp, nu = nu, range = range, sigma = sigma, sigma_e = sigma_e, samples = samples, print = print)
        for(i in 1:length(m_rat)){
            if(print){
                print(paste("Starting calibration for m =",m_rat[i]))
            }
            count <- 0
            if(times_tmp >= times_rat[i]){
                ret_m[i] <- m_tmp
            } else{
                while((times_tmp < times_rat[i]) && (count < max_it_per_m)){
                    m_tmp <- m_tmp + step_size
                    times_tmp <- timing_taper_m_vec(N = n, n_obs = n_obs, m.taper = m_tmp, nu = nu, range = range, sigma = sigma, sigma_e = sigma_e, samples = samples, print = print)
                    count <- count + 1
                }
                ret_m[i] <- m_tmp
            }
            if(print){
                print("Calibration found:")
                print(paste("m_taper =", ret_m[i]))
            }
        }
    } else{
        if(print){
            print("Computing rational times")
        }
        times_rat <- timing_rat(N = n, n_obs = n_obs, nu = nu, range = range, sigma = sigma, sigma_e = sigma_e, m_rat = m_rat, samples = samples, type = "prediction", print = print)
        if(print){
            print("Calibrating taper")
        }
        for(i in 1:length(m_rat)){
            if(print){
                print(paste("Starting calibration for m =",m_rat[i]))
            }
            count <- 0
            m_tmp <- previous_calibration[i]
            times_tmp <- timing_taper_m_vec(N = n, n_obs = n_obs, m.taper = m_tmp, nu = nu, range = range, sigma = sigma, sigma_e = sigma_e, samples = samples, print = print)

            if(times_tmp > times_rat[i]){
                while((times_tmp > times_rat[i]) &&  (count < max_it_per_m) && (m_tmp > 1)){
                    m_tmp <- m_tmp - step_size
                    times_tmp <- timing_taper_m_vec(N = n, n_obs = n_obs, m.taper = m_tmp, nu = nu, range = range, sigma = sigma, sigma_e = sigma_e, samples = samples, print = print)
                    count <- count + 1
                }
                if(times_tmp < times_rat[i]){
                    m_tmp <- m_tmp + step_size
                }
            } else{
                while((times_tmp < times_rat[i]) &&  (count < max_it_per_m)){
                    m_tmp <- m_tmp + step_size
                    times_tmp <- timing_taper_m_vec(N = n, n_obs = n_obs, m.taper = m_tmp, nu = nu, range = range, sigma = sigma, sigma_e = sigma_e, samples = samples, print = print)
                    count <- count + 1
                }
            }

            ret_m[i] <- m_tmp

            if(print){
                print("Calibration found:")
                print(paste("m_taper =", ret_m[i]))
            }
        }
    }
    return(ret_m)
}

# time1 <- Sys.time()
m_tmp <- auto_calibration_taper_rat(n=5000, n_obs=5000, nu=2.4, range=0.5, sigma=1, sigma_e=0.1, samples=200, m_rat=1:6, step_size = 10, previous_calibration = NULL, max_it_per_m = 20, print=TRUE) 
# m_tmp <- c(1, 1, 1, 41, 61, 81)
# m_tmp <- auto_calibration_taper_rat(n=5000, n_obs=5000, nu=2.4, range=0.5, sigma=1, sigma_e=0.1, samples=200, step_size = 10, m_rat=1:6, previous_calibration = m_tmp, max_it_per_m = 20, print=TRUE) 

# Calibration: Taper, nu = 1.4 : 1   1   1  46 111 156

# time2 <- Sys.time()
# print(time2-time1)

# time1 <- Sys.time()
# m_tmp2 <- auto_calibration_taper_rat(n=500, n_obs=250, nu=1.4, range=0.5, sigma=1, sigma_e=0.1, samples=50, m_rat=1:6, previous_calibration = m_tmp, max_it_per_m = 20, print=TRUE)
# time2 <- Sys.time()
# print(time2-time1)
