library(rSPDE)
source("aux_functions/aux_functions_cov.R")

compute_likelihood_rat <- function(y, loc, obs_ind, m_vec, kappa, sigma, nu, sigma_e, est_nu = TRUE, method = "L-BFGS-B"){
    post_mean <- list()
    timings <- list()
    for(m in m_vec){
        timings[[as.character(m)]] <- list()
        if(est_nu){
            t1 <- Sys.time()
            theta0 <- c(log(kappa), log(sigma), log(nu), log(sigma_e))
            par <- optim(theta0, rat.like, method = "L-BFGS-B", loc = loc[obs_ind], Y = y, m = m)
            t2 <- Sys.time()
            kappa_est = exp(par$par[1])
            sigma_est = exp(par$par[2])
            nu_est = exp(par$par[3])
            sigma_e_est = exp(par$par[4])    
            Qrat <- matern.rational.ldl(loc = loc, order = m, nu = nu_est, kappa = kappa_est, 
                                    sigma = sigma_est, type_rational = "brasil", type_interp =  "spline")   
            t3 <- Sys.time()                                           
        } else{
            t1 <- Sys.time()
            theta0 <- c(log(kappa), log(sigma), log(sigma_e))
            par <- optim(theta0, rat.like, method = "L-BFGS-B", loc = loc[obs_ind], Y = y, m = m, nu = nu)
            t2 <- Sys.time()    
            kappa_est = exp(par$par[1])
            sigma_est = exp(par$par[2])
            sigma_e_est = exp(par$par[3])          
            Qrat <- matern.rational.ldl(loc = loc, order = m, nu = nu, kappa = kappa_est, 
                                sigma = sigma_est, type_rational = "brasil", type_interp =  "spline")
            t3 <- Sys.time()          
        }

        Q <- t(Qrat$L)%*%Qrat$D%*%Qrat$L

        Qhat.rat <- Q + t(Qrat$A[obs_ind,])%*%Qrat$A[obs_ind,]/sigma_e_est^2        
        mu.rat <- Qrat$A%*%solve(Qhat.rat, t(Qrat$A[obs_ind,])%*%y/sigma_e_est^2)
        t4 <- Sys.time()
        timings[[as.character(m)]][["Optimization"]] <- as.numeric(t2-t1, units = "secs")
        timings[[as.character(m)]][["Precision_computation"]] <- as.numeric(t3-t2, units = "secs")
        timings[[as.character(m)]][["Posterior_mean"]] <- as.numeric(t4-t3, units = "secs")
        timings[[as.character(m)]][["Optimization_counts"]] <- par$counts[1]
        post_mean[[as.character(m)]] <- mu.rat
    }
    return(list(timings = timings, post_mean = post_mean))
}


compute_likelihood_nn <- function(y, loc, obs_ind, m_vec, kappa, sigma, nu, sigma_e, est_nu = TRUE, method = "L-BFGS-B"){
    post_mean <- list()
    timings <- list()
    for(m in m_vec){
        timings[[as.character(m)]] <- list()
        if(est_nu){
            t1 <- Sys.time()
            theta0 <- c(log(kappa), log(sigma), log(nu), log(sigma_e))
            par <- optim(theta0, nn.like, method = "L-BFGS-B", loc = loc[obs_ind], Y = y, n.nbr = m)
            t2 <- Sys.time()
            kappa_est = exp(par$par[1])
            sigma_est = exp(par$par[2])
            nu_est = exp(par$par[3])
            sigma_e_est = exp(par$par[4])    
            Qnn <- get.nnQ(loc = loc,kappa = kappa_est,nu = nu_est,sigma = sigma_est, n.nbr = m, S = obs_ind)
            t3 <- Sys.time()                                           
        } else{
            t1 <- Sys.time()
            theta0 <- c(log(kappa), log(sigma), log(sigma_e))
            par <- optim(theta0, nn.like, method = "L-BFGS-B", loc = loc[obs_ind], Y = y, n.nbr = m, nu = nu)
            t2 <- Sys.time()    
            kappa_est = exp(par$par[1])
            sigma_est = exp(par$par[2])
            sigma_e_est = exp(par$par[3])          
            Qnn <- get.nnQ(loc = loc,kappa = kappa_est,nu = nu,sigma = sigma_est, n.nbr = m, S = obs_ind) 
            t3 <- Sys.time()          
        }

        Qnn <- get.nnQ(loc = loc[obs_ind],kappa = kappa_est,nu = nu,sigma = sigma_est, n.nbr = m)

        A = Diagonal(n.obs)
        Qhat <- Qnn + t(A)%*%A/sigma_e^2        
        mu.nn <- solve(Qhat, t(A)%*%y/sigma_e^2)
        Bp <- get.nn.pred(loc = loc, kappa = kappa_est, nu = nu, sigma = sigma_est, n.nbr = m, S = obs_ind)$B
        mu.nn <- Bp%*%mu.nn
        t4 <- Sys.time()
        timings[[as.character(m)]][["Optimization"]] <- as.numeric(t2-t1, units = "secs")
        timings[[as.character(m)]][["Precision_computation"]] <- as.numeric(t3-t2, units = "secs")
        timings[[as.character(m)]][["Posterior_mean"]] <- as.numeric(t4-t3, units = "secs")
        timings[[as.character(m)]][["Optimization_counts"]] <- par$counts[1]        
        post_mean[[as.character(m)]] <- mu.nn
    }
    return(list(timings = timings, post_mean = post_mean))
}


compute_timings_likelihood <- function(loc, obs_ind, nu_vec, range, sigma, sigma_e, m, m_nngp_fun, est_nu){
    res <- list()
    res[["rat"]] <- list()
    res[["nngp"]] <- list()
    res[["true"]] <- list()
    res[["Y"]] <- list()
    for(nu in nu_vec){
      alpha <- nu + 1/2
      kappa = sqrt(8*nu)/range  
      Y <- sample_y(loc[obs_ind], nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e)
    
      res[["rat"]][[as.character(nu)]] <- compute_likelihood_rat(y=Y, loc=loc, obs_ind=obs_ind, m_vec=m, 
                          kappa=kappa, sigma=sigma, nu = nu, sigma_e = sigma_e, est_nu = est_nu)

      m_nngp <- m_nngp_fun(m, alpha)
      res[["nngp"]][[as.character(nu)]] <- compute_likelihood_nn(y=Y, loc=loc, obs_ind=obs_ind, m_vec=m_nngp, 
                          kappa=kappa, sigma=sigma, nu = nu,sigma_e = sigma_e, est_nu = est_nu)
      res[["true"]][[as.character(nu)]] <-  true_pred(Y, loc = loc[obs_ind], loc_pred = loc, nu = nu, 
      kappa = kappa, sigma = sigma, sigma_e = sigma_e)  
      res[["Y"]][[as.character(nu)]] <- Y 
    }
    return(res)
}
