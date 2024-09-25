rm(list=ls())

m_nngp_fun <- function(m, alpha, n, n.obs){
    if(alpha<1) {
        mn <- m - 1
        if(mn < 1){
            mn <- 1
        }
    } else if (alpha < 2) {
        if(n == 1000) {
            m_vec = c(1, 2, 6, 14, 20, 24)
        } else if(n == 5000){
            m_vec <- c(1, 2, 13, 21, 27, 31)
        } else if(n.obs == 5000){
            m_vec <- c(1, 2, 3, 4, 14, 20)
        } else if(n.obs == 10000){
            m_vec <- c(1, 2, 7, 18, 24, 29)
        } else{
            stop("not implemented")
        }
        mn <- m_vec[m]
    } else {
        if(n == 1000){
            m_vec <- c(4, 17, 24, 29, 32, 36) 
        } else if(n == 5000){
            m_vec <- c(15, 30, 37, 45, 51, 54)
        } else if(n.obs == 5000){
            m_vec <- c(1, 22, 31, 39, 47, 54)
        } else if(n.obs == 10000){
            m_vec <- c(14, 28, 37, 44, 51, 57)
        } else{
            stop("not implemented")
        }
        mn <- m_vec[m]
    }
    return(mn)
}


source("markov_approx/aux_functions/aux_functions_cov.R")
source("markov_approx/aux_functions/probability_computations.R")
library(rSPDE)
library(excursions)
library(mvtnorm)

use.excursions <- TRUE
folder_to_save <- "markov_approx"

range = 1
sigma = 1
sigma.e <- 0.1
n <- 1000 
n.obs <- seq(from=10,to=1000,by=10)
n.rep <- 100
nu <- 1.9 - 1/2
alpha <- nu + 1/2
coverage <- 0.9
m_vec = 2:6   # rational order
mn_vec <- rep(0,length(m_vec))
for(i in 1:length(m_vec)){
    mn_vec[i] = m_nngp_fun(m_vec[i], alpha, n, 1) #nngp order    
}


## Generate data
loc <- seq(0,n/100,length.out=n)
kappa = sqrt(8*nu)/range
true_cov <- get_cov_mat(loc = loc, m = NULL, method = "true", nu = nu, kappa = kappa, sigma = sigma, samples = NULL, L=NULL)

prob_nngp <- prob_true <- prob_rat <- matrix(0, nrow = length(n.obs), ncol = n.rep)
err_nngp <- err_rat <- matrix(0, length(m_vec), length(n.obs))
for(i in 1:length(n.obs)) {
    time0 <- Sys.time()
    for(j in 1:n.rep){
        time1 <- Sys.time()
        cat(i,j, 'generate data\n')
        obs.ind <- sort(sample(1:n)[1:n.obs[i]])
        z <- matrix(rnorm(n.obs[i]), ncol = n.obs[i], nrow = 1)
        L <- chol(true_cov[obs.ind,obs.ind])
        y <- t(L)%*%t(z) + sigma.e * rnorm(n.obs[i])
        
        ## Compute credible band and true probabilitiy
        cat(i,j, 'compute truth\n')
        post_true <- posterior_constructor_cov(true_cov, y, sigma.e, obs.ind, 1:n, type = "all")
        post_mean_true <- as.vector(post_true$post_mean)
        post_cov_true <- as.matrix(post_true$post_cov)
        
        Q.true <- solve(post_cov_true)
        conf <- simconf(alpha = 1-coverage, mu = post_mean_true, Q = Q.true, vars  = diag(post_cov_true))
        lb_prob <- conf$a
        ub_prob <- conf$b
        if(use.excursions) {
            prob_true <- gaussint(a=lb_prob,b=ub_prob,mu=post_mean_true,Q = Q.true, n = 1e6)$P        
        } else {
            prob_true <- pmvnorm(lower=lb_prob,upper=ub_prob,mean=post_mean_true,sigma = post_cov_true)
        }
        
        for(k in 1:length(m_vec)){
            m <- m_vec[k]
            ## Rational approx
            cat(i,j,k, 'rational\n')
            Qrat <-rSPDE:::matern.rational.precision(loc = loc, order = m, nu = nu, kappa = kappa, sigma = sigma, 
                                                     type_rational = "brasil", type_interp =  "spline")    
            A_obs <- Qrat$A[obs.ind,]
            A_mat = t(A_obs)
            Q_xgiveny <-(A_mat%*% (A_obs))/sigma.e^2 + Qrat$Q
            post_y <- (A_mat%*% y)/sigma.e^2
            R <- Matrix::Cholesky(Q_xgiveny, perm = FALSE)         
            mu_xgiveny <- solve(R, post_y, system = "A")
            mu.rat <-  as.vector(Qrat$A %*% mu_xgiveny)
            Sigma.rat <- as.matrix(Qrat$A%*%solve(Q_xgiveny, t(Qrat$A)))
            
            if(use.excursions) {
                prob_rat <- gaussint(a=lb_prob,b=ub_prob,mu=mu.rat,Q = solve(Sigma.rat), n.iter = 1e6)$P
            } else {
                prob_rat <- pmvnorm(lower=lb_prob,upper=ub_prob,mean=mu.rat,sigma = Sigma.rat)
            }
            err_rat[k,i] <- err_rat[k,i] + prob_rat-prob_true
            
            # nngp
            mn <- mn_vec[k]
            cat(i,j,k, 'nngp\n')
            prec_nngp <- get_cov_mat(loc = loc[obs.ind], m = mn, method = "nngp", nu = nu, kappa = kappa, sigma = sigma, samples = NULL, L=NULL)
            post_nngp <- posterior_constructor_nngp(prec_mat = prec_nngp, y=y, sigma_e = sigma.e, 1:n, obs.ind, loc, mn, nu, kappa, sigma, type = "full")
            post_mean_nngp <- as.vector(post_nngp$post_mean)
            post_cov_nngp <- as.matrix(post_nngp$post_cov)
            
            if(use.excursions) {
                prob_nngp <-  gaussint(a=lb_prob,b=ub_prob,mu=post_mean_nngp,Q = solve(post_cov_nngp), n.iter = 1e6)$P
            } else {
                prob_nngp <- pmvnorm(lower=lb_prob,upper=ub_prob,mean=post_mean_nngp,sigma = post_cov_nngp)
            }
            
            err_nngp[k,i] <- err_nngp[k,i] + prob_nngp-prob_true
            
            cat("rat: ", err_rat[k,i]/j, '\n')
            cat("nngp: ", err_nngp[k,i]/j, '\n')
        }
                time2 <- Sys.time()
                print("Time for replicate")
                print(time2-time1)
    }
    time3 <- Sys.time()
    print("Time for observation")
    print(time3-time0)
    res_tmp <- list(err_rat = err_rat,
    err_nngp = err_nngp,
    nu = nu,
    n.obs = n.obs,
    n.rep = n.rep,
    iteration = i)
    dir.create(file.path(folder_to_save, "prob_tables"), showWarnings = FALSE)
    dir.create(file.path(paste0(folder_to_save, "/prob_tables/"), "partials"), showWarnings = FALSE)
    dir.create(file.path(paste0(folder_to_save, "/prob_tables/partials"),paste0("nu_",nu)), showWarnings = FALSE)
    saveRDS(res_tmp, paste0(folder_to_save,"/prob_tables/res_","_range",range,"_nu_",nu,".RDS"))
}
err_rat <- err_rat/n.rep
err_nngp <- err_nngp/n.rep

res <- list(err_rat = err_rat,
    err_nngp = err_nngp,
    nu = nu,
    n.obs = n.obs,
    n.rep = n.rep)

dir.create(file.path(folder_to_save, "prob_tables"), showWarnings = FALSE)
saveRDS(res, paste0(folder_to_save, "/prob_tables/res_","_range",range,"_nu_",nu,".RDS"))