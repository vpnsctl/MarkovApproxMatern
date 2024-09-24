rm(list=ls())
source("aux_functions/calibration_functions.R")

m_nngp_fun <- function(m, alpha = NULL, n, m_vec = NULL){
    return(m_vec[m])

    # if(alpha<1) {
    #     mn <- m - 1
    #     if(mn < 1){
    #         mn <- 1
    #     }
    # } else if (alpha < 2) {
    #     # if(n == 1000) {
    #     m_vec = c(1, 2, 6, 14, 20, 24)
    #     # } else if(n == 5000){
    #     #     m_vec <- c(1, 2, 13, 21, 27, 31)
    #     # } else if(n.obs == 5000){
    #     #     m_vec <- c(1, 2, 3, 4, 14, 20)
    #     # } else if(n.obs == 10000){
    #     #     m_vec <- c(1, 2, 7, 18, 24, 29)
    #     # } else{
    #     #     stop("not implemented")
    #     # }
    #     mn <- m_vec[m]
    # } else {
    #     # if(n == 1000){
    #     m_vec <- c(4, 17, 24, 29, 32, 36) 
    #     # } else if(n == 5000){
    #     #     m_vec <- c(15, 30, 37, 45, 51, 54)
    #     # } else if(n.obs == 5000){
    #     #     m_vec <- c(1, 22, 31, 39, 47, 54)
    #     # } else if(n.obs == 10000){
    #     #     m_vec <- c(14, 28, 37, 44, 51, 57)
    #     # } else{
    #     #     stop("not implemented")
    #     # }
    #     mn <- m_vec[m]
    # }
    # return(mn)
}


source("aux_functions/aux_functions_cov.R")
source("aux_functions/probability_computations.R")
library(rSPDE)
library(excursions)
library(mvtnorm)

use.excursions <- TRUE
folder_to_save <- "markov_approx"

range = 1
sigma = 1
sigma.e <- 0.1
n <- c(0, 25, 50, 100, 150, 250,  300, 400, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000)
n.obs <- 250
n.rep <- 10
max_it_per_m <- 20 # for calibration
samples_calibration <- 50
nu <- 1.9 - 1/2
alpha <- nu + 1/2
coverage <- 0.9
m_vec = 1:6   # rational order
mn_vec <- rep(0,length(m_vec))

prob_nngp <- prob_true <- prob_rat <- matrix(0, nrow = length(n), ncol = n.rep)
err_nngp <- err_rat <- matrix(0, length(m_vec), length(n))

locs <- list()
for(i in 1:length(n)) {
    locs[[i]] = seq(from=1,to=25,length.out = n[i])
}

kappa = sqrt(8*nu)/range

previous_calibration <- NULL

for(j in 1:n.rep){
    obs_loc <- sort(25*runif(n.obs))
    obs_cov <- get_cov_mat(loc = obs_loc, m = NULL, method = "true", nu = nu, kappa = kappa, sigma = sigma, 
                            samples = NULL, L=NULL)
    z <- matrix(rnorm(n.obs), ncol = n.obs, nrow = 1)
    L <- chol(obs_cov)
    y <- t(L)%*%t(z) + sigma.e * rnorm(n.obs)
    
    
    for(i in 1:length(n)) {
        # add mesh and observation locations and sort them
        if(n[i]>0){ 
            loc <- c(obs_loc, locs[[i]])
            obs.ind <- 1:length(loc)
            tmp <- sort(loc, index.return = TRUE)
            obs.ind[tmp$ix] <- 1:length(loc)
            
            loc <- tmp$x
            obs.ind <- obs.ind[1:n.obs]    
        } else {
            loc <- obs_loc
            obs.ind <- 1:n.obs
        }

        print("Calibrating...")
        # Get calibration
        previous_calibration <- auto_calibration_nngp_rat(n=n[i]+n.obs, n_obs=n.obs, nu=nu, range=range, sigma=sigma, sigma_e=sigma.e, 
                    samples=samples_calibration, m_rat=m_vec, previous_calibration = previous_calibration, max_it_per_m = max_it_per_m, print=FALSE)
        
        print(paste("Calibration:", previous_calibration))
        
        
        ## Compute credible band and true probabilitiy
        cat("rep = ", j,"n =", n[i], 'compute truth\n')
        
        true_cov <- get_cov_mat(loc = loc, m = NULL, method = "true", nu = nu, kappa = kappa, sigma = sigma, 
                                samples = NULL, L=NULL)
        post_true <- posterior_constructor_cov(true_cov, y, sigma.e, obs.ind, 1:(n[i]), type = "all")
        post_mean_true <- as.vector(post_true$post_mean)
        post_cov_true <- as.matrix(post_true$post_cov)
        
        Q.true <- solve(post_cov_true)
        conf <- simconf(alpha = 1-coverage, mu = post_mean_true, Q = Q.true, vars  = diag(post_cov_true), n.iter=1e5)
        lb_prob <- conf$a
        ub_prob <- conf$b
        if(use.excursions) {
            prob_true <- gaussint(a=lb_prob,b=ub_prob,mu=post_mean_true,Q = Q.true, n.iter = 1e5)$P        
        } else {
            prob_true <- pmvnorm(lower=lb_prob,upper=ub_prob,mean=post_mean_true,sigma = post_cov_true)
        }
        
        for(k in 1:length(m_vec)){
            m <- m_vec[k]
            ## Rational approx
            cat("rep = ", j,"n =", n[i], 'rational, m = ', m, '\n')
            
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
                prob_rat <- gaussint(a=lb_prob,b=ub_prob,mu=mu.rat,Q = solve(Sigma.rat), n.iter = 1e5)$P
            } else {
                prob_rat <- pmvnorm(lower=lb_prob,upper=ub_prob,mean=mu.rat,sigma = Sigma.rat)
            }
            err_rat[k,i] <- err_rat[k,i] + prob_rat-prob_true
            
            cat("prob rat: ", prob_rat, '\n')
            
            # nngp
            
            mn <- previous_calibration[k]
            # mn <- m_nngp_fun(m, alpha, n[i])
            cat("rep = ", j,"n =", n[i], 'nngp, m = ', mn, '\n')
            
            prec_nngp <- get_cov_mat(loc = loc[obs.ind], m = mn, method = "nngp", 
                                     nu = nu, kappa = kappa, sigma = sigma, samples = NULL, L=NULL)
            post_nngp <- posterior_constructor_nngp(prec_mat = prec_nngp, y=y, sigma_e = sigma.e, 
                                                    1:n[i], obs.ind, loc, mn, nu, kappa, sigma, type = "full")
            post_mean_nngp <- as.vector(post_nngp$post_mean)
            post_cov_nngp <- as.matrix(post_nngp$post_cov)
            
            if(use.excursions) {
                prob_nngp <-  gaussint(a=lb_prob,b=ub_prob,mu=post_mean_nngp,Q = solve(post_cov_nngp), n.iter = 1e5)$P
            } else {
                prob_nngp <- pmvnorm(lower=lb_prob,upper=ub_prob,mean=post_mean_nngp,sigma = post_cov_nngp)
            }
            
            cat("prob nngp: ", prob_nngp, '\n')
            
            err_nngp[k,i] <- err_nngp[k,i] + prob_nngp-prob_true
            
            cat("rat: ", err_rat[k,i]/j, '\n')
            cat("nngp: ", err_nngp[k,i]/j, '\n')
        }
    }
}
err_rat <- err_rat/n.rep
err_nngp <- err_nngp/n.rep

res_list <- list(err_rat = err_rat, err_nngp = err_nngp, nu = nu, range = range)

saveRDS(res_list, paste0("prob_errors_range",range,"_nu",nu,".RDS"))


plot(n,abs(err_rat[1,]), type="l", ylim = c(0.00001,1), log="y", ylab = "Absolute error")
lines(n,abs(err_rat[2,]),col=2)
lines(n,abs(err_rat[3,]),col=3)
lines(n,abs(err_rat[4,]),col=4)
lines(n,abs(err_rat[5,]),col=5)
lines(n,abs(err_nngp[1,]),col=1, lty=2)
lines(n,abs(err_nngp[2,]),col=2, lty=2)
lines(n,abs(err_nngp[3,]),col=3, lty=2)
lines(n,abs(err_nngp[4,]),col=4, lty=2)
lines(n,abs(err_nngp[5,]),col=5, lty=2)