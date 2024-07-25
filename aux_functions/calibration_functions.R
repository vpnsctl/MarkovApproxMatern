library(rSPDE)
library(Matrix)

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
