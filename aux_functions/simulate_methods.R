## Simulate Toeplitz

simulate_Toeplitz <- function(loc, nu, kappa, sigma, nsim, print=TRUE, only_time = FALSE, samples = 1){
N <- length(loc)
samp_mat <- matrix(nrow = nsim, ncol = N)
    # Assuming loc is sorted
    loc <- loc - loc[1]
    acf <- rSPDE::matern.covariance(h=loc,kappa=kappa,nu=nu,sigma=sigma)
    acf <- as.vector(acf)
    y <- matrix(rnorm(N * nsim), ncol = N, nrow = nsim)
    time_run2_approx <- 0
    for(ii in 1:samples){
        start <- Sys.time()
        sim <- SuperGauss::cholZX(Z = t(y), acf = acf)
        end <- Sys.time()
        time_temp <- as.numeric(end-start, units = "secs")    
        time_run2_approx <- time_run2_approx + time_temp
    }
    time_run2_approx <- time_run2_approx/samples
    if(print){
        print(time_run2_approx)
    }
    if(only_time){
        return(time_run2_approx)
    }
    return(list(sim = sim, time = time_run2_approx))
}

# # Example:
# s <- seq(0,1,by=0.001)
# nu <- 0.6
# kappa <- 10
# sigma <- 2
# sim <- simulate_Toeplitz(loc = s, nu = nu, kappa = kappa, sigma = sigma, nsim = 10000)
# library(rSPDE)
# c.true <- matern.covariance(0.5-s, kappa=kappa, nu=nu, sigma=sigma)
# plot(s, c.true,
#      type = "l", ylab = "C(|s-0.5|)", xlab = "s", ylim = c(0, 5),
#      cex.main = 0.8, cex.axis = 0.8, cex.lab = 0.8
#    )
# lines(s, cov(t(sim$sim))[(length(s)-1)/2+1,], col = 2)



## Simulate rational markov (our method)
library(Rcpp)

sourceCpp("aux_functions/sample_chol_eigen.cpp")

simulate_rat_markov <- function(loc, m, nu, kappa, sigma, nsim, type_rat = "chebfun", print=TRUE, only_time = FALSE, eigen_chol = FALSE, equally_spaced = FALSE, samples = 1){
N <- length(loc)
samp_mat <- matrix(nrow = nsim, ncol = N)
time_run2_approx <- list()
sim_list <- list()
for(i_m in m){
    if(nu < 0.5){
        # r <- rSPDE::matern.rational.precision(loc, order = i_m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp = "spline", equally_spaced = equally_spaced, ordering = "field")
        r <- rSPDE::matern.rational.ldl(loc, order = i_m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp = "spline", equally_spaced = equally_spaced, ordering = "field")
    } else if ( 0.5 < nu ){
        # r <- rSPDE::matern.rational.precision(loc, order = i_m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "chebfun", type_interp = "spline", equally_spaced = equally_spaced, ordering = "field")
        r <- rSPDE::matern.rational.ldl(loc, order = i_m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "chebfun", type_interp = "spline", equally_spaced = equally_spaced, ordering = "field")
    } 
        # r$Q <- (r$Q + t(r$Q))/2
        L <- r$L %*% sqrt(r$D)
        # L <- as(L, "TsparseMatrix")
        time_run2_approx[[as.character(i_m)]] <- 0
        for(ii in 1:samples){
            y <- matrix(rnorm(ncol(r$L) * nsim), ncol = ncol(r$L), nrow = nsim)
            if(!eigen_chol){
                # R <- Matrix::Cholesky(r$Q, perm=FALSE, LDL=FALSE)
                # sim_full <- solve(R, t(y), system = "Lt")
                start = Sys.time()
                sim_full <- solve(L, t(y), system = "Lt")
                sim <- (r$A)%*%sim_full
                end = Sys.time()
            } else{
                start = Sys.time()            
                sim_full <- sampleCholSparse(r$Q, t(y))
                sim <- (r$A)%*%sim_full
                end = Sys.time()                
            }
            sim_list[[as.character(i_m)]] <- as.matrix(sim)
            time_temp = as.numeric(end-start, units = "secs")
            time_run2_approx[[as.character(i_m)]] <- time_run2_approx[[as.character(i_m)]] + time_temp
        } 
        time_run2_approx[[as.character(i_m)]] <- time_run2_approx[[as.character(i_m)]]/samples 
        if(print){
            print(time_run2_approx[[as.character(i_m)]])
        }
    }
    if(only_time){
        return(time_run2_approx)
    }
    return(list(sim = sim_list, time = time_run2_approx))
}

# # Example:
# source("aux_functions//aux_functions_cov.R")
# s <- seq(0,10,length.out = 2000)
# nu <- 0.3
# kappa <- 10
# sigma <- 2
# sim <- simulate_rat_markov(loc = s, m=1:6, nu = nu, kappa = kappa, sigma = sigma, nsim = 10000, equally_spaced = TRUE)
# library(rSPDE)
# c.true <- matern.covariance(s-5, kappa=kappa, nu=nu, sigma=sigma)
# plot(s, c.true,
#      type = "l", ylab = "C(|s-5|)", xlab = "s", ylim = c(0, 5),
#      cex.main = 0.8, cex.axis = 0.8, cex.lab = 0.8
#    )
# lines(s, cov(t(sim$sim[["6"]]))[(length(s)-1)/2+1,], col = 2)



## Simulate PCA

simulate_PCA <- function(loc, m, nu, kappa, sigma, nsim, print=TRUE, only_time = FALSE, samples = 1){
N <- length(loc)
time_run2_approx <- list()
sim_list <- list()

D_loc <- dist2matR(dist(loc))
cov_mat <- rSPDE::matern.covariance(h=D_loc,kappa=kappa,nu=nu,sigma=sigma)
eigen_cov <- eigen(cov_mat)
for(i_m in m){
    K <- eigen_cov$vec[,1:i_m]    
    D <- diag(eigen_cov$val[1:i_m])    
    time_run2_approx[[as.character(i_m)]] <- 0
    for(ii in 1:samples){
        y <- matrix(rnorm(i_m * nsim), ncol = i_m, nrow = nsim)
        start = Sys.time()
        sim <- K%*%sqrt(D)%*%t(y)
        end =Sys.time()
        sim_list[[as.character(i_m)]] <- as.matrix(sim)
        temp_time <- as.numeric(end-start, units = "secs")
        time_run2_approx[[as.character(i_m)]] <- time_run2_approx[[as.character(i_m)]] + temp_time
    }
    time_run2_approx[[as.character(i_m)]] <- time_run2_approx[[as.character(i_m)]]/samples
    if(print){
            print(time_run2_approx[[as.character(i_m)]])
    }
}
    if(only_time){
        return(time_run2_approx)
    }
return(list(sim = sim_list, time = time_run2_approx))
}

# # Example:
# s <- seq(0,1,by=0.001)
# nu <- 0.6
# kappa <- 10
# sigma <- 2
# sim <- simulate_PCA(loc = s, m=c(10,20,40,60), nu = nu, kappa = kappa, sigma = sigma, nsim = 10000)
# library(rSPDE)
# c.true <- matern.covariance(rep(0.5, length(s)), abs(s), kappa=kappa, nu=nu, sigma=sigma)
# plot(s, c.true,
#      type = "l", ylab = "C(|s-0.5|)", xlab = "s", ylim = c(0, 5),
#      cex.main = 0.8, cex.axis = 0.8, cex.lab = 0.8
#    )
# lines(s, cov(t(sim$sim[["60"]]))[(length(s)-1)/2+1,], col = 2)


## Simulate KL

simulate_KL <- function(loc, m, nu, kappa, sigma, nsim, N_KL = 10000, print=TRUE, only_time=FALSE, samples){
N <- length(loc)
time_run2_approx <- list()
sim_list <- list()

large_KL <- seq(min(loc), max(loc), length.out = N_KL)
kl_loc <- c(loc, large_KL)
kl_loc <- unique(kl_loc)

D_loc <- dist2matR(dist(kl_loc))
cov_mat <- rSPDE::matern.covariance(h=D_loc,kappa=kappa,nu=nu,sigma=sigma)
eigen_cov <- eigen(cov_mat)
for(i_m in m){
    K <- eigen_cov$vec[1:N,1:i_m]    
    D <- diag(eigen_cov$val[1:i_m])    
    time_run2_approx[[as.character(i_m)]] <- 0
    for(ii in 1:samples){
        y <- matrix(rnorm(i_m * nsim), ncol = i_m, nrow = nsim)
        start = Sys.time()
        sim <- K%*%sqrt(D)%*%t(y)
        end =Sys.time()
        sim_list[[as.character(i_m)]] <- as.matrix(sim)
        temp_time <-  as.numeric(end-start, units = "secs")
        time_run2_approx[[as.character(i_m)]] <- time_run2_approx[[as.character(i_m)]] + temp_time
    }
    time_run2_approx[[as.character(i_m)]] <- time_run2_approx[[as.character(i_m)]]/samples
    if(print){
            print(time_run2_approx[[as.character(i_m)]])
    }
}
    if(only_time){
        return(time_run2_approx)
    }
return(list(sim = sim_list, time = time_run2_approx))
}

# # Example:
# s <- seq(0,1,by=0.001)
# nu <- 0.6
# kappa <- 10
# sigma <- 2
# sim <- simulate_KL(loc = s, m=c(10,20,40,60), nu = nu, kappa = kappa, sigma = sigma, nsim = 10000)
# library(rSPDE)
# c.true <- matern.covariance(rep(0.5, length(s)), abs(s), kappa=kappa, nu=nu, sigma=sigma)
# plot(s, c.true,
#      type = "l", ylab = "C(|s-0.5|)", xlab = "s", ylim = c(0, 5),
#      cex.main = 0.8, cex.axis = 0.8, cex.lab = 0.8
#    )
# lines(s, cov(t(sim$sim[["60"]]))[(length(s)-1)/2+1,], col = 2)


# simulate NN 

simulate_rat_NN <- function(loc, m, nu, kappa, sigma, nsim, samples, print=TRUE, only_time = FALSE){
N <- length(loc)
time_run2_approx <- list()
sim_list <- list()

D <- dist2matR(dist(loc))
# perm <- comp.reo.fast(N, m = 0, alpha = 0.6)    
for(i_m in m){
        prec_mat <- get.nnQ(loc=loc,kappa=kappa,nu=nu,sigma=sigma, n.nbr = i_m)
        time_run2_approx[[as.character(i_m)]] <- 0
        # R <- Matrix::Cholesky(prec_mat, perm=FALSE, LDL=FALSE)
        R <- chol(prec_mat)
        for(ii in 1:samples){
            y <- matrix(rnorm(ncol(prec_mat) * nsim), ncol = ncol(prec_mat), nrow = nsim)
            # Starting after as the method directly obtains the cholesky
            start = Sys.time()
            sim_full <- solve(R, t(y), system = "Lt")
            end = Sys.time()
            sim_list[[as.character(i_m)]] <- as.matrix(sim_full)
            time_temp = as.numeric(end-start, units = "secs")
            time_run2_approx[[as.character(i_m)]] <- time_run2_approx[[as.character(i_m)]] + time_temp
        }
        time_run2_approx[[as.character(i_m)]] <- time_run2_approx[[as.character(i_m)]]/samples
        if(print){
            print(time_run2_approx[[as.character(i_m)]])
        }
    }
    if(only_time){
        return(time_run2_approx)
    }
    return(list(sim = sim_list, time = time_run2_approx))
}
# # Example:
# s <- seq(0,1,by=0.001)
# nu <- 0.6
# kappa <- 10
# sigma <- 2
# m <- get_m(nu = 0.6, m = 1:6, method = "nngp")
# sim2 <- simulate_rat_NN(loc = s, m = m, nu = nu, kappa = kappa, sigma = sigma, nsim = 10000, samples=1)
# library(rSPDE)
# c.true <- matern.covariance(0.5-s, kappa=kappa, nu=nu, sigma=sigma)
# plot(s, c.true,
#      type = "l", ylab = "C(|s-0.5|)", xlab = "s", ylim = c(0, 5),
#      cex.main = 0.8, cex.axis = 0.8, cex.lab = 0.8
#    )
# lines(s, cov(t(sim2$sim[["22"]]))[(length(s)-1)/2+1,], col = 2)


# Compute times

compare_times_simulation <- function(N, m, range, sigma, nsim, samples){
    timings_alpha01_rat <- list()
    timings_alpha12_rat <- list()
    timings_alpha23_rat <- list()    
    print("Starting rational")
    for(n_loc in N){
        loc <- seq(0, 100, length.out = n_loc)
        nu <- 0.3
        kappa <- sqrt(8*nu)/range
        timings_alpha01_rat[[as.character(n_loc)]] <-  simulate_rat_markov(loc = loc, m = m, nu = nu, kappa = kappa, sigma = sigma, nsim=nsim, print = FALSE, only_time=TRUE, samples = samples)
        
        nu <- 1.2
        kappa <- sqrt(8*nu)/range         
        timings_alpha12_rat[[as.character(n_loc)]] <- simulate_rat_markov(loc = loc, m = m, nu = nu, kappa = kappa, sigma = sigma, nsim = nsim, print = FALSE, only_time=TRUE, samples = samples)

        nu <- 1.8
        kappa <- sqrt(8*nu)/range         
        timings_alpha23_rat[[as.character(n_loc)]] <- simulate_rat_markov(loc = loc, m = m, nu = nu, kappa = kappa, sigma = sigma, nsim = nsim, print = FALSE, only_time=TRUE, samples = samples)
    }
    print("Starting nnGP")
    timings_alpha01_nngp <- list()
    timings_alpha12_nngp <- list()
    timings_alpha23_nngp <- list()    
    for(n_loc in N){
            loc <- seq(0, 100, length.out = n_loc)
            nu <- 0.3
            kappa <- sqrt(8*nu)/range
            m_nngp <- get_m(nu = nu, m = m, method = "nngp", type = "simulation")
            timings_alpha01_nngp[[as.character(n_loc)]] <- simulate_rat_NN(loc = loc, m = m_nngp, nu = nu, kappa = kappa, sigma = sigma, nsim=nsim, print = FALSE, only_time=TRUE, samples = samples)

            nu <- 1.2
            kappa <- sqrt(8*nu)/range        
            m_nngp <- get_m(nu = nu, m= m, method = "nngp", type = "simulation")  
            timings_alpha12_nngp[[as.character(n_loc)]] <- simulate_rat_NN(loc = loc, m = m_nngp, nu = nu, kappa = kappa, sigma = sigma, nsim=nsim, print = FALSE, only_time=TRUE, samples = samples)

            nu <- 1.8
            kappa <- sqrt(8*nu)/range        
            m_nngp <- get_m(nu = nu, m= m, method = "nngp", type = "simulation")  
            timings_alpha23_nngp[[as.character(n_loc)]] <- simulate_rat_NN(loc = loc, m = m_nngp, nu = nu, kappa = kappa, sigma = sigma, nsim=nsim, print = FALSE, only_time=TRUE, samples = samples)            
    }
    print("Starting PCA")
    timings_alpha01_kl <- list()
    timings_alpha12_kl <- list()
    timings_alpha23_kl <- list()    
    for(n_loc in N){
            loc <- seq(0, 100, length.out = n_loc)
            nu <- 0.3
            kappa <- sqrt(8*nu)/range
            m_kl <- get_m(nu = nu, m = m, method = "kl", type = "simulation")
            timings_alpha01_kl[[as.character(n_loc)]] <- simulate_PCA(loc = loc, m = m_kl, nu = nu, kappa = kappa, sigma = sigma, nsim=nsim, print = FALSE, only_time=TRUE, samples = samples)   

            nu <- 1.2
            kappa <- sqrt(8*nu)/range        
            m_kl <- get_m(nu = nu, m = m, method = "kl", type = "simulation")      
            timings_alpha12_kl[[as.character(n_loc)]] <- simulate_PCA(loc = loc, m = m_kl, nu = nu, kappa = kappa, sigma = sigma, nsim=nsim, print = FALSE, only_time=TRUE, samples = samples)    

            nu <- 1.8
            kappa <- sqrt(8*nu)/range        
            m_kl <- get_m(nu = nu, m = m, method = "kl", type = "simulation")      
            timings_alpha23_kl[[as.character(n_loc)]] <- simulate_PCA(loc = loc, m = m_kl, nu = nu, kappa = kappa, sigma = sigma, nsim=nsim, print = FALSE, only_time=TRUE, samples = samples)    
    }    
    res <- list(rational01 = timings_alpha01_rat, rational12 = timings_alpha12_rat, rational23 = timings_alpha23_rat,
     nngp01 = timings_alpha01_nngp, nngp12 = timings_alpha12_nngp, 
     nngp23 = timings_alpha23_nngp, kl01 = timings_alpha01_kl, kl12 = timings_alpha12_kl, kl23 = timings_alpha23_kl)

    res_df <- data.frame(method = "rational", alpha = c(rep("01", length(m)*length(N)),rep("12", length(m)*length(N)), rep("23", length(m)*length(N))), N = rep(N, each = length(m)), Time = c(unlist(res[["rational01"]]), unlist(res[["rational12"]]), unlist(res[["rational23"]])), m = rep(m, length(N)))
    res_df_tmp <- data.frame(method = "nnGP", alpha = c(rep("01", length(m)*length(N)),rep("12", length(m)*length(N)),rep("23", length(m)*length(N))), N = rep(N, each = length(m)), Time = c(unlist(res[["nngp01"]]), unlist(res[["nngp12"]]), unlist(res[["nngp23"]])), m = rep(m, length(N)))
    res_df <- rbind(res_df, res_df_tmp)
    res_df_tmp <- data.frame(method = "PCA", alpha = c(rep("01", length(m)*length(N)),rep("12", length(m)*length(N)),rep("23", length(m)*length(N))), N = rep(N, each = length(m)), Time = c(unlist(res[["kl01"]]), unlist(res[["kl12"]]), unlist(res[["kl23"]])), m = rep(m, length(N)))
    res_df <- rbind(res_df, res_df_tmp)
    return(res_df)
}
