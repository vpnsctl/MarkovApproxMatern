## Predict Toeplitz

predict_Toeplitz <- function(loc, nu, kappa, sigma, sigma_e, samples, print=TRUE){
N <- length(loc)
time_run2_sample <- numeric(samples)
    # Assuming loc is sorted
    loc <- loc - loc[1]
    acf = rSPDE::matern.covariance(h=loc,kappa=kappa,nu=nu,sigma=sigma)
    acf <- as.vector(acf)
    acf[1] = acf[1]+sigma_e^2
    acf2 =rSPDE::matern.covariance(h=loc,kappa=kappa,nu=nu,sigma=sigma)
    acf2 =as.vector(acf2)
    for (j in (1:samples)){
        Tz <- SuperGauss::Toeplitz$new(acf = acf)
        Tz2 <- SuperGauss::Toeplitz$new(acf = acf2)
        y=rnorm(N)
        start = Sys.time()
        d = Tz$solve(y)
        post_mean = Tz2$prod(d)
        end =Sys.time()
        time_run2_sample[j] = as.numeric(end-start, units = "secs")
        }
        time_run2_approx = sum(time_run2_sample[j])/samples        
        if(print){
            print(time_run2_approx)
        }
    return(time_run2_approx)
}

## Predict rational markov (our method)

predict_rat_markov <- function(loc, m, nu, kappa, sigma, sigma_e, samples, print=TRUE, equally_spaced = FALSE){
N <- length(loc)
time_run2_sample <- list()
time_run2_approx <- list()

for(i_m in m){
    if(nu < 0.5){
        r <- rSPDE::matern.rational.precision(loc, order = i_m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp = "spline", equally_spaced = equally_spaced, ordering = "location")
    } else if ( 0.5 < nu && nu < 1.5){
        r <- rSPDE::matern.rational.precision(loc, order = i_m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "chebfun", type_interp = "spline", equally_spaced = equally_spaced, ordering = "location")
    } else{
        stop("nu must be between 0 and 1.5")
    }
        r$Q <- (r$Q + t(r$Q))/2
        A_mat = t(r$A)
        Q_xgiveny <-(A_mat%*% (r$A))/sigma_e^2 + r$Q
        for (j in (1:samples))
        {
        y <- rnorm(N)
            start = Sys.time()
            post_y <- (A_mat%*% y)/sigma_e^2
            R <- Matrix::Cholesky(Q_xgiveny, perm=FALSE)
            mu_xgiveny <- solve(R, post_y, system = "A")
            approx_mean1 <-  r$A %*% mu_xgiveny
            end = Sys.time()
        time_run2_sample[[as.character(i_m)]][j] = as.numeric(end-start, units = "secs")
        }
        time_run2_approx[[as.character(i_m)]] = sum(time_run2_sample[[as.character(i_m)]])/samples        
        if(print){
            print(time_run2_approx[[as.character(i_m)]])
        }
    }
    return(time_run2_approx)
}

## Predict PCA

predict_PCA <- function(loc, m, nu, kappa, sigma, sigma_e, samples, print=TRUE){
N <- length(loc)
time_run2_sample <- list()
time_run2_approx <- list()

for(i in m){
  time_run2_sample[[as.character(i)]]<-rep(0, samples)
}    
D_loc <- dist2matR(dist(loc))
cov_mat <- rSPDE::matern.covariance(h=D_loc,kappa=kappa,nu=nu,sigma=sigma)
eigen_cov <- eigen(cov_mat)
for(i_m in m){
    K <- eigen_cov$vec[,1:i_m]    
    D <- diag(eigen_cov$val[1:i_m])    
    cov_KL <- K%*%D%*%t(K)
    svd_K <- svd(K%*%sqrt(D))
    cov_KL_svd_U <- cov_KL %*% svd_K$u
    for (j in (1:samples))
    {
        y=rnorm(N)
        start = Sys.time()
        y_new <- t(svd_K$u) %*% y
        prec_nugget <- cov_KL_svd_U %*% Matrix::Diagonal(x = 1/(svd_K$d^2 + sigma_e^2)) 
        post_mean = prec_nugget%*%y_new
        end =Sys.time()
        time_run2_sample[[as.character(i_m)]][j] = as.numeric(end-start, units = "secs")
        }
        time_run2_approx[[as.character(i_m)]] = sum(time_run2_sample[[as.character(i_m)]])/samples        
    if(print){
            print(time_run2_approx[[as.character(i_m)]])
    }
}
return(time_run2_approx)
}



## Predict KL

predict_KL <- function(loc, m, nu, kappa, sigma, sigma_e, samples, print=TRUE){
N <- length(loc)
time_run2_sample <- list()
time_run2_approx <- list()

large_KL <- seq(min(loc), max(loc), length.out = N_KL)
kl_loc <- c(loc, large_KL)
kl_loc <- unique(kl_loc)

for(i in m){
  time_run2_sample[[as.character(i)]]<-rep(0, samples)
}    
D_loc <- dist2matR(dist(loc))
cov_mat <- rSPDE::matern.covariance(h=D_loc,kappa=kappa,nu=nu,sigma=sigma)
eigen_cov <- eigen(cov_mat)
for(i_m in m){
    K <- eigen_cov$vec[1:N,1:i_m]    
    D <- diag(eigen_cov$val[1:i_m])    
    cov_KL <- K%*%D%*%t(K)
    svd_K <- svd(K%*%sqrt(D))
    cov_KL_svd_U <- cov_KL %*% svd_K$u    
    for (j in (1:samples))
    {
        y=rnorm(N)
        start = Sys.time()
        y_new <- t(svd_K$u) %*% y
        prec_nugget <- cov_KL_svd_U %*% Matrix::Diagonal(x = 1/(svd_K$d^2 + sigma_e^2)) 
        post_mean = prec_nugget%*%y_new
        end =Sys.time()
        time_run2_sample[[as.character(i_m)]][j] = as.numeric(end-start, units = "secs")
        }
        time_run2_approx[[as.character(i_m)]] = sum(time_run2_sample[[as.character(i_m)]])/samples        
    if(print){
            print(time_run2_approx[[as.character(i_m)]])
    }
}
return(time_run2_approx)
}


# predict NN 

predict_rat_NN <- function(loc, m, nu, kappa, sigma, sigma_e, samples, print=TRUE){
N <- length(loc)
time_run2_sample <- list()
time_run2_approx <- list()

for(i in m){
  time_run2_sample[[as.character(i)]]<-rep(0, samples)
}    
D <- dist2matR(dist(loc))
# perm <- comp.reo.fast(N, m = 0, alpha = 0.6)
Sigma <- rSPDE::matern.covariance(h=D,kappa=kappa,nu=nu,sigma=sigma)        
for(i_m in m){
        prec_mat <- get.nnQ(Sigma, i_m)
        I <- Matrix::Diagonal(n = ncol(prec_mat), x = 1)
        Q_xgiveny <- I * 1/sigma_e^2 + prec_mat
        for (j in (1:samples))
        {
        y <- rnorm(N)
        start = Sys.time()
              post_y <- y/sigma_e^2
              R <- Matrix::Cholesky(Q_xgiveny, perm=FALSE)
              mu_xgiveny <- solve(R, post_y, system = "A")
        end = Sys.time()
        time_run2_sample[[as.character(i_m)]][j] = as.numeric(end-start, units = "secs")
        }
        time_run2_approx[[as.character(i_m)]] = sum(time_run2_sample[[as.character(i_m)]])/samples
        if(print){
            print(time_run2_approx[[as.character(i_m)]])
        }
    }
    return(time_run2_approx)
}



compare_times <- function(N, m, range, sigma, samples){
    sigma_e <- 0.1
    timings_alpha01_rat <- list()
    timings_alpha12_rat <- list()
    print("Starting rational")
    for(n_loc in N){
        loc <- seq(0, 1, length.out = n_loc)
        nu <- 0.3
        kappa <- sqrt(8*nu)/range
        timings_alpha01_rat[[as.character(n_loc)]] <- predict_rat_markov(loc = loc, m = m, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e, samples = samples, print = FALSE)
        nu <- 1.2
        kappa <- sqrt(8*nu)/range            
        timings_alpha12_rat[[as.character(n_loc)]] <- predict_rat_markov(loc = loc, m = m, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e, samples = samples, print = FALSE)
    }
    print("Starting nnGP")
    timings_alpha01_nngp <- list()
    timings_alpha12_nngp <- list()
    for(n_loc in N){
            loc <- seq(0, 1, length.out = n_loc)
            nu <- 0.3
            kappa <- sqrt(8*nu)/range
            m_nngp <- get_m(nu = nu, m = m, method = "nngp")
            timings_alpha01_nngp[[as.character(n_loc)]] <- predict_rat_NN(loc = loc, m = m_nngp, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e, samples = samples, print = FALSE)
            nu <- 1.2
            kappa <- sqrt(8*nu)/range        
            m_nngp <- get_m(nu = nu, m= m, method = "nngp")    
            timings_alpha12_nngp[[as.character(n_loc)]] <- predict_rat_NN(loc = loc, m = m_nngp, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e, samples = samples, print = FALSE)
    }
    print("Starting PCA")
    timings_alpha01_kl <- list()
    timings_alpha12_kl <- list()
    for(n_loc in N){
            loc <- seq(0, 1, length.out = n_loc)
            nu <- 0.3
            kappa <- sqrt(8*nu)/range
            m_kl <- get_m(nu = nu, m = m, method = "kl")
            timings_alpha01_kl[[as.character(n_loc)]] <- predict_PCA(loc = loc, m = m_kl, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e, samples = samples, print = FALSE)
            nu <- 1.2
            kappa <- sqrt(8*nu)/range        
            m_kl <- get_m(nu = nu, m = m, method = "kl")                
            timings_alpha12_kl[[as.character(n_loc)]] <- predict_PCA(loc = loc, m = m_kl, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e, samples = samples, print = FALSE)
    }    
    res <- list(rational01 = timings_alpha01_rat, rational12 = timings_alpha12_rat, nngp01 = timings_alpha01_nngp, nngp12 = timings_alpha12_nngp, kl01 = timings_alpha01_kl, kl12 = timings_alpha12_kl)

    res_df <- data.frame(method = "rational", alpha = c(rep("01", length(m)*length(N)),rep("12", length(m)*length(N))), N = rep(N, each = length(m)), Time = c(unlist(res[["rational01"]]), unlist(res[["rational12"]])), m = rep(m, length(N)))
    res_df_tmp <- data.frame(method = "nnGP", alpha = c(rep("01", length(m)*length(N)),rep("12", length(m)*length(N))), N = rep(N, each = length(m)), Time = c(unlist(res[["nngp01"]]), unlist(res[["nngp12"]])), m = rep(m, length(N)))
    res_df <- rbind(res_df, res_df_tmp)
    res_df_tmp <- data.frame(method = "KL", alpha = c(rep("01", length(m)*length(N)),rep("12", length(m)*length(N))), N = rep(N, each = length(m)), Time = c(unlist(res[["kl01"]]), unlist(res[["kl12"]])), m = rep(m, length(N)))
    res_df <- rbind(res_df, res_df_tmp)
    return(res_df)
}

# only with nngp

compare_times_nngp <- function(N, m, range, sigma, samples){
    sigma_e <- 0.1
    timings_alpha01_rat <- list()
    timings_alpha12_rat <- list()
    print("Starting rational")
    for(n_loc in N){
        loc <- seq(0, 1, length.out = n_loc)
        nu <- 0.3
        kappa <- sqrt(8*nu)/range
        timings_alpha01_rat[[as.character(n_loc)]] <- predict_rat_markov(loc = loc, m = m, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e, samples = samples, print = FALSE)
        nu <- 1.2
        kappa <- sqrt(8*nu)/range            
        timings_alpha12_rat[[as.character(n_loc)]] <- predict_rat_markov(loc = loc, m = m, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e, samples = samples, print = FALSE)
    }
    print("Starting nnGP")
    timings_alpha01_nngp <- list()
    timings_alpha12_nngp <- list()
    for(n_loc in N){
            loc <- seq(0, 1, length.out = n_loc)
            nu <- 0.3
            kappa <- sqrt(8*nu)/range
            m_nngp <- get_m(nu = nu, m = m, method = "nngp")
            timings_alpha01_nngp[[as.character(n_loc)]] <- predict_rat_NN(loc = loc, m = m_nngp, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e, samples = samples, print = FALSE)
            nu <- 1.2
            kappa <- sqrt(8*nu)/range        
            m_nngp <- get_m(nu = nu, m = m, method = "nngp")    
            timings_alpha12_nngp[[as.character(n_loc)]] <- predict_rat_NN(loc = loc, m = m_nngp, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e, samples = samples, print = FALSE)
    }
    res <- list(rational01 = timings_alpha01_rat, rational12 = timings_alpha12_rat, nngp01 = timings_alpha01_nngp, nngp12 = timings_alpha12_nngp)
    res_df <- data.frame(method = "rational", alpha = c(rep("01", length(m)*length(N)),rep("12", length(m)*length(N))), N = rep(N, each = length(m)), Time = c(unlist(res[["rational01"]]), unlist(res[["rational12"]])), m = rep(m, length(N)))
    res_df_tmp <- data.frame(method = "nnGP", alpha = c(rep("01", length(m)*length(N)),rep("12", length(m)*length(N))), N = rep(N, each = length(m)), Time = c(unlist(res[["nngp01"]]), unlist(res[["nngp12"]])), m = rep(m, length(N)))
    res_df <- rbind(res_df, res_df_tmp)
    return(res_df)
}


# only with PCA

compare_times_pca <- function(N, m, range, sigma, samples){
    sigma_e <- 0.1
    timings_alpha01_rat <- list()
    timings_alpha12_rat <- list()
    print("Starting rational")
    for(n_loc in N){
        loc <- seq(0, 1, length.out = n_loc)
        nu <- 0.3
        kappa <- sqrt(8*nu)/range
        timings_alpha01_rat[[as.character(n_loc)]] <- predict_rat_markov(loc = loc, m = m, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e, samples = samples, print = FALSE)
        nu <- 1.2
        kappa <- sqrt(8*nu)/range            
        timings_alpha12_rat[[as.character(n_loc)]] <- predict_rat_markov(loc = loc, m = m, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e, samples = samples, print = FALSE)
    }
    print("Starting PCA")
    timings_alpha01_pca <- list()
    timings_alpha12_pca <- list()
    for(n_loc in N){
            loc <- seq(0, 1, length.out = n_loc)
            nu <- 0.3
            kappa <- sqrt(8*nu)/range
            m_nngp <- get_m(nu = nu, m = m, method = "kl")
            timings_alpha01_pca[[as.character(n_loc)]] <- predict_PCA(loc = loc, m = m_nngp, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e, samples = samples, print = FALSE)
            nu <- 1.2
            kappa <- sqrt(8*nu)/range        
            m_nngp <- get_m(nu = nu, m = m, method = "kl")    
            timings_alpha12_pca[[as.character(n_loc)]] <- predict_PCA(loc = loc, m = m_nngp, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e, samples = samples, print = FALSE)
    }
    res <- list(rational01 = timings_alpha01_rat, rational12 = timings_alpha12_rat, pca01 = timings_alpha01_pca, pca12 = timings_alpha12_pca)
    res_df <- data.frame(method = "rational", alpha = c(rep("01", length(m)*length(N)),rep("12", length(m)*length(N))), N = rep(N, each = length(m)), Time = c(unlist(res[["rational01"]]), unlist(res[["rational12"]])), m = rep(m, length(N)))
    res_df_tmp <- data.frame(method = "PCA", alpha = c(rep("01", length(m)*length(N)),rep("12", length(m)*length(N))), N = rep(N, each = length(m)), Time = c(unlist(res[["pca01"]]), unlist(res[["pca12"]])), m = rep(m, length(N)))
    res_df <- rbind(res_df, res_df_tmp)
    return(res_df)
}

