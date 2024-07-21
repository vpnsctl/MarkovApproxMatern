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
# Predict on loc_full based on observations on obs_ind

predict_rat_markov <- function(true_pred, y, loc_full, obs_ind, m, nu, kappa, sigma, sigma_e,  print=TRUE, equally_spaced = FALSE){
time_run2_approx <- list()
error <- list()

for(i_m in m){
        time_run2_approx[[as.character(i_m)]] <- list()
        error[[as.character(i_m)]] <- list()
        start = Sys.time()
        if(nu < 0.5) {
            Qrat <- matern.rational.ldl(loc = loc_full, order = i_m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp =  "spline", equally_spaced = equally_spaced)    
        } else {
            Qrat <- matern.rational.ldl(loc = loc_full, order = i_m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp =  "spline", equally_spaced = equally_spaced)    
        }
        
        Q <- t(Qrat$L)%*%Qrat$D%*%Qrat$L
        end1 <- Sys.time()
        A_obs <- Qrat$A[obs_ind,]
        A_mat <- t(A_obs)
        Q_xgiveny <- Q + (A_mat%*% A_obs)/sigma_e^2
        post_y <- (A_mat %*% y)/sigma_e^2
        approx_mean <-  Qrat$A%*%solve(Q_xgiveny, A_mat%*%y/sigma_e^2)
        end2 = Sys.time()
        time_run2_approx[[as.character(i_m)]][["Build_Q"]] <- as.numeric(end1-start, units = "secs")
        time_run2_approx[[as.character(i_m)]][["Get_pred"]] <- as.numeric(end2-end1, units = "secs")
        d_loc <- diff(loc_full)
        d_loc <- c(d_loc[1], d_loc)
        error[[as.character(i_m)]][["L2"]] <- sqrt(sum(d_loc * (true_pred-approx_mean)^2))
        error[[as.character(i_m)]][["Sup"]] <- max(abs(true_pred-approx_mean))
        if(print){
            print(time_run2_approx[[as.character(i_m)]])
        }
    }
    result <- list(error = error, timing = time_run2_approx)
    return(result)
}

## Predict PCA

predict_PCA <- function(true_pred, y, loc_full, obs_ind, m, nu, kappa, sigma, sigma_e,  print=TRUE){
time_run2_approx <- list()
error <- list()

D_loc <- as.matrix(dist(loc_full))
cov_mat <- rSPDE::matern.covariance(h=D_loc,kappa=kappa,nu=nu,sigma=sigma)
eigen_cov <- eigen(cov_mat)
for(i_m in m){
    time_run2_approx[[as.character(i_m)]] <- list()
    error[[as.character(i_m)]] <- list()
    start = Sys.time()
    K <- eigen_cov$vec[,1:i_m]    
    D <- diag(eigen_cov$val[1:i_m])    
    cov_KL <- K%*%D%*%t(K)
    end1 <- Sys.time()
    cov_KL <- cov_KL[, obs_ind]
    K <- K[obs_ind, ]
    svd_K <- svd(K%*%sqrt(D))
    cov_KL_svd_U <- cov_KL %*% svd_K$u
    y_new <- t(svd_K$u) %*% y
    prec_nugget <- cov_KL_svd_U %*% Matrix::Diagonal(x = 1/(svd_K$d^2 + sigma_e^2)) 
    post_mean = prec_nugget%*%y_new
    end2 =Sys.time()
    time_run2_approx[[as.character(i_m)]][["Build_Q"]] <- as.numeric(end1-start, units = "secs")
    time_run2_approx[[as.character(i_m)]][["Get_pred"]] <- as.numeric(end2-end1, units = "secs")
    d_loc <- diff(loc_full)
    d_loc <- c(d_loc[1], d_loc)
    error[[as.character(i_m)]][["L2"]] <- sqrt(sum(d_loc * (true_pred-post_mean)^2))
    error[[as.character(i_m)]][["Sup"]] <- max(abs(true_pred-post_mean))
    if(print){
            print(time_run2_approx[[as.character(i_m)]])
    }
}
result <- list(error = error, timing = time_run2_approx)
return(result)
}



## Predict KL

predict_KL <- function(true_pred, y, loc_full, obs_ind, m, N_KL, nu, kappa, sigma, sigma_e,  print=TRUE){
time_run2_approx <- list()
error <- list()

tmp_loc <- rep(FALSE, length(loc_full))
tmp_loc[obs_ind] <- TRUE

large_KL <- seq(min(loc_full), max(loc_full), length.out = N_KL)
kl_loc <- c(loc_full, large_KL)
kl_loc <- unique(kl_loc)

D_loc <- as.matrix(dist(loc_full))
cov_mat <- rSPDE::matern.covariance(h=D_loc,kappa=kappa,nu=nu,sigma=sigma)
eigen_cov <- eigen(cov_mat)
for(i_m in m){
    time_run2_approx[[as.character(i_m)]] <- list()
    error[[as.character(i_m)]] <- list()
    start = Sys.time()    
    K <- eigen_cov$vec[1:N,1:i_m]    
    D <- diag(eigen_cov$val[1:i_m])    
    cov_KL <- K%*%D%*%t(K)
    end1 <- Sys.time()
    cov_KL <- cov_KL[, obs_ind]
    K <- K[obs_ind, ]
    svd_K <- svd(K%*%sqrt(D))
    cov_KL_svd_U <- cov_KL %*% svd_K$u    
        y_new <- t(svd_K$u) %*% y
        prec_nugget <- cov_KL_svd_U %*% Matrix::Diagonal(x = 1/(svd_K$d^2 + sigma_e^2)) 
        post_mean = prec_nugget%*%y_new
    end2 =Sys.time()
    time_run2_approx[[as.character(i_m)]][["Build_Q"]] <- as.numeric(end1-start, units = "secs")
    time_run2_approx[[as.character(i_m)]][["Get_pred"]] <- as.numeric(end2-end1, units = "secs")
    d_loc <- diff(loc_full)
    d_loc <- c(d_loc[1], d_loc)
    error[[as.character(i_m)]][["L2"]] <- sqrt(sum(d_loc * (true_pred-post_mean)^2))
    error[[as.character(i_m)]][["Sup"]] <- max(abs(true_pred-post_mean))
    if(print){
            print(time_run2_approx[[as.character(i_m)]])
    }
}
result <- list(error = error, timing = time_run2_approx)
return(result)
}


# predict NN 

predict_rat_NN <- function(true_pred, y, loc_full, obs_ind, m, nu, kappa, sigma, sigma_e,  print=TRUE){
time_run2_approx <- list()
error <- list()

sigma.e <- sigma_e
obs.ind <- obs_ind
n.obs <- length(obs.ind)

D <- as.matrix(dist(loc_full))     
for(i_m in m){
        start = Sys.time()    
        Qnn <- get.nnQ(loc = loc_full[obs.ind],kappa = kappa,nu = nu,sigma = sigma, n.nbr = i_m)
        end1 <- Sys.time()
        Qhat <- Qnn + Diagonal(n.obs)/sigma.e^2        
        mu.nn <- solve(Qhat, y/sigma.e^2)
        Bp <- get.nn.pred(loc = loc_full, kappa = kappa, nu = nu, sigma = sigma, n.nbr = i_m, S = obs.ind)$B
        mu.nn <- Bp%*%mu.nn
        end2 <- Sys.time()
        time_run2_approx[[as.character(i_m)]][["Build_Q"]] <- as.numeric(end1-start, units = "secs")
        time_run2_approx[[as.character(i_m)]][["Get_pred"]] <- as.numeric(end2-end1, units = "secs")
        d_loc <- diff(loc_full)
        d_loc <- c(d_loc[1], d_loc)
        error[[as.character(i_m)]][["L2"]] <- sqrt(sum(d_loc * (true_pred-mu.nn)^2))
        error[[as.character(i_m)]][["Sup"]] <- max(abs(true_pred-mu.nn))   
        if(print){
                print(time_run2_approx[[as.character(i_m)]])
        }
}
result <- list(error = error, timing = time_run2_approx)
return(result)
}


# L is length of the interval

compare_times <- function(N, n_obs, L, range, sigma, sigma_e, m_rat, m_nngp_fun, m_pca_fun){
    sigma_e <- 0.1
    timings_alpha01_rat <- list()
    timings_alpha12_rat <- list()
    timings_alpha23_rat <- list()

    timings_alpha01_nngp <- list()
    timings_alpha12_nngp <- list()
    timings_alpha23_nngp <- list()

    timings_alpha01_pca <- list()
    timings_alpha12_pca <- list()
    timings_alpha23_pca <- list()    


    if(length(N) != length(n_obs)){
        stop("N and n_obs must have the same length!")
    }

    for(ii in 1:length(N)){
        n_loc <- N[ii]
        n_loc_obs <- n_obs[ii]
        loc <- seq(0, L, length.out = n_loc)
        nu <- 0.3
        kappa <- sqrt(8*nu)/range
        obs.ind <- sort(sample(1:n_loc)[1:n_loc_obs])

        Y <- y <- sample_y(loc[obs.ind], nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e)
        mu <- true_pred(y, loc=loc[obs.ind], loc_pred = loc, nu = nu, kappa=kappa, sigma=sigma, sigma_e = sigma_e)

        #

        #our method
        timings_alpha01_rat[[as.character(n_loc)]] <- predict_rat_markov(true_pred = mu, y = Y, loc_full = loc, obs_ind = obs.ind, m = m_rat, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e,   print = FALSE)

        # nngp
        m_nngp <- m_nngp_fun(m_rat, nu+0.5)

        timings_alpha01_nngp[[as.character(n_loc)]] <- predict_rat_NN(true_pred = mu, y = Y, loc_full = loc, obs_ind = obs.ind, m = m_nngp, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e,   print = FALSE)

        # PCA
        m_pca <- m_pca_fun(m_rat, nu+0.5)        
        timings_alpha01_pca[[as.character(n_loc)]] <- predict_PCA(true_pred = mu, y = Y, loc_full = loc, obs_ind = obs.ind, m = m_pca, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e,   print = FALSE)       

        nu <- 1.2

        Y <- y <- sample_y(loc[obs.ind], nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e)
        mu <- true_pred(y, loc=loc[obs.ind], loc_pred = loc, nu = nu, kappa=kappa, sigma=sigma, sigma_e = sigma_e)

        kappa <- sqrt(8*nu)/range            
        timings_alpha12_rat[[as.character(n_loc)]] <- predict_rat_markov(true_pred = mu, y = Y, loc_full = loc, obs_ind = obs.ind, m = m_rat, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e,   print = FALSE)

        m_nngp <- m_nngp_fun(m_rat, nu+0.5)
        timings_alpha12_nngp[[as.character(n_loc)]] <- predict_rat_NN(true_pred = mu, y = Y, loc_full = loc, obs_ind = obs.ind, m = m_nngp, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e,   print = FALSE)

        m_pca <- m_pca_fun(m_rat, nu+0.5)
        timings_alpha12_pca[[as.character(n_loc)]] <- predict_PCA(true_pred = mu, y = Y, loc_full = loc, obs_ind = obs.ind, m = m_pca,nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e,   print = FALSE)        

        nu <- 2.2

        Y <- y <- sample_y(loc[obs.ind], nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e)
        mu <- true_pred(y, loc=loc[obs.ind], loc_pred = loc, nu = nu, kappa=kappa, sigma=sigma, sigma_e = sigma_e)

        kappa <- sqrt(8*nu)/range            
        timings_alpha23_rat[[as.character(n_loc)]] <- predict_rat_markov(true_pred = mu, y = Y, loc_full = loc, obs_ind = obs.ind, m = m_rat, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e,   print = FALSE)

        m_nngp <- m_nngp_fun(m_rat, nu+0.5)
        timings_alpha23_nngp[[as.character(n_loc)]] <- predict_rat_NN(true_pred = mu, y = Y, loc_full = loc, obs_ind = obs.ind, m = m_nngp, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e,   print = FALSE)

        m_pca <- m_pca_fun(m_rat, nu+0.5)
        timings_alpha23_pca[[as.character(n_loc)]] <- predict_PCA(true_pred = mu, y = Y, loc_full = loc, obs_ind = obs.ind, m = m_pca,  nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e,   print = FALSE)  

    }

    res <- list(rational01 = timings_alpha01_rat, rational12 = timings_alpha12_rat, rational23 = timings_alpha23_rat, nngp01 = timings_alpha01_nngp, nngp12 = timings_alpha12_nngp, nngp23 = timings_alpha23_nngp, pca01 = timings_alpha01_pca, pca12 = timings_alpha12_pca, pca23 = timings_alpha23_pca)

    # res_df <- data.frame(method = "rational", alpha = c(rep("01", length(m)*length(N)),rep("12", length(m)*length(N))), N = rep(N, each = length(m)), Time = c(unlist(res[["rational01"]]), unlist(res[["rational12"]])), m = rep(m, length(N)))
    # res_df_tmp <- data.frame(method = "nnGP", alpha = c(rep("01", length(m)*length(N)),rep("12", length(m)*length(N))), N = rep(N, each = length(m)), Time = c(unlist(res[["nngp01"]]), unlist(res[["nngp12"]])), m = rep(m, length(N)))
    # res_df <- rbind(res_df, res_df_tmp)
    # res_df_tmp <- data.frame(method = "KL", alpha = c(rep("01", length(m)*length(N)),rep("12", length(m)*length(N))), N = rep(N, each = length(m)), Time = c(unlist(res[["kl01"]]), unlist(res[["kl12"]])), m = rep(m, length(N)))
    # res_df <- rbind(res_df, res_df_tmp)
    # return(res_df)
    return(res)
}

# only with nngp

compare_times_nngp <- function(N, n_obs, L, range, sigma, sigma_e, m_rat, m_nngp_fun){
    sigma_e <- 0.1
    timings_alpha01_rat <- list()
    timings_alpha12_rat <- list()
    timings_alpha23_rat <- list()

    timings_alpha01_nngp <- list()
    timings_alpha12_nngp <- list()
    timings_alpha23_nngp <- list()


    if(length(N) != length(n_obs)){
        stop("N and n_obs must have the same length!")
    }

    for(ii in 1:length(N)){
        n_loc <- N[ii]
        n_loc_obs <- n_obs[ii]
        loc <- seq(0, L, length.out = n_loc)
        nu <- 0.3
        kappa <- sqrt(8*nu)/range
        obs.ind <- sort(sample(1:n_loc)[1:n_loc_obs])

        Y <- y <- sample_y(loc[obs.ind], nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e)
        mu <- true_pred(y, loc=loc[obs.ind], loc_pred = loc, nu = nu, kappa=kappa, sigma=sigma, sigma_e = sigma_e)

        #

        #our method

        timings_alpha01_rat[[as.character(n_loc)]] <- predict_rat_markov(true_pred = mu, y = Y, loc_full = loc, obs_ind = obs.ind, m = m_rat, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e,   print = FALSE)

        # nngp
        m_nngp <- m_nngp_fun(m_rat, nu+0.5)
        timings_alpha01_nngp[[as.character(n_loc)]] <- predict_rat_NN(true_pred = mu, y = Y, loc_full = loc, obs_ind = obs.ind, m = m_nngp, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e,   print = FALSE)
           
        nu <- 1.2

        Y <- y <- sample_y(loc[obs.ind], nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e)
        mu <- true_pred(y, loc=loc[obs.ind], loc_pred = loc, nu = nu, kappa=kappa, sigma=sigma, sigma_e = sigma_e)


        kappa <- sqrt(8*nu)/range            
        timings_alpha12_rat[[as.character(n_loc)]] <- predict_rat_markov(true_pred = mu, y = Y, loc_full = loc, obs_ind = obs.ind, m = m_rat, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e,   print = FALSE)

        m_nngp <- m_nngp_fun(m_rat, nu+0.5)
        timings_alpha12_nngp[[as.character(n_loc)]] <- predict_rat_NN(true_pred = mu, y = Y, loc_full = loc, obs_ind = obs.ind, m = m_nngp, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e,   print = FALSE)      

        nu <- 2.2

        Y <- y <- sample_y(loc[obs.ind], nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e)
        mu <- true_pred(y, loc=loc[obs.ind], loc_pred = loc, nu = nu, kappa=kappa, sigma=sigma, sigma_e = sigma_e)


        kappa <- sqrt(8*nu)/range            
        timings_alpha23_rat[[as.character(n_loc)]] <- predict_rat_markov(true_pred = mu, y = Y, loc_full = loc, obs_ind = obs.ind, m = m_rat, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e,   print = FALSE)

        m_nngp <- m_nngp_fun(m_rat, nu+0.5)
        timings_alpha23_nngp[[as.character(n_loc)]] <- predict_rat_NN(true_pred = mu, y = Y, loc_full = loc, obs_ind = obs.ind, m = m_nngp, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e,   print = FALSE)

    }

    res <- list(rational01 = timings_alpha01_rat, rational12 = timings_alpha12_rat, rational23 = timings_alpha23_rat, nngp01 = timings_alpha01_nngp, nngp12 = timings_alpha12_nngp, nngp23 = timings_alpha23_nngp)

    # res_df <- data.frame(method = "rational", alpha = c(rep("01", length(m)*length(N)),rep("12", length(m)*length(N))), N = rep(N, each = length(m)), Time = c(unlist(res[["rational01"]]), unlist(res[["rational12"]])), m = rep(m, length(N)))
    # res_df_tmp <- data.frame(method = "nnGP", alpha = c(rep("01", length(m)*length(N)),rep("12", length(m)*length(N))), N = rep(N, each = length(m)), Time = c(unlist(res[["nngp01"]]), unlist(res[["nngp12"]])), m = rep(m, length(N)))
    # res_df <- rbind(res_df, res_df_tmp)
    # res_df_tmp <- data.frame(method = "KL", alpha = c(rep("01", length(m)*length(N)),rep("12", length(m)*length(N))), N = rep(N, each = length(m)), Time = c(unlist(res[["kl01"]]), unlist(res[["kl12"]])), m = rep(m, length(N)))
    # res_df <- rbind(res_df, res_df_tmp)
    # return(res_df)
    return(res)
}



# only with KL

compare_times_pca <- function(N, n_obs, L, range, sigma, sigma_e, m_rat, m_nngp){
    sigma_e <- 0.1
    timings_alpha01_rat <- list()
    timings_alpha12_rat <- list()
    timings_alpha23_rat <- list()

    timings_alpha01_kl <- list()
    timings_alpha12_kl <- list()
    timings_alpha23_kl <- list()


    if(length(N) != length(n_obs)){
        stop("N and n_obs must have the same length!")
    }

    for(ii in 1:length(N)){
        n_loc <- N[ii]
        n_loc_obs <- n_obs[ii]
        loc <- seq(0, L, length.out = n_loc)
        nu <- 0.3
        kappa <- sqrt(8*nu)/range
        obs.ind <- sort(sample(1:n_loc)[1:n_loc_obs])

        Y <- y <- sample_y(loc[obs.ind], nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e)
        mu <- true_pred(y, loc=loc[obs.ind], loc_pred = loc, nu = nu, kappa=kappa, sigma=sigma, sigma_e = sigma_e)


        #

        #our method

        timings_alpha01_rat[[as.character(n_loc)]] <- predict_rat_markov(true_pred = mu, y = Y, loc_full = loc, obs_ind = obs.ind, m = m_rat, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e,   print = FALSE)

        # #KL
        m_pca <- m_pca_fun(m_rat, nu+0.5)
        timings_alpha01_pca[[as.character(n_loc)]] <- predict_PCA(true_pred = mu, y = Y, loc_full = loc, obs_ind = obs.ind, m = m_pca, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e,   print = FALSE)
           
        nu <- 1.2

        Y <- y <- sample_y(loc[obs.ind], nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e)
        mu <- true_pred(y, loc=loc[obs.ind], loc_pred = loc, nu = nu, kappa=kappa, sigma=sigma, sigma_e = sigma_e)


        kappa <- sqrt(8*nu)/range            
        timings_alpha12_rat[[as.character(n_loc)]] <- predict_rat_markov(true_pred = mu, y = Y, loc_full = loc, obs_ind = obs.ind, m = m_rat, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e,   print = FALSE)

        m_pca <- m_pca_fun(m_rat, nu+0.5)
        timings_alpha12_pca[[as.character(n_loc)]] <- predict_PCA(true_pred = mu, y = Y, loc_full = loc, obs_ind = obs.ind, m = m_pca, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e,   print = FALSE)        

        nu <- 2.2

        Y <- y <- sample_y(loc[obs.ind], nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e)
        mu <- true_pred(y, loc=loc[obs.ind], loc_pred = loc, nu = nu, kappa=kappa, sigma=sigma, sigma_e = sigma_e)


        kappa <- sqrt(8*nu)/range            
        timings_alpha23_rat[[as.character(n_loc)]] <- predict_rat_markov(true_pred = mu, y = Y, loc_full = loc, obs_ind = obs.ind, m = m_rat, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e,   print = FALSE)

        m_pca <- m_pca_fun(m_rat, nu+0.5)
        timings_alpha23_pca[[as.character(n_loc)]] <- predict_PCA(true_pred = mu, y = Y, loc_full = loc, obs_ind = obs.ind, m = m_pca, nu = nu, kappa = kappa, sigma = sigma, sigma_e = sigma_e,   print = FALSE)  

    }

    res <- list(rational01 = timings_alpha01_rat, rational12 = timings_alpha12_rat, rational23 = timings_alpha23_rat, pca01 = timings_alpha01_pca, pca12 = timings_alpha12_pca, pca23 = timings_alpha23_pca)

    # res_df <- data.frame(method = "rational", alpha = c(rep("01", length(m)*length(N)),rep("12", length(m)*length(N))), N = rep(N, each = length(m)), Time = c(unlist(res[["rational01"]]), unlist(res[["rational12"]])), m = rep(m, length(N)))
    # res_df_tmp <- data.frame(method = "nnGP", alpha = c(rep("01", length(m)*length(N)),rep("12", length(m)*length(N))), N = rep(N, each = length(m)), Time = c(unlist(res[["nngp01"]]), unlist(res[["nngp12"]])), m = rep(m, length(N)))
    # res_df <- rbind(res_df, res_df_tmp)
    # res_df_tmp <- data.frame(method = "KL", alpha = c(rep("01", length(m)*length(N)),rep("12", length(m)*length(N))), N = rep(N, each = length(m)), Time = c(unlist(res[["kl01"]]), unlist(res[["kl12"]])), m = rep(m, length(N)))
    # res_df <- rbind(res_df, res_df_tmp)
    # return(res_df)
    return(res)
}


