source("aux_functions/aux_functions_cov.R")

sample_matern <- function(loc, nu, kappa, sigma, nsim = 1){
    loc <- loc - loc[1]
    N <- length(loc)
    acf <- rSPDE::matern.covariance(h=loc,kappa=kappa,nu=nu,sigma=sigma)
    z <- matrix(rnorm(N * nsim), ncol = N, nrow = nsim)
    cov_mat <- toeplitz(as.vector(acf), symmetric=TRUE)
    L <- chol(cov_mat)
    return(t(L)%*%t(z))
}

# # Example:
# s <- seq(0,1,by = 0.001)
# nu <- 0.6
# kappa <- 10
# sigma <- 1
# sim <- sample_matern(loc = s, nu = nu, kappa = kappa, sigma = sigma, nsim = 10000)
# library(rSPDE)
# c.true <- matern.covariance(0.5-s, kappa=kappa, nu=nu, sigma=sigma)
# plot(s, c.true,
#      type = "l", ylab = "C(|s-0.5|)", xlab = "s", ylim = c(0, 5),
#      cex.main = 0.8, cex.axis = 0.8, cex.lab = 0.8
#    )
# lines(s, cov(t(sim))[(length(s)-1)/2+1,], col = 2)


sample_y <- function(loc, nu, kappa, sigma, sigma_e, seed=123){
    set.seed(seed)
    z <- sample_matern(loc = loc, nu = nu, kappa = kappa, sigma = sigma, nsim = 1)
    return(z + sigma_e * rnorm(length(z)))
}

# # Example:
# y <- sample_y(s,nu,kappa,sigma,0.1, 1)


true_pred <- function(y, loc, nu, kappa, sigma, sigma_e){
    loc <- loc - loc[1]
    acf = rSPDE::matern.covariance(h=loc,kappa=kappa,nu=nu,sigma=sigma)
    acf <- as.vector(acf)
    acf[1] = acf[1]+sigma_e^2
    acf2 =rSPDE::matern.covariance(h=loc,kappa=kappa,nu=nu,sigma=sigma)
    acf2 =as.vector(acf2)
    cov_mat_nugget <-  toeplitz(as.vector(acf), symmetric=TRUE)
    cov_mat <-  toeplitz(as.vector(acf2), symmetric=TRUE)
    d <- solve(cov_mat_nugget, y)
    return(cov_mat%*%d)
}

# # Example:
# post_mean_true <- true_pred(y, loc=s, nu = nu, kappa=kappa, sigma=sigma, sigma_e = 0.1)

#### Predict error computation


## Predict rational markov (our method)

pred_rat_markov <- function(y, loc, m, nu, kappa, sigma, sigma_e, equally_spaced = FALSE){
N <- length(loc)
pred <- list()
    for(i_m in m){
        if(nu < 0.5){
            # r <- rSPDE::matern.rational.precision(loc, order = i_m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp = "spline", equally_spaced = equally_spaced, ordering = "location")
            r <- rSPDE::matern.rational.precision(loc, order = i_m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp = "spline", equally_spaced = equally_spaced, ordering = "field")            
        } else {
            # r <- rSPDE::matern.rational.precision(loc, order = i_m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "chebfun", type_interp = "spline", equally_spaced = equally_spaced, ordering = "location")
            r <- rSPDE::matern.rational.precision(loc, order = i_m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp = "spline", equally_spaced = equally_spaced, ordering = "field")            
        } 
            r$Q <- (r$Q + t(r$Q))/2
            A_mat = t(r$A)
            Q_xgiveny <-(A_mat%*% (r$A))/sigma_e^2 + r$Q
            post_y <- (A_mat%*% y)/sigma_e^2
            # R <- Matrix::Cholesky(Q_xgiveny, perm=FALSE)
            R <- Matrix::Cholesky(Q_xgiveny)            
            mu_xgiveny <- solve(R, post_y, system = "A")
            approx_mean1 <-  r$A %*% mu_xgiveny
            pred[[as.character(i_m)]] <- approx_mean1
        }
    return(pred)
}

# # Example:
# start <- Sys.time()
# post_mean_rat <- pred_rat_markov(y, loc=s, m = 1:6, nu = nu, kappa=kappa, sigma=sigma, sigma_e = 0.1, equally_spaced = TRUE)
# end <- Sys.time()
# end - start



## Predict rational markov (our method) - LDL version

pred_rat_markov_ldl <- function(y, loc, m, nu, kappa, sigma, sigma_e, equally_spaced = FALSE){
N <- length(loc)
pred <- list()
    for(i_m in m){
        if(nu < 0.5){
            r <- rSPDE::matern.rational.ldl(loc, order = i_m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp = "spline", equally_spaced = equally_spaced, ordering = "location")
        } else {
            r <- rSPDE::matern.rational.ldl(loc, order = i_m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "chebfun", type_interp = "spline", equally_spaced = equally_spaced, ordering = "location")
        } 
            r$Q <- t(r$L) %*% t(r$D) %*% r$L
            r$Q <- (r$Q + t(r$Q))/2
            A_mat = t(r$A)
            Q_xgiveny <-(A_mat%*% (r$A))/sigma_e^2 + r$Q
            post_y <- (A_mat%*% y)/sigma_e^2
            R <- Matrix::Cholesky(Q_xgiveny, perm=FALSE)
            mu_xgiveny <- solve(R, post_y, system = "A")
            approx_mean1 <-  r$A %*% mu_xgiveny
            pred[[as.character(i_m)]] <- approx_mean1
        }
    return(pred)
}

# # Example:
# start <- Sys.time()
# post_mean_rat <- pred_rat_markov(y, loc=s, m = 1:6, nu = nu, kappa=kappa, sigma=sigma, sigma_e = 0.1, equally_spaced = TRUE)
# end <- Sys.time()
# end - start

## Predict PCA

pred_PCA <- function(y, loc, m, nu, kappa, sigma, sigma_e, method = c("woodbury", "approx")){
N <- length(loc)
method <- method[[1]]
pred <- list()
D_loc <- dist2matR(dist(loc))
cov_mat <- rSPDE::matern.covariance(h=D_loc,kappa=kappa,nu=nu,sigma=sigma)
eigen_cov <- eigen(cov_mat)
rat_m <- m
m <- get_m(nu = nu, m = m, method = "kl", type = "prediction")
count <- 1
for(i_m in m){
    i_m <- min(i_m, N/2)
    K <- eigen_cov$vec[,1:i_m]    
    D <- diag(eigen_cov$val[1:i_m])    
    if(method == "woodbury"){
        K <- K%*%sqrt(D)
        tK <- t(K)
        I_mat_low <- Matrix::Diagonal(x = 1, n = ncol(K))
        I_mat_high <- Matrix::Diagonal(x = 1, n = nrow(K))
        diag_eps_inv <- Matrix::Diagonal(x = sigma_e^(-2), n = nrow(K))
        tKSig <- tK %*% diag_eps_inv
        SigK <- t(tKSig)
        inv_Part <- solve(I_mat_low + tK%*%SigK)
        nugget_part <- inv_Part %*% tKSig
        cov_mat_nugget_inv <- SigK %*% nugget_part
        solve_nugget <- diag_eps_inv - cov_mat_nugget_inv
        post_mean <- tK %*% solve_nugget %*% y
        post_mean <- K%*%post_mean
    } else{
        cov_KL <- K%*%D%*%t(K)
        svd_K <- svd(K%*%sqrt(D))
        cov_KL_svd_U <- cov_KL %*% svd_K$u
        y_new <- t(svd_K$u) %*% y
        prec_nugget <- cov_KL_svd_U %*% Matrix::Diagonal(x = 1/(svd_K$d^2 + sigma_e^2)) 
        post_mean = prec_nugget%*%y_new
    }
    pred[[as.character(rat_m[count])]] <- post_mean
    count <- count + 1
}
return(pred)
}

# # Example:
# post_mean_pca <- pred_PCA(y, loc=s, m = 1:6, nu = nu, kappa=kappa, sigma=sigma, sigma_e = 0.1)


# predict NN 

pred_rat_NN <- function(y, loc, m, nu, kappa, sigma, sigma_e){
N <- length(loc)

pred <- list()     
rat_m <- m
m <- get_m(nu = nu, m = m, method = "nngp", type = "prediction")
count <- 1
for(i_m in m){
        # prec_mat <- get.nnQ(Sigma, i_m)
        prec_mat <- get.nnQ(loc=loc,kappa=kappa,nu=nu,sigma=sigma, n.nbr = i_m)
        I <- Matrix::Diagonal(n = ncol(prec_mat), x = 1)
        Q_xgiveny <- I * 1/sigma_e^2 + prec_mat
        post_y <- y/sigma_e^2
        R <- Matrix::Cholesky(Q_xgiveny, perm=FALSE)
        mu_xgiveny <- solve(R, post_y, system = "A")
        pred[[as.character(rat_m[count])]] <- mu_xgiveny
        count <- count + 1
    }
    return(pred)
}

# # Example:
# start <- Sys.time()
# post_mean_nn <- pred_rat_NN(y, loc=s, m = 1:6, nu = nu, kappa=kappa, sigma=sigma, sigma_e = 0.1)
# end <- Sys.time()
# end-start


# Predict Fourier
# Implemented for sigma = 1

pred_Fourier <- function(y, loc, m, nu, kappa, sigma_e,samples = 100){
N <- length(loc)
pred <- list()
D_loc <- dist2matR(dist(loc))
rat_m <- m
m <- get_m(nu = nu, m = m, method = "kl", type = "prediction")
count <- 1
for(i_m in m){
    post_mean <- rep(0, length(y))
    for(jj in 1:samples){
        K <-  ff.comp(m = i_m, kappa = kappa, alpha = nu + 0.5, loc = loc)
        tK <- t(K)
        I_mat_low <- Matrix::Diagonal(x = 1, n = ncol(K))
        I_mat_high <- Matrix::Diagonal(x = 1, n = nrow(K))
        diag_eps_inv <- Matrix::Diagonal(x = sigma_e^(-2), n = nrow(K))
        tKSig <- tK %*% diag_eps_inv
        SigK <- t(tKSig)
        inv_Part <- solve(I_mat_low + tK%*%SigK)
        nugget_part <- inv_Part %*% tKSig
        cov_mat_nugget_inv <- SigK %*% nugget_part
        solve_nugget <- diag_eps_inv - cov_mat_nugget_inv
        post_mean_tmp <- tK %*% solve_nugget %*% y
        post_mean <- post_mean + K%*%post_mean_tmp
    }
    pred[[as.character(rat_m[count])]] <- post_mean/samples
    count <- count + 1
}
return(pred)
}

# # Example:
# post_mean_fourier <- pred_Fourier(y, loc=s, m = 1:6, nu = nu, kappa=kappa, sigma_e = 0.1)


# Predict SS
# loc in [0,1]

pred_statespace <- function(y, loc, m, nu, kappa, sigma_e, L=1, flim = 2, fact = 100){
N <- length(loc)
ind = 1 + fact*(0:(N-1))
h2 = seq(from=0,to=L,length.out=fact*(N-1)+1)
pred <- list()     
rat_m <- m
m <- get_m(nu = nu, m = m, method = "statespace", type = "prediction")
count <- 1
for(i_m in m){
        coeff <- spec.coeff(kappa = kappa,alpha = nu + 0.5,i_m)
        S1 <- ab2spec(coeff$a,coeff$b,h2, flim = flim)
        r1 <- S2cov(S1,h2,flim = flim)
        acf <- r1[ind]
        cov_mat <- toeplitz(acf, symmetric=TRUE)
        acf2 <- acf
        acf2[1] <- acf2[1] + sigma_e^2
        cov_mat_nugget <-  toeplitz(as.vector(acf2), symmetric=TRUE)
        d <- solve(cov_mat_nugget, y)
        pred[[as.character(rat_m[count])]] <- cov_mat%*%d
        count <- count + 1
    }
    return(pred)
}

# # Example:
# post_mean_ss <- pred_statespace(y, loc=s, m = 1:6, nu = nu, kappa=kappa, sigma_e = 0.1)


# Generate complete table
# sigma = 1
# L is the length of the interval

compute_pred_errors <- function(N, range, nu.vec, m.vec, sigma_e, L = 1, seed = 123, print = FALSE){
    post_mean_true <- list()
    post_mean_rat <- list()
    post_mean_rat_ldl <- list()    
    post_mean_PCA <- list()
    post_mean_nnGP <- list()
    post_mean_fourier <- list()
    post_mean_statespace <- list()

    #Rational

    for(n_loc in N){
        loc <- seq(0,L,length.out = n_loc)

        post_mean_true[[as.character(n_loc)]] <- list()
        post_mean_rat[[as.character(n_loc)]] <- list()
        post_mean_rat_ldl[[as.character(n_loc)]] <- list()        
        post_mean_PCA[[as.character(n_loc)]] <- list()
        post_mean_nnGP[[as.character(n_loc)]] <- list()
        post_mean_fourier[[as.character(n_loc)]] <- list()
        post_mean_statespace[[as.character(n_loc)]] <- list()

        for(i in 1:length(nu.vec)) {
            cat(i/length(nu.vec)," ")
            nu <- nu.vec[i]
            alpha <- nu + 0.5  
            kappa <- sqrt(8*nu)/range
            y <- sample_y(loc = loc,nu = nu,kappa = kappa ,sigma = 1, sigma_e = sigma_e, seed = seed)
            
            if(print){
                message("Starting true posterior")
            }
            post_mean_true[[as.character(n_loc)]][[as.character(nu)]] <- true_pred(y=y, loc=loc, nu=nu, kappa=kappa, sigma=1, sigma_e=sigma_e)
            if(print){
                message("Starting rational posterior")
            }            
            post_mean_rat[[as.character(n_loc)]][[as.character(nu)]] <- pred_rat_markov(y=y, loc=loc, m=m.vec, nu=nu, kappa=kappa, sigma=1, sigma_e=sigma_e, equally_spaced = TRUE)
            if(print){
                message("Starting rational ldl posterior")
            }            
            post_mean_rat_ldl[[as.character(n_loc)]][[as.character(nu)]] <- pred_rat_markov(y=y, loc=loc, m=m.vec, nu=nu, kappa=kappa, sigma=1, sigma_e=sigma_e, equally_spaced = TRUE)            
            if(print){
                message("Starting PCA posterior")
            }
            post_mean_PCA[[as.character(n_loc)]][[as.character(nu)]] <- pred_PCA(y=y, loc=loc, m=m.vec, nu=nu, kappa=kappa, sigma=1, sigma_e = sigma_e)  
            if(print){
                message("Starting nnGP posterior")
            }
            post_mean_nnGP[[as.character(n_loc)]][[as.character(nu)]] <- pred_rat_NN(y = y, loc=loc, m=m.vec, nu=nu, kappa=kappa, sigma=1, sigma_e=sigma_e)
            if(print){
                message("Starting Fourier posterior")
            }
            post_mean_fourier[[as.character(n_loc)]][[as.character(nu)]] <- pred_Fourier(y=y, loc=loc, m=m.vec, nu=nu, kappa=kappa, sigma_e=sigma_e,samples = 100)
            if(print){
                message("Starting state-space posterior")
            }
            post_mean_statespace[[as.character(n_loc)]][[as.character(nu)]] <- pred_statespace(y=y, loc=loc, m=m.vec, nu=nu, kappa=kappa, sigma_e=sigma_e, flim = 2, fact = 100)
        }
    }

    df_pred <- data.frame(Method = "Rational", nu = nu.vec[1], Norm = "l2", N = N[1], m = m.vec)
    tmp_true <- post_mean_true[[as.character(N[[1]])]][[as.character(nu[[1]])]]
    error_l2_tmp <- unlist(lapply(post_mean_rat[[as.character(N[1])]][[as.character(nu.vec[1])]], function(y_hat){norm(as.vector(y_hat)-tmp_true)}))
    error_max_tmp <- unlist(lapply(post_mean_rat[[as.character(N[1])]][[as.character(nu.vec[1])]], function(y_hat){max(abs(as.vector(y_hat)-tmp_true))}))    
    df_pred[["Error"]] <- error_l2_tmp
    df_tmp <- data.frame(Method = "Rational", nu = nu.vec[1], Norm = "max", N = N[1], m = m.vec)
    df_tmp[["Error"]] <- error_max_tmp
    df_pred <- rbind(df_pred, df_tmp)
    if(length(nu.vec)>1){
        for(j in 2:length(nu.vec)){
                df_tmp <- data.frame(Method = "Rational", nu = nu.vec[j], Norm = "l2", N = N[1], m = m.vec)
                tmp_true <- post_mean_true[[as.character(N[1])]][[as.character(nu.vec[j])]]
                error_l2_tmp <- unlist(lapply(post_mean_rat[[as.character(N[1])]][[as.character(nu.vec[j])]], function(y_hat){norm(as.vector(y_hat)-tmp_true)}))
                error_max_tmp <- unlist(lapply(post_mean_rat[[as.character(N[1])]][[as.character(nu.vec[j])]], function(y_hat){max(abs(as.vector(y_hat)-tmp_true))}))    
                df_tmp[["Error"]] <- error_l2_tmp
                df_pred <- rbind(df_pred, df_tmp)
                df_tmp <- data.frame(Method = "Rational", nu = nu.vec[j], Norm = "max", N = N[1], m = m.vec)
                df_tmp[["Error"]] <- error_max_tmp
                df_pred <- rbind(df_pred, df_tmp)            
        }
    }
    if(length(N)>1){
        for(i in 2:length(N)){
            for(j in 1:length(nu.vec)){
                    df_tmp <- data.frame(Method = "Rational", nu = nu.vec[j], Norm = "l2", N = N[i], m = m.vec)
                    tmp_true <- post_mean_true[[as.character(N[i])]][[as.character(nu.vec[j])]]
                    error_l2_tmp <- unlist(lapply(post_mean_rat[[as.character(N[i])]][[as.character(nu.vec[j])]], function(y_hat){norm(as.vector(y_hat)-tmp_true)}))
                    error_max_tmp <- unlist(lapply(post_mean_rat[[as.character(N[i])]][[as.character(nu.vec[j])]], function(y_hat){max(abs(as.vector(y_hat)-tmp_true))}))    
                    df_tmp[["Error"]] <- error_l2_tmp
                    df_pred <- rbind(df_pred, df_tmp)
                    df_tmp <- data.frame(Method = "Rational", nu = nu.vec[j], Norm = "max", N = N[i], m = m.vec)
                    df_tmp[["Error"]] <- error_max_tmp
                    df_pred <- rbind(df_pred, df_tmp)            
            }
        }
    }

    ## LDL
    for(i in 1:length(N)){
        for(j in 1:length(nu.vec)){
                df_tmp <- data.frame(Method = "Rational_LDL", nu = nu.vec[j], Norm = "l2", N = N[i], m = m.vec)
                tmp_true <- post_mean_true[[as.character(N[i])]][[as.character(nu.vec[j])]]
                error_l2_tmp <- unlist(lapply(post_mean_rat_ldl[[as.character(N[i])]][[as.character(nu.vec[j])]], function(y_hat){norm(as.vector(y_hat)-tmp_true)}))
                error_max_tmp <- unlist(lapply(post_mean_rat_ldl[[as.character(N[i])]][[as.character(nu.vec[j])]], function(y_hat){max(abs(as.vector(y_hat)-tmp_true))}))    
                df_tmp[["Error"]] <- error_l2_tmp
                df_pred <- rbind(df_pred, df_tmp)
                df_tmp <- data.frame(Method = "Rational_LDL", nu = nu.vec[j], Norm = "max", N = N[i], m = m.vec)
                df_tmp[["Error"]] <- error_max_tmp
                df_pred <- rbind(df_pred, df_tmp)            
        }
    }    

    ## PCA
    for(i in 1:length(N)){
        for(j in 1:length(nu.vec)){
                df_tmp <- data.frame(Method = "PCA", nu = nu.vec[j], Norm = "l2", N = N[i], m = m.vec)
                tmp_true <- post_mean_true[[as.character(N[i])]][[as.character(nu.vec[j])]]
                error_l2_tmp <- unlist(lapply(post_mean_PCA[[as.character(N[i])]][[as.character(nu.vec[j])]], function(y_hat){norm(as.vector(y_hat)-tmp_true)}))
                error_max_tmp <- unlist(lapply(post_mean_PCA[[as.character(N[i])]][[as.character(nu.vec[j])]], function(y_hat){max(abs(as.vector(y_hat)-tmp_true))}))    
                df_tmp[["Error"]] <- error_l2_tmp
                df_pred <- rbind(df_pred, df_tmp)
                df_tmp <- data.frame(Method = "PCA", nu = nu.vec[j], Norm = "max", N = N[i], m = m.vec)
                df_tmp[["Error"]] <- error_max_tmp
                df_pred <- rbind(df_pred, df_tmp)            
        }
    }

    ## nnGP
    for(i in 1:length(N)){
        for(j in 1:length(nu.vec)){
                df_tmp <- data.frame(Method = "nnGP", nu = nu.vec[j], Norm = "l2", N = N[i], m = m.vec)
                tmp_true <- post_mean_true[[as.character(N[i])]][[as.character(nu.vec[j])]]
                error_l2_tmp <- unlist(lapply(post_mean_nnGP[[as.character(N[i])]][[as.character(nu.vec[j])]], function(y_hat){norm(as.vector(y_hat)-tmp_true)}))
                error_max_tmp <- unlist(lapply(post_mean_nnGP[[as.character(N[i])]][[as.character(nu.vec[j])]], function(y_hat){max(abs(as.vector(y_hat)-tmp_true))}))    
                df_tmp[["Error"]] <- error_l2_tmp
                df_pred <- rbind(df_pred, df_tmp)
                df_tmp <- data.frame(Method = "nnGP", nu = nu.vec[j], Norm = "max", N = N[i], m = m.vec)
                df_tmp[["Error"]] <- error_max_tmp
                df_pred <- rbind(df_pred, df_tmp)            
        }
    }

    ## State-space
    for(i in 1:length(N)){
        for(j in 1:length(nu.vec)){
                df_tmp <- data.frame(Method = "Fourier", nu = nu.vec[j], Norm = "l2", N = N[i], m = m.vec)
                tmp_true <- post_mean_true[[as.character(N[i])]][[as.character(nu.vec[j])]]
                error_l2_tmp <- unlist(lapply(post_mean_statespace[[as.character(N[i])]][[as.character(nu.vec[j])]], function(y_hat){norm(as.vector(y_hat)-tmp_true)}))
                error_max_tmp <- unlist(lapply(post_mean_statespace[[as.character(N[i])]][[as.character(nu.vec[j])]], function(y_hat){max(abs(as.vector(y_hat)-tmp_true))}))    
                df_tmp[["Error"]] <- error_l2_tmp
                df_pred <- rbind(df_pred, df_tmp)
                df_tmp <- data.frame(Method = "Fourier", nu = nu.vec[j], Norm = "max", N = N[i], m = m.vec)
                df_tmp[["Error"]] <- error_max_tmp
                df_pred <- rbind(df_pred, df_tmp)            
        }
    }

    ## Fourier
    for(i in 1:length(N)){
        for(j in 1:length(nu.vec)){
                df_tmp <- data.frame(Method = "State-Space", nu = nu.vec[j], Norm = "l2", N = N[i], m = m.vec)
                tmp_true <- post_mean_true[[as.character(N[i])]][[as.character(nu.vec[j])]]
                error_l2_tmp <- unlist(lapply(post_mean_fourier[[as.character(N[i])]][[as.character(nu.vec[j])]], function(y_hat){norm(as.vector(y_hat)-tmp_true)}))
                error_max_tmp <- unlist(lapply(post_mean_fourier[[as.character(N[i])]][[as.character(nu.vec[j])]], function(y_hat){max(abs(as.vector(y_hat)-tmp_true))}))    
                df_tmp[["Error"]] <- error_l2_tmp
                df_pred <- rbind(df_pred, df_tmp)
                df_tmp <- data.frame(Method = "State-Space", nu = nu.vec[j], Norm = "max", N = N[i], m = m.vec)
                df_tmp[["Error"]] <- error_max_tmp
                df_pred <- rbind(df_pred, df_tmp)            
        }
    }
    df_pred[["m"]] <- as.character(df_pred[["m"]])
    return(df_pred)
}




# Function to plot the L2 and Linfinity distances between the covariance functions
# If methods is null, all methods will be considered
# If n_loc is null, the first value of n_loc will be considered. n_loc must be a character
# n_loc is the number of locations

library(ggplot2)

plot_pred <- function(df_pred, n_loc = NULL, norm = c("l2", "max"), methods = NULL, logscale = TRUE){
    if(is.null(n_loc)){
        n_loc <- df_pred[["N"]][1]
    }
    n_loc <- as.character(n_loc)
    df_pred <- df_pred |> dplyr::filter(N == n_loc)
    if(is.null(methods)){
        methods <- unique(df_pred[["Method"]])
    }
    norm <- norm[[1]]
    if(logscale){
        return(df_pred |> dplyr::filter(Norm == norm, Method %in% methods) |> ggplot2::ggplot() + 
            geom_line(aes(x = nu, y = Error, col = m,linetype=Method)) + scale_y_log10())
    } else{
        return(df_pred |> dplyr::filter(Norm == norm, Method %in% methods) |> ggplot2::ggplot() + 
            geom_line(aes(x = nu, y = Error, col = m,linetype=Method)))
    }
}
