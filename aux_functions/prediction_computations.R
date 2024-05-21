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

# Example:
s <- seq(0,1,by=0.001)
nu <- 0.6
kappa <- 10
sigma <- 2
sim <- sample_matern(loc = s, nu = nu, kappa = kappa, sigma = sigma, nsim = 10000)
library(rSPDE)
c.true <- matern.covariance(0.5-s, kappa=kappa, nu=nu, sigma=sigma)
plot(s, c.true,
     type = "l", ylab = "C(|s-0.5|)", xlab = "s", ylim = c(0, 5),
     cex.main = 0.8, cex.axis = 0.8, cex.lab = 0.8
   )
lines(s, cov(t(sim))[(length(s)-1)/2+1,], col = 2)


sample_y <- function(loc, nu, kappa, sigma, sigma_e, seed=123){
    set.seed(seed)
    z <- sample_matern(loc = loc, nu = nu, kappa = kappa, sigma = sigma, nsim = 1)
    return(z + sigma_e * rnorm(length(z)))
}

# Example:
y <- sample_y(s,nu,kappa,sigma,0.1, 1)


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

# Example:
post_mean_true <- true_pred(y, loc=s, nu = nu, kappa=kappa, sigma=sigma, sigma_e = 0.1)

#### Predict error computation


## Predict rational markov (our method)

pred_rat_markov <- function(y, loc, m, nu, kappa, sigma, sigma_e, equally_spaced = FALSE){
N <- length(loc)
pred <- list()
    for(i_m in m){
        if(nu < 0.5){
            r <- rSPDE::matern.rational.precision(loc, order = i_m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp = "spline", equally_spaced = equally_spaced, ordering = "location")
        } else {
            r <- rSPDE::matern.rational.precision(loc, order = i_m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "chebfun", type_interp = "spline", equally_spaced = equally_spaced, ordering = "location")
        } 
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

# Example:
start <- Sys.time()
post_mean_rat <- pred_rat_markov(y, loc=s, m = 1:6, nu = nu, kappa=kappa, sigma=sigma, sigma_e = 0.1, equally_spaced = TRUE)
end <- Sys.time()
end - start

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
    print(i_m)
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

# Example:
post_mean_pca <- pred_PCA(y, loc=s, m = 1:6, nu = nu, kappa=kappa, sigma=sigma, sigma_e = 0.1)


# predict NN 

pred_rat_NN <- function(loc, m, nu, kappa, sigma, sigma_e, samples, print=TRUE){
N <- length(loc)

D <- dist2matR(dist(loc))
Sigma <- rSPDE::matern.covariance(h=D,kappa=kappa,nu=nu,sigma=sigma)   
pred <- list()     
rat_m <- m
m <- get_m(nu = nu, m = m, method = "nngp", type = "prediction")
count <- 1
for(i_m in m){
        # prec_mat <- get.nnQ(Sigma, i_m)
        prec_mat <- get.nnQ2(loc=loc,kappa=kappa,nu=nu,sigma=sigma, n.nbr = i_m)
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

# Example:
start <- Sys.time()
post_mean_nn <- pred_rat_NN(y, loc=s, m = 1:6, nu = nu, kappa=kappa, sigma=sigma, sigma_e = 0.1)
end <- Sys.time()
end-start
