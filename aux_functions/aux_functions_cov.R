# rSPDE functions and coeffs

library(rSPDE)

# Sparse matrices

library(Matrix)

# Toeplitz

library(SuperGauss)

# get corresponding m to other method based on ours:

get_m <- function(nu, m, method = c("nngp", "statespace", "kl"), type = c("prediction", "simulation")){
  type <- type[[1]]
  alpha <- nu + 0.5
  if(type == "prediction"){
    if(method == "statespace"){
      return((m-1)*ceil(alpha) + 1)
    } else if(method == "nngp"){
      return(2*round((m+3)*sqrt(m) * ceil(alpha)^(3/2)))
    } else if(method == "kl"){
      return(9*round((m+1)*sqrt(m) * ceil(alpha)^(3/2)))
    } else{
      stop("method not implemented.")
    }
  } else{
    if(method == "statespace"){
      return(round(ceil(alpha) * (m^(1/3)-1)) + 1)
    } else if(method == "nngp"){
      m_nngp <- round(sqrt(9*m) * ceil(alpha+2)^(3/2))
      d_m <- diff(m_nngp)
      d_m <- ifelse(d_m == 0, 1, d_m)
      m_nngp <- c(m_nngp[1], m_nngp[1]+cumsum(d_m))
      return(m_nngp)
    } else if(method == "kl"){
      return(9*round((m+1)*sqrt(m) * ceil(alpha)^(3/2)))
    } else{
      stop("method not implemented.")
    }
  } 
}



## NN GP - nearest neighbor

get.nn <- function(loc,kappa,nu,sigma, n.nbr) {
    k <- length(loc)
    Bs <- Fs <- Fsi <- Diagonal(k)
    Fs[1,1] <- sigma^2
    Fsi[1,1] <- 1/sigma^2
    for(i in 2:k) {
        nbrs <- get.neighbor(i,n.nbr)
        Sigma.in <- matern.covariance(h = abs(loc[nbrs]-loc[i]),
                                      kappa = kappa, nu = nu, sigma = sigma)
        Sigma.nn <- matern.covariance(h = as.matrix(dist(loc[nbrs])),
                                      kappa = kappa, nu = nu, sigma = sigma)
        tmp <- solve(Sigma.nn, Sigma.in)
        Bs[i,nbrs] <- -t(tmp)
        Fs[i,i] <- sigma^2 - t(Sigma.in)%*%tmp
        Fsi[i,i] <- 1/Fs[i,i]
    }
    return(list(Bs = Bs, Fs = Fs, Fi = Fsi))
}

get.nnQ3 <- function(loc,kappa,nu,sigma, n.nbr) {
    tmp <- get.nn(loc,kappa,nu,sigma,n.nbr)
    
    return(t(tmp$Bs) %*% tmp$Fi %*%tmp$Bs)
}


get.neighbor <- function(i, n.nbr) {
    if(i-n.nbr <= 0) {
        if(i==1) {
            return(numeric())
        } else {
            return(1:(i-1))    
        }
        
    } else {
        return((i-n.nbr):(i-1))
    }
}

get.Fi <- function(i, Sigma, n.nbr) {
    nbrs <- get.neighbor(i,n.nbr)
    if(length(nbrs) == 0) {
        return(Sigma[i,i])
    } else  {
        return(Sigma[i,i] - Sigma[i,nbrs]%*%solve(Sigma[nbrs,nbrs], Sigma[i,nbrs]))
    }
    
}

get.Fs <- function(Sigma, n.nbr) {
    k <- dim(Sigma)[1]
    Fs.d <- rep(0,k)
    for(i in 1:k) {
        Fs.d[i] <- get.Fi(i, Sigma, n.nbr)
    }
    return(Diagonal(k, Fs.d))
}


get.Bi <- function(i, Sigma, n.nbr) {
    nbrs <- get.neighbor(i,n.nbr)
    return(t(solve(Sigma[nbrs,nbrs], Sigma[nbrs,i])))
}


get.Bij <- function(i,j, n.nbr, Sigma) {
    if(i==j) {
        return(1)
    } else if (j %in% get.neighbor(i,n.nbr)) {
        l <- which(j == get.neighbor(i,n.nbr))
        Bi <- get.Bi(i,Sigma, n.nbr)
        return(-Bi[l])
    } else {
        return(0)
    }
}


get.Bs <- function(Sigma, n.nbr) {
    k <- dim(Sigma)[1]
    Bs <- diag(k)
    for(i in 2:k) {
        for(j in max(i-1-n.nbr,1):(i-1)) {
            Bs[i,j] <- get.Bij(i,j,n.nbr, Sigma)    
        }
    }
    return(as(Bs,"TsparseMatrix"))
}

get.nnCov <- function(Sigma, n.nbr) {
    Bs <- solve(get.Bs(Sigma, n.nbr))
    Fs <- get.Fs(Sigma, n.nbr)
    return(Bs %*% Fs %*%t(Bs))
}

get.nnQ <- function(Sigma, n.nbr) {
    Bs <- get.Bs(Sigma, n.nbr)
    Fs <- solve(get.Fs(Sigma, n.nbr))
    return(t(Bs) %*% Fs %*%Bs)
}

get.Bi2 <- function(loc,kappa,nu,sigma,i, nbrs) {
    Sigma.in <- matern.covariance(h = abs(loc[nbrs]-loc[i]),
                                  kappa = kappa, nu = nu, sigma = sigma)
    Sigma.nn <- matern.covariance(h = as.matrix(dist(loc[nbrs])),
                                  kappa = kappa, nu = nu, sigma = sigma)
    
    return(t(solve(Sigma.nn, Sigma.in)))
}


get.Bij2 <- function(loc,kappa,nu,sigma,i,j, n.nbr) {
    nbrs <- get.neighbor(i,n.nbr)
    if(i==j) {
        return(1)
    } else if (j %in% nbrs) {
        l <- which(j == nbrs)
        Bi <- get.Bi2(loc,kappa,nu,sigma,i, nbrs)
        return(-Bi[l])
    } else {
        return(0)
    }
}


get.Bs2 <- function(loc,kappa,nu,sigma, n.nbr) {
    k <- length(loc)
    Bs <- Diagonal(k)
    for(i in 2:k) {
        for(j in max(i-1-n.nbr,1):(i-1)) {
            Bs[i,j] <- get.Bij2(loc,kappa,nu,sigma,i,j,n.nbr)    
        }
    }
    return(Bs)
}

get.Fi2 <- function(loc,kappa,nu,sigma,i, n.nbr) {
    nbrs <- get.neighbor(i,n.nbr)
    
    if(length(nbrs) == 0) {
        return(sigma^2)
    } else  {
        Sigma.in <- matern.covariance(h = abs(loc[nbrs]-loc[i]),
                                      kappa = kappa, nu = nu, sigma = sigma)
        Sigma.nn <- matern.covariance(h = as.matrix(dist(loc[nbrs])),
                                      kappa = kappa, nu = nu, sigma = sigma)
        return(sigma^2 - t(Sigma.in)%*%solve(Sigma.nn, Sigma.in))
    }
    
}

get.Fs2 <- function(loc,kappa,nu,sigma, n.nbr) {
    k <- length(loc)
    Fs.d <- rep(0,k)
    
    for(i in 1:k) {
        Fs.d[i] <- get.Fi2(loc,kappa,nu,sigma, i, n.nbr)
    }
    return(Diagonal(k, Fs.d))
}

get.nnCov2 <- function(loc,kappa,nu,sigma, n.nbr) {
    Bs <- solve(get.Bs2(loc,kappa,nu,sigma, n.nbr))
    Fs <- get.Fs2(loc,kappa,nu,sigma, n.nbr)
    return(Bs %*% Fs %*%t(Bs))
}

get.nnQ2 <- function(loc,kappa,nu,sigma, n.nbr) {
    Bs <- get.Bs2(loc,kappa,nu,sigma, n.nbr)
    Fs <- solve(get.Fs2(loc,kappa,nu,sigma, n.nbr))
    return(t(Bs) %*% Fs %*%Bs)
}

## Fourier 


library(pracma)
library(gsignal)
library(evmix)
mat.spec <- function(x,kappa,alpha) {
    A <- gamma(alpha)*sqrt(4*pi)*kappa^(2*(alpha-0.5))/(2*pi*gamma(alpha-0.5))
    return(A/((kappa^2 + x^2)^alpha))
}


sample.mat <- function(n,kappa,alpha) {
    c = mat.spec(0,kappa = kappa, alpha = alpha)/dcauchy(0, location = 0, scale = kappa)
    k = 0
    out <- rep(0,n)
    while(k<n) {
        X <- rcauchy(1,location = 0, scale = kappa)
        U1 <- runif(1)
        fx <- dcauchy(X, location = 0, scale = kappa)
        fy <- mat.spec(X, kappa = kappa, alpha = alpha)
        Y <- c*fx*U1     
        if(Y< fy) {
            k = k + 1
            out[k] <- X
        }
    }
    return(out)
}

ff.approx <- function(m, kappa, alpha, loc) {
    w <- sample.mat(m,kappa,alpha)
    b <- runif(m,0,2*pi)
    ZX <- matrix(0, nrow = m, ncol = length(loc))
    for(i in 1:m){
        ZX[i,] <- sqrt(2)*cos(w[i]*loc + b[i])/sqrt(m)
    }
    return(t(ZX)%*%ZX)
}



## state space 


spec.coeff <- function(kappa,alpha, n, nm = 0) {
    A <- gamma(alpha)*sqrt(4*pi)*kappa^(2*(alpha-0.5))/(2*pi*gamma(alpha-0.5))
    ca <- ceiling(alpha)
    a <- rep(0, ca + n + 1)
    for(k in 0:(ca + n)) {
        if(k==0) {
            a[k+1] <- (kappa^(2*(alpha+n)-2*k)/factorial(k))
        } else {
            a[k+1] <- (kappa^(2*(alpha+n)-2*k)/factorial(k))*prod(alpha + n + c(0:(-(k-1))))
        }
        
    }
    b <- rep(0,n+1-nm)
    for(i in 0:(n-nm)) {
        b[i+1] <- nchoosek(n-nm,i)*kappa^(2*(n-nm-i))
    }
    return(list(a=a, b=b*A))
}


ab2spec <- function(a,b,x, flim = 2) {
    nx = length(x)
    x_max =x[nx]
    
    x = seq(from = 0,to = x_max*flim,length.out = flim*nx-1)
    n = flim*nx-1
    ind = 0:(n-1)
    dx = x[2]-x[1]
    
    ds = 2*pi/(n*dx)
    c = -(n/2)*ds
    s = c + ind*ds
    Snum  <- Sden <- 0
    for (i in 1:length(b)) {
        Snum <- Snum + b[i]*s^(2*(i-1))    
    }
    for (i in 1:length(a)) {
        Sden <- Sden + a[i]*s^(2*(i-1))    
    }
    
    return(Snum/Sden)
}

mat.spec_ss <- function(kappa,alpha,x, flim = 2) {
    A <- gamma(alpha)*sqrt(4*pi)*kappa^(2*(alpha-0.5))/(2*pi*gamma(alpha-0.5))
    nx = length(x)
    x_max =x[nx]
    
    x = seq(from = 0,to = x_max*flim,length.out = flim*nx-1)
    n = flim*nx-1
    ind = 0:(n-1)
    dx = x[2]-x[1]
    
    ds = 2*pi/(n*dx)
    c = -(n/2)*ds
    s = c + ind*ds
    S <- A/((kappa^2 + s^2)^alpha)
    return(S)
}

S2cov <- function(S,x, flim = 2) {
    nx = length(x)
    x_max =x[nx]
    
    x = seq(from = 0,to = x_max*flim,length.out = flim*nx-1)
    n = flim*nx-1
    ind = 0:(n-1)
    dx = x[2]-x[1]
    
    ds = 2*pi/(n*dx)
    c = -(n/2)*ds
    s = c + ind*ds
    
    fact_s = exp(-1i*(s-min(s))*x[1])
    
    phi_fft = fft(fact_s*S)
    C = Real(ds*exp(-1i*c*x)*phi_fft)
    #C = Real(exp(-1i*c*x)*phi_fft)/ds
    return(C[1:nx])
}


## KL

library(Rcpp)
sourceCpp("aux_functions/dist2mat.cpp")

dist2matR <- function(dist, buffer = 128){
    dist2mat(x = dist, bf = buffer)
}

build_KL <- function(covmat = NULL, order, eigen_cov = NULL){
    if(is.null(eigen_cov)){
        if(is.null(covmat)){
            stop("if eigen_cov is null, covmat must be nonnull!")
        }
        eigen_cov <- eigen(covmat)
        if(order>nrow(covmat)){
          order <- nrow(covmat)
        }
        K <- eigen_cov$vec[,1:order]
        D <- diag(eigen_cov$val[1:order])
        return(K%*%D%*%t(K))
    } else {
        if(order>nrow(eigen_cov$vec)){
          order <- nrow(eigen_cov$vec)
        }      
        K <- eigen_cov$vec[,1:order]
        D <- diag(eigen_cov$val[1:order])
        return(K%*%D%*%t(K))
    }
}

KL_matern <- function(m, loc = NULL, nu = NULL, kappa = NULL, sigma = NULL, eigen_cov = NULL){
    if(is.null(eigen_cov)){
        if(any(c(is.null(loc), is.null(nu),is.null(kappa),is.null(sigma)))){
            stop("if eigen_cov is null, nu, kappa, sigma and loc must be nonnull!")
        }
        D <- dist2matR(dist(loc))
        return(build_KL(rSPDE::matern.covariance(h=D,kappa=kappa,nu=nu,sigma=sigma), eigen_cov = eigen_cov))
    } else{
        return(build_KL(order = m, eigen_cov = eigen_cov))
    }
}

## Taper


wendland_K1 <- function(h, gamma){
    tmp <- (1-h/gamma)
    return(tmp^4 * (tmp>0) * (1+4*h/gamma))
}

wendland_K2 <- function(h, gamma){
    tmp <- (1-h/gamma)
    return(tmp^6 * (tmp>0) * (1+6*h/gamma + 35*h^2 /(3*gamma^2)))
}


## Assuming equally spaced mesh
taper_matern <- function(m, loc = NULL, delta_loc = NULL, nu, kappa = NULL, sigma = NULL, Sigma = NULL){
    if(is.null(loc)){
        if(is.null(delta_loc)){
            stop("delta_loc must be provided if loc is null.")
        }
    } else{
        sorted_loc <- sort(loc)
        delta_loc <- sorted_loc[2] - sorted_loc[1]
        rm(sorted_loc)
    }
    gamma <- (m+1) * delta_loc
    if(nu >= 2){
        stop("nu must be less than 2.")
    }
    if(is.null(Sigma)){
        if(any(c(is.null(nu), is.null(kappa),is.null(sigma)))){
            stop("if Sigma is null, nu, kappa and sigma must be nonnull!")
        }
        D <- dist2matR(dist(loc))
        Sigma <- rSPDE::matern.covariance(h=D,kappa=kappa,nu=nu,sigma=sigma)
    } 
    if(nu < 1){
        Sigma_tap <- wendland_K1(D, gamma = gamma)
        return(Sigma*Sigma_tap)
    } else{
        Sigma_tap <- wendland_K2(D, gamma = gamma)
        return(Sigma*Sigma_tap)
    }
}
