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
      return(max((m-1)*ceil(alpha) + 1,1))
    } else if(method == "nngp"){
      return(max(2*round((m+3)*sqrt(m) * ceil(alpha)^(3/2)),1))
    } else if(method == "kl"){
      return(max(9*round((m+1)*sqrt(m) * ceil(alpha)^(3/2)),2))
    } else{
      stop("method not implemented.")
    }
  } else{
    if(method == "statespace"){
      return(max(round(ceil(alpha) * (m^(1/3)-1)) + 1,1))
    } else if(method == "nngp"){
      m_nngp <- round(sqrt(9*m) * ceil(alpha+2)^(3/2))
      d_m <- diff(m_nngp)
      d_m <- ifelse(d_m == 0, 1, d_m)
      m_nngp <- c(m_nngp[1], m_nngp[1]+cumsum(d_m))
      return(max(m_nngp,1))
    } else if(method == "kl"){
      return(max(9*round((m+1)*sqrt(m) * ceil(alpha)^(3/2)),2))
    } else{
      stop("method not implemented.")
    }
  } 
}

## NN GP - nearest neighbor

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

get.nn <- function(loc,kappa,nu,sigma, n.nbr) {
    k <- length(loc)
    Bs <- Fs <- Fsi <- Matrix::Diagonal(k)
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

get.nnQ <- function(loc,kappa,nu,sigma, n.nbr) {
    tmp <- get.nn(loc,kappa,nu,sigma,n.nbr)
    
    return(t(tmp$Bs) %*% tmp$Fi %*%tmp$Bs)
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

ff.comp <- function(m, kappa, alpha, loc) {
    w <- sample.mat(m,kappa,alpha)
    b <- runif(m,0,2*pi)
    ZX <- matrix(0, nrow = m, ncol = length(loc))
    for(i in 1:m){
        ZX[i,] <- sqrt(2)*cos(w[i]*loc + b[i])/sqrt(m)
    }
    return(t(ZX))
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


## Parsimonious rational approximation

library(gsignal)

kappa_integral <- function(n, beta, kappa) {
    y <- 0
    for (k in 0:n) {
        y <- y + (-1)^(k + 1) * choose(n, k) / (n - k - beta + 1)
    }
    return(kappa^(2 * (n - beta + 1)) * y)
}

lapprox <- function(beta, kappa) {
    nu <- 2 * beta - 1 / 2
    alpha <- nu + 1 / 2
    L <- alpha - floor(alpha)
    const <- gamma(nu)/(gamma(alpha)*(4*pi)^(1/2)*kappa^(2*nu))
    p <- ceiling(alpha)
    B <- matrix(0, nrow = p + 1, ncol = p + 1)
    c <- numeric(p + 1)
    for (i in 0:p) {
        c[i + 1] <- 2 * kappa_integral(i, -alpha + 2 * p + 1 + L, kappa)
        for (j in 0:p) {
            B[j + 1, i + 1] <- 2 * kappa_integral(i + j, 2 * p + 1 + L, kappa)
        }
    }
    b <- solve(solve(diag(sqrt(diag(B))), B), solve(diag(sqrt(diag(B))), c))
    return (const*b)
}

lindgren_cov <- function(x, kappa, beta) {
    b <- lapprox(beta, kappa)
    b <- rev(b)
    rpk <- residue(1, b)
    r = rev(rpk$r)
    p = rev(rpk$p)
    length_r <- length(r)
    nu <- 2 * beta - 1 / 2
    C <- 0
    for (i in 1:length_r) {
        a <- sqrt(-p[i])
        Ci <- r[i]  * exp(-a * x) /(2 * a) 
        C <- C + Ci
    }
    C <- Re(C)
    if (is.nan(C[1])) {
        C[1] <- C[2]
    }
    return(C)
}