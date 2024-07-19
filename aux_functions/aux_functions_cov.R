# rSPDE functions and coeffs

library(rSPDE)

# Sparse matrices

library(Matrix)

# Toeplitz

library(SuperGauss)

# Our method

rat.like <- function(theta, loc, m, Y,nu = NULL) {
    
    kappa = exp(theta[1])
    sigma = exp(theta[2])
    if(length(theta)==4) { 
        nu = exp(theta[3])
        sigma.e = exp(theta[4])
    } else {
        sigma.e = exp(theta[3])    
    }
    tmp <- matern.rational.ldl(loc = loc, order = m, nu = nu, kappa = kappa, sigma = sigma, 
                               type_rational = "brasil", type_interp =  "spline")    
    n <- length(Y)
    
    prior.ld <- 0.5 * sum(log(diag(tmp$D)))
    
    Q.post <- t(tmp$L)%*%tmp$D%*%tmp$L + t(tmp$A) %*% tmp$A / sigma.e^2
    
    R.post <- Matrix::Cholesky(Q.post)
    
    posterior.ld <- c(determinant(R.post, logarithm = TRUE, sqrt = TRUE)$modulus)
    
    AtY <- t(tmp$A) %*% Y / sigma.e^2
    
    mu.post <- solve(R.post, AtY, system = "A")
    
    lik <- prior.ld - posterior.ld - n * log(sigma.e)
    v <- tmp$L %*% mu.post
    lik <- lik -  0.5 * (t(v) %*% tmp$D %*% v +
                             t(Y - tmp$A %*% mu.post) %*% (Y - tmp$A %*% mu.post) / sigma.e^2)
    
    return(-as.double(lik))
}

# get corresponding m to other method based on ours:

get_m <- function(nu, m, method = c("nngp", "statespace", "kl"), type = c("prediction", "simulation")){
  type <- type[[1]]
  alpha <- nu + 0.5
  if(type == "prediction"){
    if(method == "statespace"){
      return(pmax((m-1)*ceil(alpha) + 1,1))
    } else if(method == "nngp"){
        if(alpha < 1.5){
            return(m+1)
        }
      return(pmax(2*round((m+3)*sqrt(m) * ceil(alpha)^(3/2)),1))
    } else if(method == "kl"){
      return(pmax(9*round((m+1)*sqrt(m) * ceil(alpha)^(3/2)),2))
    } else{
      stop("method not implemented.")
    }
  } else{
    if(method == "statespace"){
    #   return(pmax(round(ceil(alpha) * (m^(1/3)-1)) + 1,1))
    return(m+floor(alpha) + 1)
    } else if(method == "nngp"){
      m_nngp <- round(8*(m+1)* ceil(alpha+1))
      d_m <- diff(m_nngp)
      d_m <- ifelse(d_m == 0, 1, d_m)
      m_nngp <- c(m_nngp[1], m_nngp[1]+cumsum(d_m))
      return(pmax(m_nngp,1))
    } else if(method == "kl"){
        m_kl <- round(8*(m+1) * ceil(alpha+1))
        d_m <- diff(m_kl)
        d_m <- ifelse(d_m == 0, 1, d_m)
        m_kl <- c(m_kl[1], m_kl[1]+cumsum(d_m))
        m_kl <- m_kl + 1
        return(pmax(m_kl,1))
    } else{
      stop("method not implemented.")
    }
  } 
}

## NN GP - nearest neighbor

get.neighbor <- function(i, n.nbr, S = NULL) {
    if(is.null(S)) {
        S = 1:i
        S.i <- sort(S[S<i], decreasing = TRUE)
        if(length(S.i) < n.nbr) {
            return(sort(S.i))
        } else {
            return(sort(S.i[1:n.nbr]))
        }    
    } else {
        if(i %in% S) { 
            S.i <- sort(S[S<i], decreasing = TRUE)
            if(length(S.i) < n.nbr) {
                return(sort(S.i))
            } else {
                return(sort(S.i[1:n.nbr]))
            }    
        } else {
             dists <- abs(i - S)
             return(sort(sort(dists, index.return = TRUE)$ix[1:n.nbr]))
        }
    }
    
    
}

get.nn <- function(loc,kappa,nu,sigma, n.nbr, S = NULL) {
    k <- length(loc)
    ii <- numeric((2*k-n.nbr)*(n.nbr+1)/2)
    jj <- numeric((2*k-n.nbr)*(n.nbr+1)/2)
    val <- numeric((2*k-n.nbr)*(n.nbr+1)/2)
    
    Fs.d <- Fsi.d <- numeric(k)
    
    Fs.d[1] <- sigma^2
    Fsi.d[1] <- 1/sigma^2
    val[1] <- 1
    ii[1] <- 1
    jj[1] <- 1
    counter <- 1
    for(i in 2:k) {
        nbrs <- get.neighbor(i,n.nbr, S)
        if(isempty(nbrs)) {
            Fs.d[i] <- sigma^2
            Fsi.d[i] <- 1/sigma^2
            val[counter +1] <- 1
            ii[counter+1] <- i
            jj[counter+1] <- i
            counter = counter + 1
        } else {
            Sigma.in <- matern.covariance(h = abs(loc[nbrs]-loc[i]),
                                          kappa = kappa, nu = nu, sigma = sigma)
            Sigma.nn <- matern.covariance(h = as.matrix(dist(loc[nbrs])),
                                          kappa = kappa, nu = nu, sigma = sigma)
            N <- length(nbrs) + 1
            tmp <- solve(Sigma.nn, Sigma.in)
            val[counter + (1:N)] <- c(-t(tmp),1)
            ii[counter + (1:N)] <- rep(i,N)
            jj[counter + (1:N)] <- c(nbrs,i)
            Fs.d[i] <- sigma^2 - t(Sigma.in)%*%tmp
            Fsi.d[i] <- 1/Fs.d[i]    
            counter <- counter + N
        }
        
        
        
        
    }
    ii <- ii[1:counter]
    jj <- jj[1:counter]
    val <- val[1:counter]
    Bs <-  Matrix::sparseMatrix(i   = ii,
                                j    = jj,
                                x    = val,
                                dims = c(k, k))
    Fs <-  Matrix::Diagonal(k,Fs.d)
    Fsi <-  Matrix::Diagonal(k,Fsi.d)
    return(list(Bs = Bs, Fs = Fs, Fi = Fsi))
}

get.nnQ <- function(loc,kappa,nu,sigma, n.nbr,S=NULL) {
    tmp <- get.nn(loc,kappa,nu,sigma,n.nbr,S)
    
    return(t(tmp$Bs) %*% tmp$Fi %*%tmp$Bs)
}


get.nn.pred <- function(loc,kappa,nu,sigma, n.nbr, S = NULL) {
    
    n.S <- length(S)
    n.loc <- length(loc)
    k <- n.loc - n.S
    N <- k*n.nbr + n.S
    ii <- numeric(N)
    jj <- numeric(N)
    val <- numeric(N)
    Fs.d <- numeric(n.loc)
    counter <- 0
    ii[1:n.S] <- S
    jj[1:n.S] <- 1:n.S
    val[1:n.S] <- 1
    counter <- n.S
    not.observed <- setdiff(1:length(loc), S)
    reo <- c(S, not.observed)
    for(i in not.observed) {
            dists <- abs(loc[i] - loc[S])
            nbrs <- sort(sort(dists, index.return = TRUE)$ix[1:n.nbr])
            Sigma.in <- matern.covariance(h = abs(loc[S[nbrs]]-loc[i]),
                                          kappa = kappa, nu = nu, sigma = sigma)
            Sigma.nn <- matern.covariance(h = as.matrix(dist(loc[S[nbrs]])),
                                          kappa = kappa, nu = nu, sigma = sigma)
            
            tmp <- solve(Sigma.nn, Sigma.in)
            val[counter + (1:n.nbr)] <- t(tmp)
            ii[counter + (1:n.nbr)] <- rep(i,n.nbr)
            jj[counter + (1:n.nbr)] <- nbrs
            counter <- counter + n.nbr    
            Fs.d[i] <- sigma^2 - t(Sigma.in)%*%tmp
    }
    Bs <-  Matrix::sparseMatrix(i   = ii,
                                j    = jj,
                                x    = val,
                                dims = c(n.loc, n.S))
    Fs <-  Matrix::Diagonal(n.loc,Fs.d)
    return(list(B = Bs, F = Fs))
}


nn.like <- function(theta, loc, n.nbr, Y,nu = NULL) {
    kappa = exp(theta[1])
    sigma = exp(theta[2])
    if(length(theta)==4) { 
        nu = exp(theta[3])
        sigma.e = exp(theta[4])
    } else {
        sigma.e = exp(theta[3])    
    }

    tmp <- get.nn(loc,kappa,nu,sigma,n.nbr)
    
    n <- length(Y)
    
    prior.ld <- 0.5 * sum(log(diag(tmp$Fi)))
    
    Q.post <- t(tmp$Bs)%*%tmp$Fi%*%tmp$Bs + Diagonal(n) / sigma.e^2
    
    
    R.post <- Matrix::Cholesky(Q.post)
    
    posterior.ld <-  c(determinant(R.post, logarithm = TRUE, sqrt = TRUE)$modulus)
    
    AtY <- Y / sigma.e^2
    
    mu.post <- solve(R.post, AtY, system = "A")
    
    lik <- prior.ld - posterior.ld - n * log(sigma.e)
    v <- tmp$Bs %*% mu.post
    lik <- lik - 0.5 * (t(v) %*% tmp$Fi %*% v +
                            t(Y - mu.post) %*% (Y - mu.post) / sigma.e^2)
    
    return(-as.double(lik))
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
