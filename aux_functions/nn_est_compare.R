library(rSPDE)
source("aux_functions/aux_functions_cov.R")

nn.like <- function(theta, loc, n.nbr, Y) {
    kappa = exp(theta[1])
    sigma = exp(theta[2])
    nu = exp(theta[3])
    sigma.e = exp(theta[4])
    
    cat(kappa,sigma,nu,sigma.e,"\n")
    
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

rat.like <- function(theta, loc, m, Y) {
    
    kappa = exp(theta[1])
    sigma = exp(theta[2])
    nu = exp(theta[3])
    sigma.e = exp(theta[4])
    cat(kappa,sigma,nu,sigma.e,"\n")
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

range = 1
sigma = 2
sigma.e <- 0.2
n <- 1000
n.rep <- 1
loc <- seq(0,n/100,length.out=n)
D <- as.matrix(dist(loc))

nu <- 1.6
m = 3

alpha <- nu + 1/2
kappa = sqrt(8*nu)/range
Sigma <- matern.covariance(h=D,kappa=kappa,nu=nu,sigma=sigma)
R <- chol(Sigma)
X <- t(R)%*%rnorm(n)
Y <- X + sigma.e*rnorm(n)
Sigma.hat <- Sigma + sigma.e^2*diag(n)
mu <- Sigma%*%solve(Sigma.hat,Y)
        
theta0 <- c(log(kappa), log(sigma), log(nu), log(sigma.e))
par <- optim(theta0, rat.like, method = "L-BFGS-B", loc = loc, Y = Y, m = m)
kappa_est = exp(par$par[1])
sigma_est = exp(par$par[2])
nu_est = exp(par$par[3])
sigma.e_est = exp(par$par[4])
Qrat <- matern.rational.ldl(loc = loc, order = m, nu = nu_est, kappa = kappa_est, 
                            sigma = sigma_est, type_rational = "brasil", type_interp =  "spline")    

Q <- t(Qrat$L)%*%Qrat$D%*%Qrat$L

Qhat.rat <- Q + t(Qrat$A)%*%Qrat$A/sigma.e_est^2        
mu.rat <- Qrat$A%*%solve(Qhat.rat, t(Qrat$A)%*%Y/sigma.e_est^2)

#mn <- round(sqrt((m+1)^3*(ceil(alpha)+1)^2))
mn <- round(((m+2)*(ceil(alpha)+1)^(1.5)))

theta0 <- c(log(kappa), log(sigma), log(nu), log(sigma.e))
par <- optim(theta0, nn.like, method = "L-BFGS-B", loc = loc, Y = Y, n.nbr = mn)
kappa_est = exp(par$par[1])
sigma_est = exp(par$par[2])
nu_est = exp(par$par[3])
sigma.e_est = exp(par$par[4])

Qnn <- get.nnQ(loc = loc,kappa = kappa_est,nu = nu_est,sigma = sigma_est, n.nbr = mn)
Qhat <- Qnn + Diagonal(n)/sigma.e^2        
mu.nn <- solve(Qhat, Y/sigma.e^2)

cat(sqrt((loc[2]-loc[1])*sum((mu-mu.nn)^2)), sqrt((loc[2]-loc[1])*sum((mu-mu.rat)^2)))
