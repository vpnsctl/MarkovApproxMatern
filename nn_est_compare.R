library(rSPDE)
source("aux_functions/aux_functions_cov.R")

range = 5
sigma = 2
sigma.e <- 0.2
n <- 1000
n.obs <- 200
n.rep <- 1
loc <- seq(0,n/100,length.out=n)
obs.ind <- sort(sample(1:n)[1:n.obs])
D <- as.matrix(dist(loc))

nu <- 1.4
m = 6

est.nu <- FALSE


alpha <- nu + 1/2
kappa = sqrt(8*nu)/range

#mn <- round(sqrt((m+1)^3*(ceil(alpha)+1)^2))
#mn <- round(((m+2)*(ceil(alpha)+1)^(1.5)))
mn <- round(((m+2)*(ceil(alpha)+1)^(1.5)))

Sigma <- matern.covariance(h=D,kappa=kappa,nu=nu,sigma=sigma)
R <- chol(Sigma[obs.ind,obs.ind])
X <- t(R)%*%rnorm(n.obs)
Y <- X + sigma.e*rnorm(n.obs)
Sigma.hat <- Sigma[obs.ind,obs.ind] + sigma.e^2*diag(n.obs)
mu <- Sigma[,obs.ind]%*%solve(Sigma.hat,Y)


t1 <- Sys.time()
if(est.nu) {
    theta0 <- c(log(kappa), log(sigma), log(nu), log(sigma.e))
    par <- optim(theta0, rat.like, method = "L-BFGS-B", loc = loc[obs.ind], Y = Y, m = m, nu = nu)
    kappa_est = exp(par$par[1])
    sigma_est = exp(par$par[2])
    nu_est = exp(par$par[3])
    sigma.e_est = exp(par$par[4])
    Qrat <- matern.rational.ldl(loc = loc, order = m, nu = nu_est, kappa = kappa_est, 
                                sigma = sigma_est, type_rational = "brasil", type_interp =  "spline")    
    
} else {
    theta0 <- c(log(kappa), log(sigma), log(sigma.e))
    par <- optim(theta0, rat.like, method = "L-BFGS-B", loc = loc[obs.ind], Y = Y, m = m, nu = nu)
    kappa_est = exp(par$par[1])
    sigma_est = exp(par$par[2])
    sigma.e_est = exp(par$par[3])
    Qrat <- matern.rational.ldl(loc = loc, order = m, nu = nu, kappa = kappa_est, 
                                sigma = sigma_est, type_rational = "brasil", type_interp =  "spline")    
    
}

Q <- t(Qrat$L)%*%Qrat$D%*%Qrat$L

Qhat.rat <- Q + t(Qrat$A[obs.ind,])%*%Qrat$A[obs.ind,]/sigma.e_est^2        
mu.rat <- Qrat$A%*%solve(Qhat.rat, t(Qrat$A[obs.ind,])%*%Y/sigma.e_est^2)
t2 <- Sys.time()
t.rat <- t2 - t1

t1 <- Sys.time()
if(est.nu) {
    par2 <- optim(theta0, nn.like, method = "L-BFGS-B", loc = loc[obs.ind], Y = Y, n.nbr = mn, nu = nu)
    kappa_est = exp(par2$par[1])
    sigma_est = exp(par2$par[2])
    nu_est = exp(par2$par[3])
    sigma.e_est = exp(par2$par[4])
    Qnn <- get.nnQ(loc = loc,kappa = kappa_est,nu = nu_est,sigma = sigma_est, n.nbr = mn, S = obs.ind)
} else {
    par2 <- optim(theta0, nn.like, method = "L-BFGS-B", loc = loc[obs.ind], Y = Y, n.nbr = mn, nu = nu)
    kappa_est = exp(par2$par[1])
    sigma_est = exp(par2$par[2])
    sigma.e_est = exp(par2$par[3])
    Qnn <- get.nnQ(loc = loc,kappa = kappa_est,nu = nu,sigma = sigma_est, n.nbr = mn, S = obs.ind)
    
}
Qnn <- get.nnQ(loc = loc[obs.ind],kappa = kappa_est,nu = nu,sigma = sigma_est, n.nbr = mn)

A = Diagonal(n.obs)
Qhat <- Qnn + t(A)%*%A/sigma.e^2        
mu.nn <- solve(Qhat, t(A)%*%Y/sigma.e^2)
Bp <- get.nn.pred(loc = loc, kappa = kappa_est, nu = nu, sigma = sigma_est, n.nbr = mn, S = obs.ind)$B
mu.nn <- Bp%*%mu.nn
t2 <- Sys.time()
t.nn <- t2 - t1

res <- data.frame(error = c(sqrt((loc[2]-loc[1])*sum((mu-mu.nn)^2)), sqrt((loc[2]-loc[1])*sum((mu-mu.rat)^2))),
                  time = c(t.nn, t.rat),
                  count = c(par2$counts[1], par$counts[1]),
                  row.names = c("NN", "Rat"))
print(res)
par(mfrow = c(1,1))
plot(loc,mu, type = "l", main = "black true, red nn, green markov")
points(loc[obs.ind],Y)
lines(loc,mu.nn,col=2)
lines(loc,mu.rat,col=3)