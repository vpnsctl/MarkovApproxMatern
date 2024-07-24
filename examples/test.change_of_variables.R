rm(list=ls())
library(Matrix)
range = 2
sigma = 2
sigma.e <- 0.2
n <- 10000
n.obs <- 10000
loc <- seq(0,n/100,length.out=n)

obs.ind <- sort(sample(1:n)[1:n.obs])

nu = 2.4
alpha <- nu + 1/2
kappa = sqrt(8*nu)/range
#Sigma <- matern.covariance(h=D,kappa=kappa,nu=nu,sigma=sigma)
#R <- chol(Sigma[obs.ind,obs.ind])
#X <- t(R)%*%rnorm(n.obs)
#Y <- X + sigma.e*rnorm(n.obs)
Y <- matrix(rnorm(n.obs),n.obs,1)


m = 1
t1 <- Sys.time()
Qrat <- rSPDE:::matern.rational.precision(loc = loc, order = m, nu = nu, kappa = kappa, 
                                          sigma = sigma, type_rational = "brasil", 
                                          type_interp =  "spline")    
Q <- Qrat$Q
t2 <- Sys.time()
times.con1 <- t2 - t1

t1 <- Sys.time()
Qhat.rat <- Q + t(Qrat$A[obs.ind,])%*%Qrat$A[obs.ind,]/sigma.e^2        
mu.rat1 <- Qrat$A%*%solve(Qhat.rat, t(Qrat$A[obs.ind,])%*%Y/sigma.e^2)
t2 <- Sys.time()
times1 <- t2 - t1

t1 <- Sys.time()
Qrat <- rSPDE:::matern.rational.precision(loc = loc, order = m, nu = nu, kappa = kappa,
                                          cumsum = TRUE,
                                          sigma = sigma, type_rational = "brasil", 
                                          type_interp =  "spline")    
Q2 <- Qrat$Q
A2 <- Qrat$A
t2 <- Sys.time()
times.con2 <- t2 - t1

t3 <- Sys.time()
A2o <- A2[obs.ind,]
t4 <- Sys.time()
Qhat.rat2 <- Q2 + t(A2o)%*%A2o/sigma.e^2        
mu.rat2 <- A2%*%solve(Qhat.rat2, t(A2o)%*%Y/sigma.e^2)
t2 <- Sys.time()
times2 <- t4-t3

res <- data.frame(Q.build = c(times.con1,times.con2),
                  solve = c(times1,times2),
                  total = c(times.con1+times1, times.con2+times2),
                  row.names = c("Original", "New"))
print(res)
cat(max(abs(mu.rat1-mu.rat2)))


