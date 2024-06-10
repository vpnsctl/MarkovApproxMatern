library(rSPDE)
source("aux_functions/aux_functions_cov.R")
range = 2
sigma = 2
sigma.e <- 0.2
n <- 3000
n.obs <- 2500
n.rep <- 1
loc <- seq(0,n/100,length.out=n)

D <- as.matrix(dist(loc))
nu.vec <- seq(from = 0.1, to = 2.45, length.out = 10)
m.vec <- 2:6
t.nn <- t.rat <- t.nn1 <- t.rat1 <- err.nn <- err.rat <- matrix(0,nrow= length(nu.vec), ncol = length(m.vec))
for(kk in 1:n.rep) {
    obs.ind <- sort(sample(1:n)[1:n.obs])
    for(i in 1:length(nu.vec)) {
        cat(i/length(nu.vec),"\n")
        nu <- nu.vec[i]    
        alpha <- nu + 1/2
        kappa = sqrt(8*nu)/range
        Sigma <- matern.covariance(h=D,kappa=kappa,nu=nu,sigma=sigma)
        R <- chol(Sigma[obs.ind,obs.ind])
        X <- t(R)%*%rnorm(n.obs)
        Y <- X + sigma.e*rnorm(n.obs)
        Sigma.hat <- Sigma[obs.ind,obs.ind] + sigma.e^2*diag(n.obs)
        mu <- Sigma[,obs.ind]%*%solve(Sigma.hat,Y)
        
        for(j in 1:length(m.vec)) { 
            m <- m.vec[j]
            
            #mn <- round(sqrt((m+1)^3*(ceil(alpha)+1)^2))
            if(alpha<1) {
                mn <- m
            } else if (alpha < 2) {
                mn <- round(m*(ceil(alpha)+1)^1.5)    
            } else {
                mn <- round(m*(ceil(alpha)+1)^1.5)
            }
            
            
            t1 <- Sys.time()
            Qnn <- get.nnQ(loc = loc[obs.ind],kappa = kappa,nu = nu,sigma = sigma, n.nbr = mn)
            t2 <- Sys.time()
            Qhat <- Qnn + Diagonal(n.obs)/sigma.e^2        
            mu.nn <- solve(Qhat, Y/sigma.e^2)
            Bp <- get.nn.pred(loc = loc, kappa = kappa, nu = nu, sigma = sigma, n.nbr = mn, S = obs.ind)
            mu.nn <- Bp%*%mu.nn
            t3 <- Sys.time()
            t.nn[i,j] <- t.nn[i,j] + (t3-t2)/n.rep
            t.nn1[i,j] <- t.nn1[i,j] + (t2-t1)/n.rep
            
            t1 <- Sys.time()
            if(nu < 0.5) {
                Qrat <- matern.rational.ldl(loc = loc, order = m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp =  "spline")    
            } else {
                Qrat <- matern.rational.ldl(loc = loc, order = m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp =  "spline")    
            }
            
            Q <- t(Qrat$L)%*%Qrat$D%*%Qrat$L
            t2 <- Sys.time()
            Qhat.rat <- Q + t(Qrat$A[obs.ind,])%*%Qrat$A[obs.ind,]/sigma.e^2        
            mu.rat <- Qrat$A%*%solve(Qhat.rat, t(Qrat$A[obs.ind,])%*%Y/sigma.e^2)
            t3 <- Sys.time()
            t.rat[i,j] <- t.rat[i,j] + (t3-t2)/n.rep
            t.rat1[i,j] <- t.rat1[i,j] + (t3-t1)/n.rep
            err.nn[i,j] <- err.nn[i,j] + sqrt((loc[2]-loc[1])*sum((mu-mu.nn)^2))/n.rep
            err.rat[i,j] <- err.rat[i,j] + sqrt((loc[2]-loc[1])*sum((mu-mu.rat)^2))/n.rep
        }
    }
    
}
    

par(mfrow=c(2,2))
plot(nu.vec,err.rat[,1],type="l", log= "y", ylim = c(min(c(min(err.rat),min(err.nn))),
                                                     max(c(max(err.rat),max(err.nn)))),
     main = "error")
lines(nu.vec,err.rat[,2], col = 2)
lines(nu.vec,err.rat[,3], col = 3)
lines(nu.vec,err.rat[,4], col = 4)
lines(nu.vec,err.rat[,5], col = 5)
lines(nu.vec,err.nn[,1], col = 1,lty=2)
lines(nu.vec,err.nn[,2], col = 2,lty=2)
lines(nu.vec,err.nn[,3], col = 3,lty=2)
lines(nu.vec,err.nn[,4], col = 4,lty=2)
lines(nu.vec,err.nn[,5], col = 5,lty=2)

plot(nu.vec,t.rat[,1],type="l", log= "y", ylim = c(min(c(min(t.rat),min(t.nn))),
                                                     max(c(max(t.rat),max(t.nn)))),
     main = "time solve")
lines(nu.vec,t.rat[,2], col = 2)
lines(nu.vec,t.rat[,3], col = 3)
lines(nu.vec,t.rat[,4], col = 4)
lines(nu.vec,t.rat[,5], col = 5)
lines(nu.vec,t.nn[,1], col = 1,lty=2)
lines(nu.vec,t.nn[,2], col = 2,lty=2)
lines(nu.vec,t.nn[,3], col = 3,lty=2)
lines(nu.vec,t.nn[,4], col = 4,lty=2)
lines(nu.vec,t.nn[,5], col = 5,lty=2)

plot(nu.vec,t.rat1[,1],type="l", log= "y", ylim = c(min(c(min(t.rat1),min(t.nn1))),
                                                   max(c(max(t.rat1),max(t.nn1)))),
     main = "time construct")
lines(nu.vec,t.rat1[,2], col = 2)
lines(nu.vec,t.rat1[,3], col = 3)
lines(nu.vec,t.rat1[,4], col = 4)
lines(nu.vec,t.rat1[,5], col = 5)
lines(nu.vec,t.nn1[,1], col = 1,lty=2)
lines(nu.vec,t.nn1[,2], col = 2,lty=2)
lines(nu.vec,t.nn1[,3], col = 3,lty=2)
lines(nu.vec,t.nn1[,4], col = 4,lty=2)
lines(nu.vec,t.nn1[,5], col = 5,lty=2)

plot(nu.vec,t.rat1[,1]+t.rat[,1],type="l", log= "y", ylim = c(min(c(min(t.rat1+t.rat),min(t.nn1+t.nn))),
                                                   max(c(max(t.rat1+t.rat),max(t.nn1+t.nn)))),
     main = "time total")
lines(nu.vec,t.rat1[,2]+t.rat[,2], col = 2)
lines(nu.vec,t.rat1[,3]+t.rat[,3], col = 3)
lines(nu.vec,t.rat1[,4]+t.rat[,4], col = 4)
lines(nu.vec,t.rat1[,5]+t.rat[,5], col = 5)
lines(nu.vec,t.nn1[,1]+t.nn[,1], col = 1,lty=2)
lines(nu.vec,t.nn1[,2]+t.nn[,2], col = 2,lty=2)
lines(nu.vec,t.nn1[,3]+t.nn[,3], col = 3,lty=2)
lines(nu.vec,t.nn1[,4]+t.nn[,4], col = 4,lty=2)
lines(nu.vec,t.nn1[,5]+t.nn[,5], col = 5,lty=2)

