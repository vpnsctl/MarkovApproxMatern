library(rSPDE)
library(fmesher)
source("aux_functions/aux_functions_cov.R")
range = 1
sigma = 2
sigma.e <- 0.1
n <- 5001
n.obs <- 5001
n.rep <- 5
loc <- seq(0,n/100,length.out=n)

D <- as.matrix(dist(loc))

nu.vec <- c(0.4,1.4,2.4)#seq(from = 0.1, to = 2.45, length.out = 10)

new.data <- TRUE
if(new.data){
    obs.ind.list <- list()
    mu.list <- list()
    Y.list <- list()
    for(kk in 1:n.rep) {
        obs.ind <- sort(sample(1:n)[1:n.obs])
        obs.ind.list[[kk]] = obs.ind
        Y.list[[kk]] = list()
        mu.list[[kk]] = list()
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
            mu.list[[kk]][[i]] = mu
            Y.list[[kk]][[i]] = Y
        }
    }
}

m.vec <- 1:5
t.nn <- t.rat <- t.rspde <- t.rspde1 <- t.nn1 <- t.rat1 <- err.nn <- err.rat <- err.rspde <- matrix(0,nrow= length(nu.vec), ncol = length(m.vec))
for(kk in 1:n.rep) {
    obs.ind <- obs.ind.list[[kk]]
    for(i in 1:length(nu.vec)) {
        nu <- nu.vec[i]    
        alpha <- nu + 1/2
        kappa = sqrt(8*nu)/range
        cat(i/length(nu.vec),"\n")
        mu <- mu.list[[kk]][[i]]
        Y <- Y.list[[kk]][[i]]
        for(j in 1:length(m.vec)) { 
            m <- m.vec[j]
            
            #rational
            t1 <- Sys.time()
            Qrat <-rSPDE:::matern.rational.ldl(loc = loc, order = m, nu = nu, kappa = kappa, 
                                               sigma = sigma, type_rational = "brasil", 
                                               type_interp =  "spline", equally_spaced = FALSE)    
            
            Q <- t(Qrat$L)%*%Qrat$D%*%Qrat$L
            t2 <- Sys.time()
            Qhat.rat <- Q + t(Qrat$A[obs.ind,])%*%Qrat$A[obs.ind,]/sigma.e^2        
            mu.rat <- Qrat$A%*%solve(Qhat.rat, t(Qrat$A[obs.ind,])%*%Y/sigma.e^2)
            t3 <- Sys.time()
            t.rat[i,j] <- t.rat[i,j] + (t3-t2)/n.rep
            t.rat1[i,j] <- t.rat1[i,j] + (t2-t1)/n.rep
            
            #rspde
            if(alpha<1) {
                mr <- m
                extension <- c(seq(from = 0, to = 4*range, length.out = 200)[-1],
                               seq(from = range, to = 4*range, length.out = 50)[-1])
                loc_mesh <- seq(0,n/100,length.out=(m+1)*(n-1)+1)
                loc_mesh <- c(-rev(extension), loc_mesh, loc[length(loc)] + extension)
            } else if (alpha < 2) {
                mr <- m
                extension <- c(seq(from = 0, to = 4*range, length.out = 200)[-1],
                               seq(from = range, to = 4*range, length.out = 50)[-1])
                if(m>3){
                    loc_mesh <- seq(0,n/100,length.out=2*(m+1)*(n-1)+1)    
                } else {
                    loc_mesh <- seq(0,n/100,length.out=(m+1)*(n-1)+1)    
                }
                
                loc_mesh <- c(-rev(extension), loc_mesh, loc[length(loc)] + extension)
            } else {
                mr <- m
                extension <- c(seq(from = 0, to = 4*range, length.out = 200)[-1],
                               seq(from = range, to = 4*range, length.out = 50)[-1])
                if(m>3){
                    loc_mesh <- seq(0,n/100,length.out=2*(m+1)*(n-1)+1)    
                } else {
                    loc_mesh <- seq(0,n/100,length.out=(m+1)*(n-1)+1)    
                }
                loc_mesh <- c(-rev(extension), loc_mesh, loc[length(loc)] + extension)
            }
            t1 <- Sys.time()
            
            mesh <- rSPDE::rSPDE.fem1d(loc_mesh)
            A  <- rSPDE::rSPDE.A1d(mesh, loc)
            
            op.cov <- matern.operators(sigma = sigma, range = range, nu = nu,
                                       loc_mesh = loc_mesh, d = 1, m = mr,
                                       parameterization = "matern", 
                                       type_rational_approximation = "brasil")
            t2 <- Sys.time()
            mu.rspde <- as.vector(predict(op.cov, A[obs.ind,], A, Y, sigma.e)$mean)
            t3 <- Sys.time()
            t.rspde[i,j] <- t.rspde[i,j] + (t3-t2)/n.rep
            t.rspde1[i,j] <- t.rspde1[i,j] + (t2-t1)/n.rep
            #err.nn[i,j] <- err.nn[i,j] + sqrt((loc[2]-loc[1])*sum((mu-mu.nn)^2))/n.rep
            err.rat[i,j] <- err.rat[i,j] + sqrt((loc[2]-loc[1])*sum((mu-mu.rat)^2))/n.rep
            err.rspde[i,j] <- err.rspde[i,j] + sqrt((loc[2]-loc[1])*sum((mu-mu.rspde)^2))/n.rep
        }
    }
    
}


par(mfrow=c(2,2))
plot(nu.vec,err.rat[,1],type="l", log= "y", ylim = c(min(c(min(err.rat),min(err.rspde))),
                                                     max(c(max(err.rat),max(err.rspde)))),
     main = "error")
lines(nu.vec,err.rat[,2], col = 2)
lines(nu.vec,err.rat[,3], col = 3)
lines(nu.vec,err.rat[,4], col = 4)
lines(nu.vec,err.rat[,5], col = 5)
lines(nu.vec,err.rspde[,1], col = 1,lty=3)
lines(nu.vec,err.rspde[,2], col = 2,lty=3)
lines(nu.vec,err.rspde[,3], col = 3,lty=3)
lines(nu.vec,err.rspde[,4], col = 4,lty=3)
lines(nu.vec,err.rspde[,5], col = 5,lty=3)

plot(nu.vec,t.rat[,1],type="l", log= "y", ylim = c(min(c(min(t.rat),min(t.rspde))),
                                                   max(c(max(t.rat),max(t.rspde)))),
     main = "time solve")
lines(nu.vec,t.rat[,2], col = 2)
lines(nu.vec,t.rat[,3], col = 3)
lines(nu.vec,t.rat[,4], col = 4)
lines(nu.vec,t.rat[,5], col = 5)
lines(nu.vec,t.rspde[,1], col = 1,lty=3)
lines(nu.vec,t.rspde[,2], col = 2,lty=3)
lines(nu.vec,t.rspde[,3], col = 3,lty=3)
lines(nu.vec,t.rspde[,4], col = 4,lty=3)
lines(nu.vec,t.rspde[,5], col = 5,lty=3)

plot(nu.vec,t.rat1[,1],type="l", log= "y", ylim = c(min(c(min(t.rat1),min(t.rspde1))),
                                                    max(c(max(t.rat1),max(t.rspde1)))),
     main = "time construct")
lines(nu.vec,t.rat1[,2], col = 2)
lines(nu.vec,t.rat1[,3], col = 3)
lines(nu.vec,t.rat1[,4], col = 4)
lines(nu.vec,t.rat1[,5], col = 5)
lines(nu.vec,t.rspde1[,1], col = 1,lty=3)
lines(nu.vec,t.rspde1[,2], col = 2,lty=3)
lines(nu.vec,t.rspde1[,3], col = 3,lty=3)
lines(nu.vec,t.rspde1[,4], col = 4,lty=3)
lines(nu.vec,t.rspde1[,5], col = 5,lty=3)

plot(nu.vec,t.rat1[,1]+t.rat[,1],type="l", log= "y", ylim = c(min(c(min(t.rat1+t.rat),min(t.rspde1+t.rspde))),
                                                              max(c(max(t.rat1+t.rat),max(t.rspde1+t.rspde)))),
     main = "time total")
lines(nu.vec,t.rat1[,2]+t.rat[,2], col = 2)
lines(nu.vec,t.rat1[,3]+t.rat[,3], col = 3)
lines(nu.vec,t.rat1[,4]+t.rat[,4], col = 4)
lines(nu.vec,t.rat1[,5]+t.rat[,5], col = 5)
#lines(nu.vec,t.nn1[,1]+t.nn[,1], col = 1,lty=2)
#lines(nu.vec,t.nn1[,2]+t.nn[,2], col = 2,lty=2)
#lines(nu.vec,t.nn1[,3]+t.nn[,3], col = 3,lty=2)
#lines(nu.vec,t.nn1[,4]+t.nn[,4], col = 4,lty=2)
#lines(nu.vec,t.nn1[,5]+t.nn[,5], col = 5,lty=2)
lines(nu.vec,t.rspde1[,1]+t.rspde[,1], col = 1,lty=3)
lines(nu.vec,t.rspde1[,2]+t.rspde[,2], col = 2,lty=3)
lines(nu.vec,t.rspde1[,3]+t.rspde[,3], col = 3,lty=3)
lines(nu.vec,t.rspde1[,4]+t.rspde[,4], col = 4,lty=3)
lines(nu.vec,t.rspde1[,5]+t.rspde[,5], col = 5,lty=3)

