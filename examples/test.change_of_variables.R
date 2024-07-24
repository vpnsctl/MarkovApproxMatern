m = 8
t1 <- Sys.time()
Qrat <- rSPDE:::matern.rational.precision(loc = loc, order = m, nu = nu, kappa = kappa, 
                                          sigma = sigma, type_rational = "brasil", 
                                          type_interp =  "spline")    
Q <- Qrat$Q
t2 <- Sys.time()
times.con1 <- t2 - t1

t3 <- Sys.time()
Qhat.rat <- Q + t(Qrat$A[obs.ind,])%*%Qrat$A[obs.ind,]/sigma.e^2        
mu.rat1 <- Qrat$A%*%solve(Qhat.rat, t(Qrat$A[obs.ind,])%*%Y/sigma.e^2)
t4 <- Sys.time()
times1 <- t4 - t3

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
Qhat.rat2 <- Q2 + t(A2o)%*%A2o/sigma.e^2        
mu.rat2 <- A2%*%solve(Qhat.rat2, t(A2o)%*%Y/sigma.e^2)
t4 <- Sys.time()
times2 <- t4-t3

res <- data.frame(Q.build = c(times.con1,times.con2),
                  solve = c(times1,times2),
                  total = c(times.con1+times1, times.con2+times2),
                  row.names = c("Original", "New"))
print(res)
cat(max(abs(mu.rat1-mu.rat2)))


