library(rSPDE)
source("aux_functions/aux_functions_cov.R")
rm(list=ls())

change.of.variables <- function(alpha,n, m, A) {
    
    #change of variables
    # u0 max(floor(alpha)-1,0) derivatives
    # ui floor(alpha) derivatives
    #alpha < 1 : 0, 0
    #alpha < 2 : 0, 1
    #alpha < 3 : 1, 2
    k1 = 1+max(c(floor(alpha)-1,0))
    k = floor(alpha)+1
    if(alpha < 1) { 
        B1 <- Diagonal(n)
    } else {
        B1 <- sparseMatrix(x=rep(1,k1*n),
                           i=(1:(k*n))[-seq(from=0,to=k*n,by=k)[-1]],
                           j=1:(k1*n),dims=c(k*n,k1*n)) #not sure if correct, map u0 to u1
    }
    
    I1 = Diagonal(k1*n)
    I2 = Diagonal(k*n)
    O1 <- Matrix(0,k*n,k1*n)
    O2 <- Matrix(0,k*n,k*n)
    if(m==1) {
        B <- rbind(cbind(I1, Matrix(0,k1*n,k*n)),
                   cbind(-B1, I2))    
    } else if(m==2) {
        B <- rbind(cbind(I1, Matrix(0,k1*n,k*n*m)),
                   cbind(-B1, I2,O2),
                   cbind(O1,-I2, I2))  
    } else if(m == 3) {
        B <- rbind(cbind(I1, Matrix(0,k1*n,k*n*m)),
                   cbind(-B1, I2,O2,O2),
                   cbind(O1,-I2, I2, O2),
                   cbind(O1,O2,-I2, I2))  
    } else if(m == 4) {
        B <- rbind(cbind( I1, Matrix(0,k1*n,k*n*m)),
                   cbind(-B1, I2, O2, O2, O2),
                   cbind( O1,-I2, I2, O2, O2),
                   cbind( O1, O2,-I2, I2, O2),
                   cbind( O1, O2, O2,-I2, I2))  
    } else if(m == 5) {
        B <- rbind(cbind( I1, Matrix(0,k1*n,k*n*m)),
                   cbind(-B1, I2, O2, O2, O2, O2),
                   cbind( O1,-I2, I2, O2, O2, O2),
                   cbind( O1, O2,-I2, I2, O2, O2),
                   cbind( O1, O2, O2,-I2, I2, O2),
                   cbind( O1, O2, O2, O2, -I2, I2))  
    } else if(m==6) {
        B <- rbind(cbind( I1, Matrix(0,k1*n,k*n*m)),
                   cbind(-B1, I2, O2, O2, O2, O2, O2),
                   cbind( O1,-I2, I2, O2, O2, O2, O2),
                   cbind( O1, O2,-I2, I2, O2, O2, O2),
                   cbind( O1, O2, O2,-I2, I2, O2, O2),
                   cbind( O1, O2, O2, O2,-I2, I2, O2),
                   cbind( O1, O2, O2, O2, O2,-I2, I2))  
    } else if(m==7){
        B <- rbind(cbind( I1, Matrix(0,k1*n,k*n*m)),
                   cbind(-B1, I2, O2, O2, O2, O2, O2, O2),
                   cbind( O1,-I2, I2, O2, O2, O2, O2, O2),
                   cbind( O1, O2,-I2, I2, O2, O2, O2, O2),
                   cbind( O1, O2, O2,-I2, I2, O2, O2, O2),
                   cbind( O1, O2, O2, O2,-I2, I2, O2, O2),
                   cbind( O1, O2, O2, O2, O2,-I2, I2, O2),
                   cbind( O1, O2, O2, O2, O2, O2,-I2, I2)) 
    } else if(m==8) {
        B <- rbind(cbind( I1, Matrix(0,k1*n,k*n*m)),
                   cbind(-B1, I2, O2, O2, O2, O2, O2, O2, O2),
                   cbind( O1,-I2, I2, O2, O2, O2, O2, O2, O2),
                   cbind( O1, O2,-I2, I2, O2, O2, O2, O2, O2),
                   cbind( O1, O2, O2,-I2, I2, O2, O2, O2, O2),
                   cbind( O1, O2, O2, O2,-I2, I2, O2, O2, O2),
                   cbind( O1, O2, O2, O2, O2,-I2, I2, O2, O2),
                   cbind( O1, O2, O2, O2, O2, O2,-I2, I2, O2),
                   cbind( O1, O2, O2, O2, O2, O2, O2,-I2, I2)) 
    } else {
        stop("only m<9 implemented")
    }
    A = cbind(Matrix(0, n.obs,k1*n + k*(m-1)*n), A[,(k1*n+k*(m-1)*n+1):(k1*n+k*m*n)])
    return(list(B = B, A = A))
}

range = 2
sigma = 2
sigma.e <- 0.2
n <- 5000
n.obs <- 5000
loc <- seq(0,n/100,length.out=n)
D <- as.matrix(dist(loc))
m=1

obs.ind <- sort(sample(1:n)[1:n.obs])

nu = 1.4
alpha <- nu + 1/2
kappa = sqrt(8*nu)/range
Sigma <- matern.covariance(h=D,kappa=kappa,nu=nu,sigma=sigma)
R <- chol(Sigma[obs.ind,obs.ind])
X <- t(R)%*%rnorm(n.obs)
Y <- X + sigma.e*rnorm(n.obs)


m = 6
t1 <- Sys.time()
Qrat <- matern.rational.ldl(loc = loc, order = m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp =  "spline", equally_spaced = FALSE)    
Q <- t(Qrat$L)%*%Qrat$D%*%Qrat$L
t2 <- Sys.time()
times.con <- t2 - t1


t1 <- Sys.time()
Qhat.rat <- Q + t(Qrat$A[obs.ind,])%*%Qrat$A[obs.ind,]/sigma.e^2        
mu.rat1 <- Qrat$A%*%solve(Qhat.rat, t(Qrat$A[obs.ind,])%*%Y/sigma.e^2)
t2 <- Sys.time()
times0 <- t2 - t1

t3 <- Sys.time()
tmp <- change.of.variables(alpha,n, m, Qrat$A)
A2 <- tmp$A
Q2 <- t(tmp$B)%*%Q%*%tmp$B
A2o <- A2[obs.ind,]
t4 <- Sys.time()
Qhat.rat2 <- Q2 + t(A2o)%*%A2o/sigma.e^2        
mu.rat2 <- A2%*%solve(Qhat.rat2, t(A2o)%*%Y/sigma.e^2)
t2 <- Sys.time()
times2 <- t4-t3

res <- data.frame(solve.orig = times0,
                  solve.new =  times2,
                  Q.build = times.con,
                  diff.est = max(abs((mu.rat1 - mu.rat2)/mu.rat1)),
                  alpha = alpha,
                  m = m)
print(res)


