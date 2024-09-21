rm(list=ls())
source("aux_functions/aux_functions_cov.R")
source("examples/error.computations.R")
library(rSPDE)
library(foreach)
library(doParallel)
library(doSNOW)

cores=14

cl <- makeCluster(cores[1]-1) 
registerDoSNOW(cl)

range = 0.5 # not relative range
sigma = 1
sigma.e <- 0.1


Dists <- as.matrix(dist(loc))

nu.vec <- seq(from = 0.01, to = 2.49, by = 0.01)
#nu.vec <- c(0.3,1.3,1.4)
nu.vec <- nu.vec[length(nu.vec):1]
m.vec <- 1:6

iterations <- length(nu.vec)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


### Full


res = foreach(i = 1:iterations, .options.snow = opts, .packages=c('Matrix', 'rSPDE', 'pracma')) %dopar% {
    res <- error.computations(range = range, sigma = sigma, sigma.e = sigma.e, 
    n = n, n.obs = n.obs, loc = loc, nu = nu.vec[i], m.vec = m.vec, Dists = Dists,
    n.rep = n.rep)
    return(res)
}

err.ss <- err.fourier <- err.nn <- err.rat <- err.pca <- nu <- NULL
for(i in 1:length(res)) {
    nu <- c(nu, res[[i]]$nu)
    err.ss <- rbind(err.ss, res[[i]]$err.ss)
    err.fourier <- rbind(err.ss, res[[i]]$err.fourier)
    err.nn <- rbind(err.ss, res[[i]]$err.nn)
    err.rat <- rbind(err.ss, res[[i]]$err.rat)
    err.pca <- rbind(err.ss, res[[i]]$err.pca)
    
}
nu.vec <- nu

min.err <- min(c(min(err.ss), min(err.fourier), min(err.nn), min(err.rat), min(err.pca)))
max.err <- max(c(max(err.ss), max(err.fourier), max(err.nn), max(err.rat), max(err.pca)))
par(mfrow=c(1,1))
plot(nu.vec,err.rat[,1],type="l", log= "y",
     main = "error", ylim =c(min.err, max.err))
lines(nu.vec,err.rat[,2], col = 2)
lines(nu.vec,err.rat[,3], col = 3)
lines(nu.vec,err.rat[,4], col = 4)
lines(nu.vec,err.rat[,5], col = 5)
lines(nu.vec,err.nn[,1], col = 1,lty=2)
lines(nu.vec,err.nn[,2], col = 2,lty=2)
lines(nu.vec,err.nn[,3], col = 3,lty=2)
lines(nu.vec,err.nn[,4], col = 4,lty=2)
lines(nu.vec,err.nn[,5], col = 5,lty=2)
lines(nu.vec,err.pca[,1], col = 1,lty=3)
lines(nu.vec,err.pca[,2], col = 2,lty=3)
lines(nu.vec,err.pca[,3], col = 3,lty=3)
lines(nu.vec,err.pca[,4], col = 4,lty=3)
lines(nu.vec,err.pca[,5], col = 5,lty=3)

lines(nu.vec,err.fourier[,1], col = 1,lty=4)
lines(nu.vec,err.fourier[,2], col = 2,lty=4)
lines(nu.vec,err.fourier[,3], col = 3,lty=4)
lines(nu.vec,err.fourier[,4], col = 4,lty=4)
lines(nu.vec,err.fourier[,5], col = 5,lty=4)

lines(nu.vec,err.ss[,1], col = 1,lty=5)
lines(nu.vec,err.ss[,2], col = 2,lty=5)
lines(nu.vec,err.ss[,3], col = 3,lty=5)
lines(nu.vec,err.ss[,4], col = 4,lty=5)
lines(nu.vec,err.ss[,5], col = 5,lty=5)




########## SETUP

rm(list=ls())
source("aux_functions/aux_functions_cov.R")
source("examples/error.computations.R")
library(rSPDE)
library(foreach)
library(doParallel)
library(doSNOW)

cores=12

cl <- makeCluster(cores[1]-1) 
registerDoSNOW(cl)

sigma = 1
sigma.e <- 0.1


nu.vec <- seq(from = 0.01, to = 2.49, by = 0.01)
nu.vec <- nu.vec[length(nu.vec):1]
m.vec <- 1:6

iterations <- length(nu.vec)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
folder_to_save <- getwd()



### NO PCA + No Fourier + No SS + n = n_obs (to use fast methods)

n <- 5000
n.obs <- 5000
n.rep <- 100
loc <- seq(0,n/100,length.out=n)

# range <- 0.2
range <- 0.5

res = foreach(i = 1:iterations, .options.snow = opts, .packages=c('Matrix', 'rSPDE', 'pracma')) %dopar% {
    res <- error.computations_nopca_nofourier_noss_n_equal_nobs(range = range, sigma = sigma, sigma.e = sigma.e, n = n, loc = loc, nu = nu.vec[i], m.vec = m.vec, n.rep = n.rep, folder_to_save = folder_to_save)
    return(res)
}

err.ss <- err.nn <- err.rat <- nu <- NULL
for(i in 1:length(res)) {
    nu <- c(nu, res[[i]]$nu)
    err.nn <- rbind(err.nn, res[[i]]$err.nn)
    err.rat <- rbind(err.rat, res[[i]]$err.rat)
    
}

res_5000_pred <- list(nu = nu, err.nn = err.nn, err.rat = err.rat)

saveRDS(res_5000_pred, paste0("pred_tables/res_5000_range",range,"_rat_nngp.RDS"))

## 


rm(list=ls())
source("aux_functions/aux_functions_cov.R")
source("examples/error.computations.R")
library(rSPDE)
library(foreach)
library(doParallel)
library(doSNOW)

cores=19

cl <- makeCluster(cores[1]-1) 
registerDoSNOW(cl)

sigma = 1
sigma.e <- 0.1


nu.vec <- seq(from = 0.01, to = 2.49, by = 0.01)
nu.vec <- nu.vec[length(nu.vec):1]
m.vec <- 1:6

iterations <- length(nu.vec)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
folder_to_save <- getwd()



n <- 5000
n.obs <- 5000
n.rep <- 100
loc <- seq(0,n/100,length.out=n)
# range <- 0.5
range <- 1

res = foreach(i = 1:iterations, .options.snow = opts, .packages=c('Matrix', 'rSPDE', 'pracma')) %dopar% {
    res <- error.computations_nopca_nofourier_noss_n_equal_nobs(range = range, sigma = sigma, sigma.e = sigma.e, n = n, loc = loc, nu = nu.vec[i], m.vec = m.vec, n.rep = n.rep, folder_to_save = folder_to_save)
    return(res)
}

err.ss <- err.nn <- err.rat <- nu <- NULL
for(i in 1:length(res)) {
    nu <- c(nu, res[[i]]$nu)
    err.nn <- rbind(err.nn, res[[i]]$err.nn)
    err.rat <- rbind(err.rat, res[[i]]$err.rat)
    
}

res_5000_pred <- list(nu = nu, err.nn = err.nn, err.rat = err.rat)

saveRDS(res_5000_pred, paste0("pred_tables/res_5000_range",range,"_rat_nngp.RDS"))



## 


rm(list=ls())
source("aux_functions/aux_functions_cov.R")
source("examples/error.computations.R")
library(rSPDE)
library(foreach)
library(doParallel)
library(doSNOW)

cores=22

cl <- makeCluster(cores[1]-1) 
registerDoSNOW(cl)

sigma = 1
sigma.e <- 0.1


nu.vec <- seq(from = 0.01, to = 2.49, by = 0.01)
nu.vec <- nu.vec[length(nu.vec):1]
m.vec <- 1:6

iterations <- length(nu.vec)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
folder_to_save <- getwd()


n <- 5000
n.obs <- 5000
n.rep <- 100
loc <- seq(0,n/100,length.out=n)
# range <- 1
range <- 2

res = foreach(i = 1:iterations, .options.snow = opts, .packages=c('Matrix', 'rSPDE', 'pracma')) %dopar% {
    res <- error.computations_nopca_nofourier_noss_n_equal_nobs(range = range, sigma = sigma, sigma.e = sigma.e, n = n, loc = loc, nu = nu.vec[i], m.vec = m.vec, n.rep = n.rep, folder_to_save = folder_to_save)
    return(res)
}

err.ss <- err.nn <- err.rat <- nu <- NULL
for(i in 1:length(res)) {
    nu <- c(nu, res[[i]]$nu)
    err.nn <- rbind(err.nn, res[[i]]$err.nn)
    err.rat <- rbind(err.rat, res[[i]]$err.rat)
    
}

res_5000_pred <- list(nu = nu, err.nn = err.nn, err.rat = err.rat)

saveRDS(res_5000_pred, "pred_tables/res_5000_range20_rat_nngp.RDS")




## 



rm(list=ls())
source("aux_functions/aux_functions_cov.R")
source("examples/error.computations.R")
library(rSPDE)
library(foreach)
library(doParallel)
library(doSNOW)

cores=12

cl <- makeCluster(cores[1]-1) 
registerDoSNOW(cl)

sigma = 1
sigma.e <- 0.1


nu.vec <- seq(from = 0.01, to = 2.49, by = 0.01)
nu.vec <- nu.vec[length(nu.vec):1]
m.vec <- 1:6

iterations <- length(nu.vec)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
folder_to_save <- getwd()




n <- 10000
n.obs <- 10000
n.rep <- 100
loc <- seq(0,n/100,length.out=n)
# range <- 0.2
range <- 0.5

res = foreach(i = 1:iterations, .options.snow = opts, .packages=c('Matrix', 'rSPDE', 'pracma')) %dopar% {
    res <- error.computations_nopca_nofourier_noss_n_equal_nobs(range = range, sigma = sigma, sigma.e = sigma.e, n = n, loc = loc, nu = nu.vec[i], m.vec = m.vec, n.rep = n.rep, folder_to_save = folder_to_save)
    return(res)
}

err.ss <- err.nn <- err.rat <- nu <- NULL
for(i in 1:length(res)) {
    nu <- c(nu, res[[i]]$nu)
    err.nn <- rbind(err.nn, res[[i]]$err.nn)
    err.rat <- rbind(err.rat, res[[i]]$err.rat)
    
}

res_10000_pred <- list(nu = nu, err.nn = err.nn, err.rat = err.rat)

saveRDS(res_10000_pred, paste0("pred_tables/res_10000_10000_range",range,"_rat_nngp.RDS"))



## 




rm(list=ls())
source("aux_functions/aux_functions_cov.R")
source("examples/error.computations.R")
library(rSPDE)
library(foreach)
library(doParallel)
library(doSNOW)

cores=12

cl <- makeCluster(cores[1]-1) 
registerDoSNOW(cl)

sigma = 1
sigma.e <- 0.1


nu.vec <- seq(from = 0.01, to = 2.49, by = 0.01)
nu.vec <- nu.vec[length(nu.vec):1]
m.vec <- 1:6

iterations <- length(nu.vec)
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
folder_to_save <- getwd()



n <- 10000
n.obs <- 10000
n.rep <- 100
loc <- seq(0,n/100,length.out=n)
# range <- 0.5
range <- 1

res = foreach(i = 1:iterations, .options.snow = opts, .packages=c('Matrix', 'rSPDE', 'pracma')) %dopar% {
    res <- error.computations_nopca_nofourier_noss_n_equal_nobs(range = range, sigma = sigma, sigma.e = sigma.e, n = n, loc = loc, nu = nu.vec[i], m.vec = m.vec, n.rep = n.rep, folder_to_save = folder_to_save)
    return(res)
}

err.ss <- err.nn <- err.rat <- nu <- NULL
for(i in 1:length(res)) {
    nu <- c(nu, res[[i]]$nu)
    err.nn <- rbind(err.nn, res[[i]]$err.nn)
    err.rat <- rbind(err.rat, res[[i]]$err.rat)
    
}

res_10000_pred <- list(nu = nu, err.nn = err.nn, err.rat = err.rat)

saveRDS(res_10000_pred, paste0("pred_tables/res_10000_10000_range",range,"_rat_nngp.RDS"))



## 

n <- 10000
n.obs <- 10000
n.rep <- 100
loc <- seq(0,n/100,length.out=n)
# range <- 1
range <- 2

res = foreach(i = 1:iterations, .options.snow = opts, .packages=c('Matrix', 'rSPDE', 'pracma')) %dopar% {
    res <- error.computations_nopca_nofourier_noss_n_equal_nobs(range = range, sigma = sigma, sigma.e = sigma.e, n = n, loc = loc, nu = nu.vec[i], m.vec = m.vec, n.rep = n.rep, folder_to_save = folder_to_save)
    return(res)
}

err.ss <- err.nn <- err.rat <- nu <- NULL
for(i in 1:length(res)) {
    nu <- c(nu, res[[i]]$nu)
    err.nn <- rbind(err.nn, res[[i]]$err.nn)
    err.rat <- rbind(err.rat, res[[i]]$err.rat)
    
}

res_10000_pred <- list(nu = nu, err.nn = err.nn, err.rat = err.rat)

saveRDS(res_10000_pred, "pred_tables/res_10000_10000_range1_rat_nngp.RDS")
