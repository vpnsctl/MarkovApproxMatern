rm(list=ls())
source("aux_functions/aux_functions_cov.R")
source("examples/error.computations.R")
library(rSPDE)
library(foreach)
library(doParallel)
library(doSNOW)

cores <- 23

cl <- makeCluster(cores[1]) 
registerDoSNOW(cl)


range <- 2
n <- 10000
n.obs <- 10000

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


n.rep <- 100
loc <- seq(0,n/100,length.out=n)

method <- "rational"

fourier_samples <- 100

res = foreach(i = 1:iterations, .options.snow = opts, .packages=c('Matrix', 'rSPDE', 'pracma', 'SuperGauss','rhdf5')) %dopar% {
    res <- error.computations_general(method = method, range = range, sigma = sigma, sigma.e = sigma.e, n = n, n.obs = n.obs, samples.fourier = fourier_samples, loc = loc, nu = nu.vec[i], m.vec = m.vec, n.rep = n.rep, folder_to_save = folder_to_save)
    return(res)
}

if(method == "rational"){
        err.rat <- nu <- NULL
        for(i in 1:length(res)) {
            nu <- c(nu, res[[i]]$nu)
            err.rat <- rbind(err.rat, res[[i]]$err.rat)
        } 
        res_pred <- list(nu = nu, err.rat = err.rat)        
    } else if(method == "nngp"){
        err.nn <- nu <- NULL
        for(i in 1:length(res)) {
            nu <- c(nu, res[[i]]$nu)
            err.nn <- rbind(err.nn, res[[i]]$err.nn)
        } 
        res_pred <- list(nu = nu, err.nn = err.nn)      
    } else if(method == "pca"){
        err.pca <- nu <- NULL
        for(i in 1:length(res)) {
            nu <- c(nu, res[[i]]$nu)
            err.pca <- rbind(err.pca, res[[i]]$err.pca)
        } 
        res_pred <- list(nu = nu, err.pca = err.pca)    
    } else if(method == "fourier"){
        err.fourier <- nu <- NULL
        for(i in 1:length(res)) {
            nu <- c(nu, res[[i]]$nu)
            err.fourier <- rbind(err.fourier, res[[i]]$err.fourier)
        } 
        res_pred <- list(nu = nu, err.fourier = err.fourier)    
    } else if(method == "statespace"){
        err.ss <- nu <- NULL
        for(i in 1:length(res)) {
            nu <- c(nu, res[[i]]$nu)
            err.ss <- rbind(err.ss, res[[i]]$err.ss)
        } 
        res_pred <- list(nu = nu, err.ss = err.ss)  
    }

saveRDS(res_pred, paste0("markov_approx/pred_tables/res_",as.character(n),"_",as.character(n.obs),"_range",range,"_",as.character(method),".RDS"))