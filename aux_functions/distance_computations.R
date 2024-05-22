# Computes L2 and Linfinity distances for the rational approximation (our method)

compute_distances_rational <- function(N, m.vec, nu.vec, range, sigma){
    L2dist <- list()
    Linfdist <- list()
    for(n_loc in N){
        l2.err <- sup.err <-matrix(0,length(nu.vec),length(m.vec))        
        loc <- seq(0, 1, length.out = n_loc+1)
        for(i in 1:length(nu.vec)) {
            cat(i/length(nu.vec)," ")
            nu <- nu.vec[i]
            alpha <- nu + 0.5  
            kappa <- sqrt(8*nu)/range
            D <- dist2matR(dist(loc))            
            Sigma.t <- matern.covariance(h=D,kappa=kappa,nu=nu,sigma=sigma)
            for(j in 1:length(m.vec)){
                m = m.vec[j]
                if(m == 0){
                    col_tmp <- lindgren_cov(loc, kappa = kappa, beta = (nu + 1/2)/2) * sigma^2
                    Sigma_rat <- toeplitz(x = drop(col_tmp), symmetric = TRUE)
                } else{
                    if(nu < 0.5) {
                        col_tmp <- rSPDE::matern.rational(h = loc, order = m, kappa = kappa, nu = nu, sigma = sigma, type_rational = "brasil", type_interp = "spline")    
                        Sigma_rat <- toeplitz(x = drop(col_tmp), symmetric = TRUE)
                    } else {
                        col_tmp <- rSPDE::matern.rational(h = loc, order = m, kappa = kappa, nu = nu, sigma = sigma, type_rational = "chebfun", type_interp = "spline")
                        Sigma_rat <- toeplitz(x = drop(col_tmp), symmetric = TRUE)
                    }
                }
                l2.err[i,j] <- sqrt(sum((Sigma.t-Sigma_rat)^2))*(loc[2]-loc[1])
                sup.err[i,j] <- max(abs(Sigma.t-Sigma_rat))
            }
        }
        L2dist[[as.character(n_loc)]] <- l2.err
        Linfdist[[as.character(n_loc)]] <- sup.err
    }
    ret <- list(L2 = L2dist, Linf = Linfdist)
    attr(ret, "type") <- "Rational"
    attr(ret, "nu.vec") <- nu.vec
    attr(ret, "m.vec") <- m.vec
    return(ret)
}

# Statespace

compute_distances_statespace <- function(N, m.vec, nu.vec, range, sigma, flim = 2, fact = 100, type = "prediction"){
    L2dist <- list()
    Linfdist <- list()
    for(n_loc in N){
        l2.err <- sup.err <-matrix(0,length(nu.vec),length(m.vec))        
        loc <- seq(0, 1, length.out = n_loc+1)
        ind = 1 + fact*(0:n_loc)
        h2 = seq(from=0,to=1,length.out= fact*n_loc+1)
        D <- dist2matR(dist(loc))
        for(i in 1:length(nu.vec)) {
            cat(i/length(nu.vec)," ")
            nu <- nu.vec[i]
            alpha <- nu + 0.5  
            kappa <- sqrt(8*nu)/range
            Sigma.t <- matern.covariance(h=D,kappa=kappa,nu=nu,sigma=1)
            for(j in 1:length(m.vec)){
                m = m.vec[j]
                m <- get_m(nu = nu, m = m, method = "statespace", type = type)
                coeff <- spec.coeff(kappa,alpha,m)
                S1 <- ab2spec(coeff$a,coeff$b,h2, flim = flim)
                r1 <- S2cov(S1,h2,flim = flim)
                Sigma_nn <- toeplitz(r1[ind] * sigma^2)
                l2.err[i,j] <- sqrt(sum((Sigma.t-Sigma_nn)^2))*(loc[2]-loc[1])
                sup.err[i,j] <- max(abs(Sigma.t-Sigma_nn))
            }
        }
        L2dist[[as.character(n_loc)]] <- l2.err
        Linfdist[[as.character(n_loc)]] <- sup.err
    }
    ret <- list(L2 = L2dist, Linf = Linfdist)
    attr(ret, "type") <- "State-Space"
    attr(ret, "nu.vec") <- nu.vec
    attr(ret, "m.vec") <- m.vec    
    return(ret)
}


# nnGP 

compute_distances_nngp <- function(N, m.vec, nu.vec, range, sigma, type = "prediction"){
    L2dist <- list()
    Linfdist <- list()
    for(n_loc in N){
        l2.err <- sup.err <-matrix(0,length(nu.vec),length(m.vec))        
        loc <- seq(0, 1, length.out = n_loc+1)
        D <- dist2matR(dist(loc))        
        for(i in 1:length(nu.vec)) {
            cat(i/length(nu.vec)," ")
            nu <- nu.vec[i]
            alpha <- nu + 0.5  
            kappa <- sqrt(8*nu)/range
            Sigma.t <- matern.covariance(h=D,kappa=kappa,nu=nu,sigma=sigma)
            for(j in 1:length(m.vec)){
                m = m.vec[j]
                m <- get_m(nu = nu, m = m, method = "nngp", type = type)
                Q.nn <- tryCatch(get.nnQ(loc=loc,kappa=kappa,nu=nu,sigma=sigma, n.nbr=m), error=function(e){NULL})
                if(!is.null(Q.nn)){
                    Sigma_nn <- tryCatch(solve(Q.nn), error=function(e){NULL})
                    if(!is.null(Sigma_nn)){
                        l2.err[i,j] <- sqrt(sum((Sigma.t-Sigma_nn)^2))*(loc[2]-loc[1])
                        sup.err[i,j] <- max(abs(Sigma.t-Sigma_nn))      
                    } else{
                        l2.err[i,j] <- NaN
                        sup.err[i,j] <- NaN
                    }
                } else{
                    l2.err[i,j] <- NaN
                    sup.err[i,j] <- NaN                    
                }
            }
        }
        L2dist[[as.character(n_loc)]] <- l2.err
        Linfdist[[as.character(n_loc)]] <- sup.err
    }
    ret <- list(L2 = L2dist, Linf = Linfdist)
    attr(ret, "type") <- "nnGP"
    attr(ret, "nu.vec") <- nu.vec
    attr(ret, "m.vec") <- m.vec    
    return(ret)
}

# PCA

compute_distances_pca <- function(N, m.vec, nu.vec, range, sigma, type = "prediction"){
    L2dist <- list()
    Linfdist <- list()
    for(n_loc in N){
        l2.err <- sup.err <-matrix(0,length(nu.vec),length(m.vec))        
        loc <- seq(0, 1, length.out = n_loc)
        D <- dist2matR(dist(loc))        
        for(i in 1:length(nu.vec)) {
            cat(i/length(nu.vec)," ")
            nu <- nu.vec[i]
            alpha <- nu + 0.5  
            kappa <- sqrt(8*nu)/range
            Sigma.t <- matern.covariance(h=D,kappa=kappa,nu=nu,sigma=sigma)
            eigen_cov <- eigen(Sigma.t)
            for(j in 1:length(m.vec)){
                m = m.vec[j]
                m <- get_m(nu = nu, m = m, method = "kl", type = type)
                Sigma_KL <- KL_matern(m=m, eigen_cov = eigen_cov)    
                l2.err[i,j] <- sqrt(sum((Sigma.t-Sigma_KL)^2))*(loc[2]-loc[1])
                sup.err[i,j] <- max(abs(Sigma.t-Sigma_KL))               
            }
        }
        L2dist[[as.character(n_loc)]] <- l2.err
        Linfdist[[as.character(n_loc)]] <- sup.err
    }
    ret <- list(L2 = L2dist, Linf = Linfdist)
    attr(ret, "type") <- "PCA"
    attr(ret, "nu.vec") <- nu.vec
    attr(ret, "m.vec") <- m.vec    
    return(ret)
}


# kl

compute_distances_kl <- function(N, m.vec, nu.vec, range, sigma, N_KL=10000, type = "prediction"){
    L2dist <- list()
    Linfdist <- list()
    for(n_loc in N){
        l2.err <- sup.err <-matrix(0,length(nu.vec),length(m.vec))        
        loc <- seq(0, 1, length.out = n_loc)
        large_KL <- seq(min(loc), max(loc), length.out = N_KL)
        kl_loc <- c(loc, large_KL)
        kl_loc <- unique(kl_loc)
        D_loc <- dist2matR(dist(kl_loc))
        for(i in 1:length(nu.vec)) {
            nu <- nu.vec[i]
            cat(i/length(nu.vec)," ")
            alpha <- nu + 0.5  
            kappa <- sqrt(8*nu)/range
            cov_mat <- rSPDE::matern.covariance(h=D_loc,kappa=kappa,nu=nu,sigma=sigma)
            eigen_cov <- eigen(cov_mat) 
            eigen_cov$vec <- eigen_cov$vec[1:N,]            
            Sigma.t <- cov_mat[1:N, 1:N]
            for(j in 1:length(m.vec)){
                m = m.vec[j]
                m <- get_m(nu = nu, m = m, method = "kl", type = type)
                Sigma_KL <- KL_matern(m=m, eigen_cov = eigen_cov)    
                l2.err[i,j] <- sqrt(sum((Sigma.t-Sigma_KL)^2))*(loc[2]-loc[1])
                sup.err[i,j] <- max(abs(Sigma.t-Sigma_KL))         
            }
        }
        L2dist[[as.character(n_loc)]] <- l2.err
        Linfdist[[as.character(n_loc)]] <- sup.err
    }
    ret <- list(L2 = L2dist, Linf = Linfdist)
    attr(ret, "type") <- "KL"
    attr(ret, "nu.vec") <- nu.vec
    attr(ret, "m.vec") <- m.vec    
    return(ret)
}


# Fourier
compute_distances_fourier <- function(N, m.vec, nu.vec, range, sigma, samples, type = "prediction"){
    L2dist <- list()
    Linfdist <- list()
    for(n_loc in N){
        l2.err <- sup.err <-matrix(0,length(nu.vec),length(m.vec))        
        loc <- seq(0, 1, length.out = n_loc+1)
        D <- dist2matR(dist(loc))        
        for(i in 1:length(nu.vec)) {
            cat(i/length(nu.vec)," ")
            nu <- nu.vec[i]
            alpha <- nu + 0.5  
            kappa <- sqrt(8*nu)/range
            Sigma.t <- matern.covariance(h=D,kappa=kappa,nu=nu,sigma=1)
            for(j in 1:length(m.vec)){
                m = m.vec[j]
                m <- get_m(nu = nu, m = m, method = "kl", type = type)                
                for(k in 1:samples){
                    Sigma_fou <- ff.approx(m=m, kappa=kappa, alpha = alpha, loc = loc) * sigma^2      
                    l2.err[i,j] <- l2.err[i,j] + sqrt(sum((Sigma.t-Sigma_fou)^2))*(loc[2]-loc[1])/samples
                    sup.err[i,j] <- sup.err[i,j] + max(abs(Sigma.t-Sigma_fou))/samples
                }
          
            }
        }
        L2dist[[as.character(n_loc)]] <- l2.err
        Linfdist[[as.character(n_loc)]] <- sup.err
    }
    ret <- list(L2 = L2dist, Linf = Linfdist)
    attr(ret, "type") <- "Fourier"
    attr(ret, "nu.vec") <- nu.vec
    attr(ret, "m.vec") <- m.vec    
    return(ret)
}