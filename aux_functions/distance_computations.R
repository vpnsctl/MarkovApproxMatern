# Computes L2 and Linfinity distances for the rational approximation (our method)
# Range is percentage of domain length


compute_distances_rational <- function(N, n_obs, m.vec, nu.vec, range, sigma){
    print("rational")
    N <- N[[1]]
    n_obs <- n_obs[[1]]
    loc <- seq(0, N/100, length.out = N)
    D <- as.matrix(dist(loc))  
    l2.err <- sup.err <-matrix(0,length(nu.vec),length(m.vec))        
    # range <- range * max(loc)
    print("Rational")
    print(paste("N = ",N))
    print(paste("n_obs = ", n_obs))
    print(paste("Domain length = ", max(loc)))
    print(paste("range = ", range))       
    for(i in 1:length(nu.vec)) {
            cat(i/length(nu.vec)," ")
            nu <- nu.vec[i]
            alpha <- nu + 0.5  
            kappa <- sqrt(8*nu)/range          
            Sigma.t <- matern.covariance(h=D,kappa=kappa,nu=nu,sigma=sigma)
            for(j in 1:length(m.vec)){
                m = m.vec[j]
                if(m == 0){
                    col_tmp <- lindgren_cov(loc, kappa = kappa, beta = (nu + 1/2)/2) * sigma^2
                    Sigma_rat <- toeplitz(x = drop(col_tmp), symmetric = TRUE)
                } else{
                    if(nu < 0.5) {
                        col_tmp <- rSPDE:::matern.rational.cov(h = loc, order = m, kappa = kappa, nu = nu, sigma = sigma, type_rational = "brasil", type_interp = "spline")    
                        Sigma_rat <- toeplitz(x = drop(col_tmp), symmetric = TRUE)
                    } else {
                        col_tmp <- rSPDE:::matern.rational.cov(h = loc, order = m, kappa = kappa, nu = nu, sigma = sigma, type_rational = "chebfun", type_interp = "spline")
                        Sigma_rat <- toeplitz(x = drop(col_tmp), symmetric = TRUE)
                    }
                }
                if(n_obs < N){
                    l2.err[i,j] <- sqrt(sum((Sigma.t[1:n_obs,]-Sigma_rat[1:n_obs,])^2))*(loc[2]-loc[1])
                    sup.err[i,j] <- max(abs(Sigma.t[1:n_obs,]-Sigma_rat[1:n_obs,]))
                } else{
                    l2.err[i,j] <- sqrt(sum((Sigma.t-Sigma_rat)^2))*(loc[2]-loc[1])
                    sup.err[i,j] <- max(abs(Sigma.t-Sigma_rat))                    
                }
            }
    }
    ret <- list(L2 = l2.err, Linf = sup.err)
    attr(ret, "type") <- "Rational"
    attr(ret, "nu.vec") <- nu.vec
    attr(ret, "m.vec") <- m.vec
    attr(ret, "N") <- N
    attr(ret, "n_obs") <- n_obs    
    return(ret)
}

# Statespace

compute_distances_statespace <- function(N, n_obs, m.vec, nu.vec, range, sigma, flim = 2, fact = 100, m_statespace_fun){
    N <- N[[1]]
    n_obs <- n_obs[[1]]
    l2.err <- sup.err <-matrix(0,length(nu.vec),length(m.vec))        
    loc <- seq(0, N/100, length.out = N)
    ind = 1 + fact*(0:(N-1))
    h2 = seq(from=0,to=max(loc),length.out= fact*(N-1)+1)
    D <- as.matrix(dist(loc))
    # range <- range * max(loc)    
    print("StateSpace")
    print(paste("N = ",N))
    print(paste("n_obs = ", n_obs))
    print(paste("Domain length = ", max(loc)))
    print(paste("range = ", range))    
    for(i in 1:length(nu.vec)) {
        cat(i/length(nu.vec)," ")
        nu <- nu.vec[i]
        alpha <- nu + 0.5  
        kappa <- sqrt(8*nu)/range
        Sigma.t <- matern.covariance(h=D,kappa=kappa,nu=nu,sigma=1)
        for(j in 1:length(m.vec)){
                m = m.vec[j]
                m <- m_statespace_fun(m, alpha)
                coeff <- spec.coeff(kappa,alpha,m)
                S1 <- ab2spec(coeff$a,coeff$b,h2, flim = flim)
                r1 <- S2cov(S1,h2,flim = flim)
                Sigma_nn <- toeplitz(r1[ind] * sigma^2)
                if(n_obs < N){
                    l2.err[i,j] <- sqrt(sum((Sigma.t[1:n_obs,]-Sigma_nn[1:n_obs,])^2))*(loc[2]-loc[1])
                    sup.err[i,j] <- max(abs(Sigma.t[1:n_obs,]-Sigma_nn[1:n_obs,]))
                } else{
                    l2.err[i,j] <- sqrt(sum((Sigma.t-Sigma_nn)^2))*(loc[2]-loc[1])
                    sup.err[i,j] <- max(abs(Sigma.t-Sigma_nn))                    
                }         
            }
        }
    ret <- list(L2 = l2.err, Linf = sup.err)
    attr(ret, "type") <- "State-Space"
    attr(ret, "nu.vec") <- nu.vec
    attr(ret, "m.vec") <- m.vec    
    attr(ret, "N") <- N
    attr(ret, "n_obs") <- n_obs
    return(ret)
}


# nnGP 

compute_distances_nngp <- function(N, n_obs, m.vec, nu.vec, range, sigma, m_nngp_fun, invert_prec = TRUE){
    N <- N[[1]]
    n_obs <- n_obs[[1]]
    l2.err <- sup.err <-matrix(0,length(nu.vec),length(m.vec))        
    loc <- seq(0, N/100, length.out = N)
    D <- as.matrix(dist(loc))        
    # range <- range * max(loc)    
    print("NNGP")
    print(paste("N = ",N))
    print(paste("n_obs = ", n_obs))
    print(paste("Domain length = ", max(loc)))
    print(paste("range = ", range))
    for(i in 1:length(nu.vec)) {
            cat(i/length(nu.vec)," ")
            nu <- nu.vec[i]
            alpha <- nu + 0.5  
            kappa <- sqrt(8*nu)/range
            Sigma.t <- matern.covariance(h=D,kappa=kappa,nu=nu,sigma=sigma)
            for(j in 1:length(m.vec)){
                m = m.vec[j]
                m <- m_nngp_fun(m, alpha)
                Sigma_nn <- tryCatch(get.nnCov(loc=loc,kappa=kappa,nu=nu,sigma=sigma, n.nbr=m, invert_prec = invert_prec), error=function(e){NULL})
                if(!is.null(Sigma_nn)){
                    if(n_obs < N){
                        l2.err[i,j] <- sqrt(sum((Sigma.t[1:n_obs,]-Sigma_nn[1:n_obs,])^2))*(loc[2]-loc[1])
                        sup.err[i,j] <- max(abs(Sigma.t[1:n_obs,]-Sigma_nn[1:n_obs,]))      
                    } else{
                        l2.err[i,j] <- sqrt(sum((Sigma.t-Sigma_nn)^2))*(loc[2]-loc[1])
                        sup.err[i,j] <- max(abs(Sigma.t-Sigma_nn))                             
                    }
                } else{
                    l2.err[i,j] <- NaN
                    sup.err[i,j] <- NaN
                }          
            }
    }
    ret <- list(L2 = l2.err, Linf = sup.err)
    attr(ret, "type") <- "nnGP"
    attr(ret, "nu.vec") <- nu.vec
    attr(ret, "m.vec") <- m.vec    
    attr(ret, "N") <- N
    attr(ret, "n_obs") <- n_obs    
    return(ret)
}

# PCA

compute_distances_pca <- function(N, n_obs, m.vec, nu.vec, range, sigma, m_pca_fun){
    N <- N[[1]]
    n_obs <- n_obs[[1]]
    l2.err <- sup.err <-matrix(0,length(nu.vec),length(m.vec))        
    loc <- seq(0, N/100, length.out = N)
    D <- as.matrix(dist(loc))     
    # range <- range * max(loc)      
    print("PCA")
    print(paste("N = ",N))
    print(paste("n_obs = ", n_obs))
    print(paste("Domain length = ", max(loc)))
    print(paste("range = ", range))        
    for(i in 1:length(nu.vec)) {
            cat(i/length(nu.vec)," ")
            nu <- nu.vec[i]      
            alpha <- nu + 0.5  
            kappa <- sqrt(8*nu)/range
            Sigma.t <- matern.covariance(h=D,kappa=kappa,nu=nu,sigma=sigma)
            eigen_cov <- eigen(Sigma.t)
            for(j in 1:length(m.vec)){
                m = m.vec[j]
                m <- m_pca_fun(m, alpha)
                Sigma_KL <- KL_matern(m=m, eigen_cov = eigen_cov)    
                if(n_obs < N){
                    l2.err[i,j] <- sqrt(sum((Sigma.t[1:n_obs,]-Sigma_KL[1:n_obs,])^2))*(loc[2]-loc[1])
                    sup.err[i,j] <- max(abs(Sigma.t[1:n_obs,]-Sigma_KL[1:n_obs,]))                           
                } else{
                    l2.err[i,j] <- sqrt(sum((Sigma.t-Sigma_KL)^2))*(loc[2]-loc[1])
                    sup.err[i,j] <- max(abs(Sigma.t-Sigma_KL))               
                }          
            }
    }
    ret <- list(L2 = l2.err, Linf = sup.err)
    attr(ret, "type") <- "PCA"
    attr(ret, "nu.vec") <- nu.vec
    attr(ret, "m.vec") <- m.vec    
    attr(ret, "N") <- N
    attr(ret, "n_obs") <- n_obs    
    return(ret)
}

# Fourier
compute_distances_fourier <- function(N, n_obs, m.vec, nu.vec, range, sigma, samples, m_fourier_fun){
    N <- N[[1]]
    n_obs <- n_obs[[1]]
    l2.err <- sup.err <-matrix(0,length(nu.vec),length(m.vec))        
    loc <- seq(0, N/100, length.out = N)
    D <- as.matrix(dist(loc))        
    # range <- range * max(loc)    
    print("Fourier")
    print(paste("N = ",N))
    print(paste("n_obs = ", n_obs))
    print(paste("Domain length = ", max(loc)))
    print(paste("range = ", range))       
    for(i in 1:length(nu.vec)) {
            cat(i/length(nu.vec)," ")
            nu <- nu.vec[i]         
            alpha <- nu + 0.5  
            kappa <- sqrt(8*nu)/range
            Sigma.t <- matern.covariance(h=D,kappa=kappa,nu=nu,sigma=sigma)
            for(j in 1:length(m.vec)){
                m = m.vec[j]
                m <- m_fourier_fun(m, alpha)
                Sigma_fou <- matrix(0, ncol=ncol(Sigma.t), nrow=nrow(Sigma.t))            
                for(k in 1:samples){
                    Sigma_fou <- Sigma_fou + ff.approx(m=m, kappa=kappa, alpha = alpha, loc = loc) * sigma^2      
                }
                Sigma_fou <- Sigma_fou/samples
                if(n_obs < N){
                    l2.err[i,j] <- sqrt(sum((Sigma.t[1:n_obs,]-Sigma_fou[1:n_obs,])^2))*(loc[2]-loc[1])
                    sup.err[i,j] <-max(abs(Sigma.t[1:n_obs,]-Sigma_fou[1:n_obs,]))                    
                } else{
                    l2.err[i,j] <- sqrt(sum((Sigma.t-Sigma_fou)^2))*(loc[2]-loc[1])
                    sup.err[i,j] <-max(abs(Sigma.t-Sigma_fou))
                }            
            }
        }
    ret <- list(L2 = l2.err, Linf = sup.err)
    attr(ret, "type") <- "Fourier"
    attr(ret, "nu.vec") <- nu.vec
    attr(ret, "m.vec") <- m.vec    
    attr(ret, "N") <- N
    attr(ret, "n_obs") <- n_obs    
    return(ret)
}


compute_distances_taper <- function(N, n_obs, m.vec, nu.vec, range, sigma, m_taper_fun){
    N <- N[[1]]
    n_obs <- n_obs[[1]]
    l2.err <- sup.err <-matrix(0,length(nu.vec),length(m.vec))        
    loc <- seq(0, N/100, length.out = N)
    D <- as.matrix(dist(loc))     
    # range <- range * max(loc)      
    print("Taper")
    print(paste("N = ",N))
    print(paste("n_obs = ", n_obs))
    print(paste("Domain length = ", max(loc)))
    print(paste("range = ", range))        
    for(i in 1:length(nu.vec)) {
            cat(i/length(nu.vec)," ")
            nu <- nu.vec[i]      
            alpha <- nu + 0.5  
            kappa <- sqrt(8*nu)/range
            Sigma.t <- matern.covariance(h=D,kappa=kappa,nu=nu,sigma=sigma)
            eigen_cov <- eigen(Sigma.t)
            for(j in 1:length(m.vec)){
                m = m.vec[j]
                m <- m_taper_fun(m = m, alpha = alpha, n = N, n.obs = n_obs)
                Sigma_taper <- taper_matern_efficient(m=m, loc = loc, nu = nu, kappa = kappa, sigma = sigma)    
                if(n_obs < N){
                    l2.err[i,j] <- sqrt(sum((Sigma.t[1:n_obs,]-Sigma_taper[1:n_obs,])^2))*(loc[2]-loc[1])
                    sup.err[i,j] <- max(abs(Sigma.t[1:n_obs,]-Sigma_taper[1:n_obs,]))                           
                } else{
                    l2.err[i,j] <- sqrt(sum((Sigma.t-Sigma_taper)^2))*(loc[2]-loc[1])
                    sup.err[i,j] <- max(abs(Sigma.t-Sigma_taper))               
                }          
            }
        saveRDS(list(L2 = l2.err, Linf = sup.err), paste0("distance_tables/raw_tables/partials/res_",N,"_",n_obs,"_nu_",nu,"_range_",range,".RDS"))

    }
    ret <- list(L2 = l2.err, Linf = sup.err)
    attr(ret, "type") <- "Taper"
    attr(ret, "nu.vec") <- nu.vec
    attr(ret, "m.vec") <- m.vec    
    attr(ret, "N") <- N
    attr(ret, "n_obs") <- n_obs    
    return(ret)
}


compute_distances_fem <- function(N, n_obs, m.vec, nu.vec, range, sigma, m_fem_fun){
    N <- N[[1]]
    n_obs <- n_obs[[1]]
    l2.err <- sup.err <-matrix(0,length(nu.vec),length(m.vec))        
    loc <- seq(0, N/100, length.out = N)
    D <- as.matrix(dist(loc))     
    # range <- range * max(loc)      
    print("FEM")
    print(paste("N = ",N))
    print(paste("n_obs = ", n_obs))
    print(paste("Domain length = ", max(loc)))
    print(paste("range = ", range))        
    for(i in 1:length(nu.vec)) {
            cat(i/length(nu.vec)," ")
            nu <- nu.vec[i]      
            alpha <- nu + 0.5  
            kappa <- sqrt(8*nu)/range
            Sigma.t <- matern.covariance(h=D,kappa=kappa,nu=nu,sigma=sigma)
            eigen_cov <- eigen(Sigma.t)
            for(j in 1:length(m.vec)){
                m_rat = m.vec[j]
                m <- m_fem_fun(m_rat, alpha, n = N, n.obs = N_obs)
                Sigma_fem <- fem_cov(m=m_rat, mesh_fem = m, loc = loc, nu = nu, kappa = kappa, sigma = sigma)    
                if(n_obs < N){
                    l2.err[i,j] <- sqrt(sum((Sigma.t[1:n_obs,]-Sigma_fem[1:n_obs,])^2))*(loc[2]-loc[1])
                    sup.err[i,j] <- max(abs(Sigma.t[1:n_obs,]-Sigma_fem[1:n_obs,]))                           
                } else{
                    l2.err[i,j] <- sqrt(sum((Sigma.t-Sigma_fem)^2))*(loc[2]-loc[1])
                    sup.err[i,j] <- max(abs(Sigma.t-Sigma_fem))               
                }          
            }
            saveRDS(list(L2 = l2.err, Linf = sup.err), paste0("distance_tables/raw_tables/partials_fem/res_",N,"_",n_obs,"_nu_",nu,"_range_",range,".RDS"))
    }
    ret <- list(L2 = l2.err, Linf = sup.err)
    attr(ret, "type") <- "Taper"
    attr(ret, "nu.vec") <- nu.vec
    attr(ret, "m.vec") <- m.vec    
    attr(ret, "N") <- N
    attr(ret, "n_obs") <- n_obs    
    return(ret)
}