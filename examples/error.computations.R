

m_nngp_fun <- function(m, alpha, n, n.obs){
            if(alpha<1) {
                mn <- m - 1
                if(mn < 1){
                    mn <- 1
                }
            } else if (alpha < 2) {
                if(n == 5000){
                    m_vec <- c(1, 2, 13, 21, 27, 31) 
                } else if(n.obs == 5000){
                    m_vec <- c(1, 2, 3, 4, 14, 20)
                } else if(n.obs == 10000){
                    m_vec <- c(1, 2, 7, 18, 24, 29)
                } else{
                    stop("not implemented")
                }
                mn <- m_vec[m]
            } else {
                if(n == 5000){
                    m_vec <- c(15, 30, 37, 45, 51, 54)
                } else if(n.obs == 5000){
                    m_vec <- c(1, 22, 31, 39, 47, 54)
                } else if(n.obs == 10000){
                    m_vec <- c(14, 28, 37, 44, 51, 57)
                } else{
                    stop("not implemented")
                }
                mn <- m_vec[m]
            }
            return(mn)
} 

m_pca_fun <- function(m, alpha, n, n.obs){
            if(alpha<1) {
                if(n == 5000){
                    m_vec <- c(268, 308, 355, 406, 433, 478)
                } else if(n.obs == 5000){
                    m_vec <- c(381, 448, 533, 588, 654, 727)
                } else if (n.obs == 10000){
                    m_vec <- c(271, 311, 367, 412, 452, 493)
                } else{
                    stop("not implemented")
                }
                mn <- m_vec[m]
            } else if (alpha < 2) {
                if(n == 5000){
                    m_vec <-  c(380, 473, 561, 651, 708, 776)
                } else if(n.obs == 5000){
                    m_vec <- c(532, 704, 844, 953, 1065, 1162)
                } else if (n.obs == 10000){
                    m_vec <- c(372, 493, 591, 672, 751, 821)
                } else{
                    stop("not implemented")
                }
                mn <- m_vec[m]
            } else {
                if(n == 5000){
                    m_vec <- c(611, 810, 945, 1082, 1205, 1325)
                } else if(n.obs == 5000){
                    m_vec <- c(904, 1202, 1431, 1622, 1808, 1965)
                } else if (n.obs == 10000){
                    m_vec <- c(640, 848, 1016, 1168, 1299, 1420)
                } else{
                    stop("not implemented")
                }
                mn <- m_vec[m]
            }
            return(mn)
}




m_taper_fun <- function(m, alpha, n, n.obs){
            if(alpha<1) {
                mn <- m - 1
                if(mn < 1){
                    mn <- 1
                }
            } else if (alpha < 2) {
                if(n == 5000){
                    m_vec <- c(1, 2, 3, 62, 124, 166) 
                } else{
                    stop("not implemented")
                }
                mn <- m_vec[m]
            } else {
                if(n == 5000){
                    m_vec <- c(31, 210, 342, 376, 405, 501)
                } else{
                    stop("not implemented")
                }
                mn <- m_vec[m]
            }
            return(mn)
} 


m_fem_fun <- function(m, alpha, n, n.obs){
            if(alpha<1) {
                if(n == 5000){
                m_vec <- c(1, 1, 2, 3, 3, 3)
                } else{
                    stop("not implemented")
                }
                mn <- m_vec[m]                
            } else if (alpha < 2) {
                if(n == 5000){
                    m_vec <- c(2, 4, 6, 7, 8, 9) 
                } else{
                    stop("not implemented")
                }
                mn <- m_vec[m]
            } else {
                if(n == 5000){
                    # m_vec <- c(7, 13, 17, 21, 21, 21)
                    m_vec <- c(5,4,3,2,2,1) # Chose to stabilize the results, not based on cost.
                } else{
                    stop("not implemented")
                }
                mn <- m_vec[m]
            }
            return(mn)
} 

sample_supergauss <- function(kappa, sigma, sigma.e, obs.ind, nu, loc, Sigma){
    loc <- loc - loc[1]
    acf = rSPDE::matern.covariance(h=loc,kappa=kappa,nu=nu,sigma=sigma)
    acf <- as.vector(acf)    
    y <- matrix(rnorm(length(loc)), ncol = length(loc), nrow = 1)
    sim <- SuperGauss::cholZX(Z = t(y), acf = acf)
    sim <- sim + sigma.e*rnorm(length(sim))
    sim <- sim[obs.ind]
    if(any(is.nan(sim))){
            R <- tryCatch(chol(Sigma[obs.ind,obs.ind]), error=function(e){NULL})
            if(!is.null(R)){
                X <- t(R)%*%rnorm(length(obs.ind))
                sim <- as.vector(X + sigma.e*rnorm(length(obs.ind)))
            } else{
                sim <- NULL
            }
    }
    return(sim)
}


## Assumes equally spaced.
error.computations_nopca_nofourier_noss_n_equal_nobs <- function(range, sigma, sigma.e, n, loc, nu, m.vec, n.rep, folder_to_save) {

    set.seed(123)

    m.vec <- 1:6
    err.nn <- err.rat <- err.ss <- matrix(0,nrow=1, ncol = length(m.vec))
    # range <- range * max(loc)

    alpha <- nu + 1/2
    kappa = sqrt(8*nu)/range
    loc2 <- loc - loc[1]
    acf = rSPDE::matern.covariance(h=loc2,kappa=kappa,nu=nu,sigma=sigma)
    acf <- as.vector(acf)
    acf[1] = acf[1]+sigma.e^2
    acf2 =rSPDE::matern.covariance(h=loc2,kappa=kappa,nu=nu,sigma=sigma)
    acf2 =as.vector(acf2)
    # Tz <- SuperGauss::Toeplitz$new(acf = acf)
    # Tz2 <- SuperGauss::Toeplitz$new(acf = acf2)    

    Sigma <- toeplitz(acf2)
    Sigma.hat <- Sigma + sigma.e^2*diag(n)
    chol_tmp <- NULL

    print("range")
    print(range)
    print("nu")
    print(nu)

    for(kk in 1:n.rep) {
        time1 <- Sys.time()
            cat(kk, "True pred\n")
            # obs.ind <- sort(sample(1:n)[1:n.obs])
            obs.ind <- 1:n
            t1 <- Sys.time()
            y <- matrix(rnorm(length(loc)), ncol = length(loc), nrow = 1)
            sim <- SuperGauss::cholZX(Z = t(y), acf = acf2)
            sim <- sim + sigma.e*rnorm(length(sim))
            sim <- sim[obs.ind]
            if(any(is.nan(sim))){
                if(is.null(chol_tmp)){
                    chol_tmp <- tryCatch(chol(Sigma[obs.ind,obs.ind]), error=function(e){NULL})
                }
                if(!is.null(chol_tmp)){
                    X <- t(chol_tmp)%*%rnorm(length(obs.ind))
                    sim <- as.vector(X + sigma.e*rnorm(length(obs.ind)))
                } else{
                    sim <- NULL
                }
            }
            Y <- sim

            # d_tmp = Tz$solve(Y)
            # mu = Tz2$prod(d_tmp)
            if(!is.null(Y)){
                           mu <- Sigma[,obs.ind]%*%solve(Sigma.hat,Y)


            for(j in 1:length(m.vec)) { 
                timem1 <- Sys.time()
                
                m <- m.vec[j]
                
                #########################
                ## Rational prediction
                #########################
                cat(kk, j, "Rational\n")
                if((nu + 0.5)%%1 > 1e-10){

                temp <- tryCatch({
                    Qrat <- tryCatch(rSPDE:::matern.rational.precision(loc = loc, order = m, nu = nu, kappa = kappa, sigma = sigma, 
                                                             cumsum = FALSE, ordering = "location",
                                                             type_rational = "brasil", type_interp =  "spline"), error=function(e){NULL})
                    if(!is.null(Qrat)){
                        # Qrat$Q <- Qrat$Q + Matrix::Diagonal(n = nrow(Qrat$Q), x=1e-8)
                        A_obs <- Qrat$A[obs.ind,]
                        A_mat = Matrix::t(A_obs)
                        Q_xgiveny <-(A_mat%*% (A_obs))/sigma.e^2 + Qrat$Q
                        post_y <- (A_mat%*% Y)/sigma.e^2
                        # Q_xgiveny <- Q_xgiveny + Matrix::Diagonal(n = nrow(Qrat$Q), x=1e-8)
                        R <- tryCatch(Matrix::Cholesky(Q_xgiveny, perm = FALSE) , error=function(e){NULL})
                        if(!is.null(R)){
                            mu_xgiveny <- tryCatch(Matrix::solve(R, post_y, system = "A"), error=function(e){NULL})
                            if(!is.null(mu_xgiveny)){
                                mu.rat <-  Qrat$A %*% mu_xgiveny
                            } else{
                                mu.rat <- NaN
                            }
                        } else{
                            mu.rat <- NaN
                        }
                    } else{
                        mu.rat <- NaN
                    }

                    err.rat[1,j] <- err.rat[1,j] + sqrt((loc[2]-loc[1])*sum((mu-mu.rat)^2))/n.rep
                    1
                }, error=function(e){NULL})
                if(is.null(temp)){
                    err.rat[1,j] <- NaN
                }
                }
                
                print("rational")
                print(err.rat[1,j])

                #########################
                ## nngp prediction
                #########################

                temp2 <- tryCatch({
                mn <- m_nngp_fun(m, alpha, n, n)
                Qnn <- tryCatch(get.nnQ(loc = loc[obs.ind],kappa = kappa,nu = nu,sigma = sigma, n.nbr = mn), error=function(e){NULL})
                if(!is.null(Qnn)){
                    Qhat <- Qnn + Diagonal(n)/sigma.e^2        
                    mu.nn <- tryCatch(solve(Qhat, Y/sigma.e^2), error=function(e){NULL})
                    Bp <- tryCatch(get.nn.pred(loc = loc, kappa = kappa, nu = nu, sigma = sigma, n.nbr = mn, S = obs.ind)$B, error=function(e){NULL})
                    if(!is.null(mu.nn) && !is.null(Bp)){
                        mu.nn <- Bp%*%mu.nn
                    } else{
                        mu.nn <- NaN
                    }
                } else{
                    mu.nn <- NaN
                }

                err.nn[1,j] <- err.nn[1,j] + sqrt((loc[2]-loc[1])*sum((mu-mu.nn)^2))/n.rep
                1
                }, error = function(e){NULL})
                if(is.null(temp2)){
                    err.nn[1,j] <- NaN
                }

                print("nngp")
                print(err.nn[1,j])
                

                # ########################
                # # Statespace prediction
                # #######################
                # cat(kk, j, "Statespace\n")
                # t1 <- Sys.time()
                # ind = 1 + 100*(0:(n-1))
                # h2 = seq(from=0,to=max(loc),length.out=100*(n-1)+1)
                
                # mn <- max(c(1,m - floor(alpha)))
                # coeff <- spec.coeff(kappa = kappa,alpha = nu + 0.5,mn)
                # S1 <- ab2spec(coeff$a,coeff$b,h2, flim = 2)
                # r1 <- S2cov(S1,h2,flim = 2)
                # acf <- r1[ind]
                # acf <- acf * sigma^2
                # cov_mat <- toeplitz(acf, symmetric=TRUE)
                # acf2 <- acf
                # acf2[1] <- acf2[1] + sigma.e^2
                # cov_mat_nugget <-  toeplitz(as.vector(acf2), symmetric=TRUE)
                # cov_mat_nugget <- cov_mat_nugget[obs.ind,obs.ind]
                # d <- solve(cov_mat_nugget, Y)
                # cov_mat <- cov_mat[, obs.ind]
                # mu.ss <- cov_mat%*%d
                # t2 <- Sys.time()
                # print("Statespace time")
                # print(t2 - t1)
                
                # err.ss[1,j] <- err.ss[1,j] + sqrt((loc[2]-loc[1])*sum((mu-mu.ss)^2))/n.rep   
            }
            } else{
                for(j in 1:length(m.vec)) { 
                    err.nn[1,j] <- NaN
                    err.rat[1,j] <- NaN
                }
                break
            }

                    time2 <- Sys.time()
                    print("Time replicate")
                    print(time2-time1)
 
    }
    err.nn <- as.data.frame(err.nn)
    err.rat <- as.data.frame(err.rat)
    err.nn[["nu"]] <- nu
    err.rat[["nu"]] <- nu
    colnames(err.nn) <- c(as.character(m.vec), "nu")
    colnames(err.rat) <- c(as.character(m.vec), "nu")
    res <- list(err.nn = err.nn, 
                err.rat = err.rat,  
                nu = nu)
    dir.create(file.path(folder_to_save, "pred_tables"), showWarnings = FALSE)
    dir.create(file.path(paste0(folder_to_save, "/pred_tables/"), as.character(n)), showWarnings = FALSE)
    dir.create(file.path(paste0(folder_to_save, "/pred_tables/", as.character(n)),paste0("range_",as.character(range))), showWarnings = FALSE)
    saveRDS(res, paste0(folder_to_save,"/pred_tables/",as.character(n),"/range_",as.character(range),"/res_",as.character(nu),"_",as.character(n),"_range_",as.character(range),"_nngp_rat.RDS"))
    return(res)    
}





error.computations <- function(range, sigma, sigma.e, n, n.obs, samples.fourier, loc, nu.vec, m.vec, Dists) {

    set.seed(123)
    m.vec <- 1:6
    err.nn <- err.rat <- err.pca <- err.fourier <- err.ss <- matrix(0,nrow= length(nu.vec), ncol = length(m.vec))

    
    for(i in 1:length(nu.vec)) {
        cat(i/length(nu.vec),"\n")
        nu <- nu.vec[i]    
        alpha <- nu + 1/2
        kappa = sqrt(8*nu)/range
        Sigma <- rSPDE::matern.covariance(h=Dists,kappa=kappa,nu=nu,sigma=sigma)
        cat("Eigen expansion\n")
        eigen_cov <- eigen(Sigma)   
        for(kk in 1:n.rep) {
            cat(i/length(nu.vec), kk, "True pred\n")
            obs.ind <- sort(sample(1:n)[1:n.obs])
            # R <- chol(Sigma[obs.ind,obs.ind])
            # X <- t(R)%*%rnorm(n.obs)
            # Y <- X + sigma.e*rnorm(n.obs)
            Y <- sample_supergauss(kappa = kappa, sigma = sigma, sigma.e = sigma.e, obs.ind = obs.ind, nu = nu, loc = loc)
            Sigma.hat <- Sigma[obs.ind,obs.ind] + sigma.e^2*diag(n.obs)
            mu <- Sigma[,obs.ind]%*%solve(Sigma.hat,Y)
            
            time <- Sys.time()
                        mu <- Sigma[,obs.ind]%*%solve(Sigma.hat,Y)
            time2 <- Sys.time()

            for(j in 1:length(m.vec)) { 
                
                m <- m.vec[j]
                
                #########################
                ## Rational prediction
                #########################
                cat(i/length(nu.vec), kk, j, "Rational\n")
                if((nu + 0.5)%%1 > 1e-10){
                    Qrat <-rSPDE:::matern.rational.precision(loc = loc, order = m, nu = nu, kappa = kappa, sigma = sigma, 
                                                         cumsum = FALSE, ordering = "location",
                                                         type_rational = "brasil", type_interp =  "spline")    
                
                    A_obs <- Qrat$A[obs.ind,]
                    A_mat = Matrix::t(A_obs)
                    Q_xgiveny <-(A_mat%*% (A_obs))/sigma.e^2 + Qrat$Q
                    post_y <- (A_mat%*% Y)/sigma.e^2
                    R <- Matrix::Cholesky(Q_xgiveny, perm = FALSE)         
                    mu_xgiveny <- Matrix::solve(R, post_y, system = "A")
                    mu.rat <-  Qrat$A %*% mu_xgiveny
                
                    err.rat[i,j] <- err.rat[i,j] + sqrt((loc[2]-loc[1])*sum((mu-mu.rat)^2))/n.rep
                } 
                
                #########################
                ## nngp prediction
                #########################
                cat(i/length(nu.vec), kk, j, "nngp\n")
                # if(alpha<1) {
                #     mn <- m - 1
                # } else if (alpha < 2) {
                #     if(m==2) {
                #         mn = 4
                #     } else if(m == 3) {
                #         mn = 14
                #     } else if(m ==4){
                #         mn = 22
                #     } else if(m ==5){
                #         mn = 26
                #     } else if(m==6){
                #         mn = 30
                #     }
                # } else {
                #     if(m==2) {
                #         mn = 30
                #     } else if(m == 3) {
                #         mn = 38
                #     } else if(m ==4){
                #         mn = 44
                #     } else if(m ==5){
                #         mn = 50
                #     } else if(m==6){
                #         mn = 54
                #     }
                # }

                mn <- m_nngp_fun(m, alpha, n, n.obs)
                
                Qnn <- get.nnQ(loc = loc[obs.ind],kappa = kappa,nu = nu,sigma = sigma, n.nbr = mn)
                Qhat <- Qnn + Diagonal(n.obs)/sigma.e^2        
                mu.nn <- solve(Qhat, Y/sigma.e^2)
                Bp <- get.nn.pred(loc = loc, kappa = kappa, nu = nu, sigma = sigma, n.nbr = mn, S = obs.ind)$B
                mu.nn <- Bp%*%mu.nn
                
                err.nn[i,j] <- err.nn[i,j] + sqrt((loc[2]-loc[1])*sum((mu-mu.nn)^2))/n.rep
                
                ########################
                # PCA prediction
                #######################
                cat(i/length(nu.vec), kk, j, "PCA\n")
                # if(alpha<1) {
                #     if(m==2) {
                #         mn = 308
                #     } else if(m == 3) {
                #         mn = 355
                #     } else if(m ==4){
                #         mn = 406
                #     } else if(m ==5){
                #         mn = 433
                #     } else if(m==6){
                #         mn = 478
                #     }
                # } else if (alpha < 2) {
                #     if(m==2) {
                #         mn = 473
                #     } else if(m == 3) {
                #         mn = 561
                #     } else if(m ==4){
                #         mn = 651
                #     } else if(m ==5){
                #         mn = 708
                #     } else if(m==6){
                #         mn = 776
                #     }
                # } else {
                #     if(m==2) {
                #         mn = 810
                #     } else if(m == 3) {
                #         mn = 945
                #     } else if(m ==4){
                #         mn = 1082
                #     } else if(m ==5){
                #         mn = 1205
                #     } else if(m==6){
                #         mn = 1325
                #     }
                # }

                mn <- m_pca_fun(m, alpha, n, n.obs)
                
                K <- eigen_cov$vec[,1:mn]    
                D <- Diagonal(mn,eigen_cov$val[1:mn]) 
                
                Bo <- K[obs.ind,] 
                Q.hat <- solve(D) + t(Bo)%*%Bo/sigma.e^2 
                mu.pca <- K%*%solve(Q.hat, t(Bo)%*%Y/sigma.e^2) 
                
                err.pca[i,j] <- err.pca[i,j] + sqrt((loc[2]-loc[1])*sum((mu-mu.pca)^2))/n.rep
                
                
                #######################
                # Fourier prediction
                #####################
                cat(i/length(nu.vec), kk, j, "Fourier\n")
                # err.tmp = foreach(i = 1:samples.fourier,.combine = '+') %dopar% {
                #     K <-  ff.comp(m = mn, kappa = kappa, alpha = nu + 0.5, loc = loc) * sigma^2
                #     D <- Matrix::Diagonal(x = 1, n = ncol(K))
                #     Bo <- K[obs.ind,] 
                #     Q.hat <- solve(D) + t(Bo)%*%Bo/sigma.e^2 
                #     mu.fourier <- K%*%solve(Q.hat, t(Bo)%*%Y/sigma.e^2)
                #     return(sqrt((loc[2]-loc[1])*sum((mu-mu.fourier)^2))/(n.rep*samples.fourier))
                # }
                # 
                err.tmp <- 0
                for(jj in 1:samples.fourier){
                    K <-  ff.comp(m = mn, kappa = kappa, alpha = nu + 0.5, loc = loc) * sigma^2
                    D <- Matrix::Diagonal(x = 1, n = ncol(K))
                    Bo <- K[obs.ind,]
                    Q.hat <- solve(D) + t(Bo)%*%Bo/sigma.e^2
                    mu.fourier <- K%*%solve(Q.hat, t(Bo)%*%Y/sigma.e^2)
                    err.tmp <-  err.tmp + sqrt((loc[2]-loc[1])*sum((mu-mu.fourier)^2))/(n.rep*samples.fourier)
                }
                err.fourier[i,j] <- err.tmp
                ########################
                # Statespace prediction
                #######################
                cat(i/length(nu.vec), kk, j, "Statespace\n")
                ind = 1 + 100*(0:(n-1))
                h2 = seq(from=0,to=max(loc),length.out=100*(n-1)+1)
                
                mn <- max(c(1,m - floor(alpha)))
                coeff <- spec.coeff(kappa = kappa,alpha = nu + 0.5,mn)
                S1 <- ab2spec(coeff$a,coeff$b,h2, flim = 2)
                r1 <- S2cov(S1,h2,flim = 2)
                acf <- r1[ind]
                acf <- acf * sigma^2
                cov_mat <- toeplitz(acf, symmetric=TRUE)
                acf2 <- acf
                acf2[1] <- acf2[1] + sigma.e^2
                cov_mat_nugget <-  toeplitz(as.vector(acf2), symmetric=TRUE)
                cov_mat_nugget <- cov_mat_nugget[obs.ind,obs.ind]
                d <- solve(cov_mat_nugget, Y)
                cov_mat <- cov_mat[, obs.ind]
                mu.ss <- cov_mat%*%d
                
                err.ss[i,j] <- err.ss[i,j] + sqrt((loc[2]-loc[1])*sum((mu-mu.ss)^2))/n.rep   
            }
        }
    }
    return(list(err.nn = err.nn, 
                err.rat = err.rat, 
                err.pca = err.pca, 
                err.fourier = err.fourier, 
                err.ss = err.ss,
                nu = nu.vec))    
}



error.computations_norat_nonngp <- function(range, sigma, sigma.e, n, n.obs, samples.fourier, loc, nu, m.vec, n.rep, folder_to_save) {

    set.seed(123)
    m.vec <- 1:6
   err.nn <- err.rat <- err.pca <- err.fourier <- err.ss <- matrix(0,nrow= 1, ncol = length(m.vec))
    
    alpha <- nu + 1/2
    kappa = sqrt(8*nu)/range
    loc2 <- loc - loc[1]
    acf = rSPDE::matern.covariance(h=loc2,kappa=kappa,nu=nu,sigma=sigma)
    acf <- as.vector(acf)
    acf[1] = acf[1]+sigma.e^2
    acf2 =rSPDE::matern.covariance(h=loc2,kappa=kappa,nu=nu,sigma=sigma)
    acf2 =as.vector(acf2)
    Sigma <- toeplitz(acf2)
    Sigma.hat <- Sigma + sigma.e^2*diag(n)

    eigen_cov <- eigen(Sigma)   
    for(kk in 1:n.rep) {
        cat(kk, "True pred\n")
        obs.ind <- sort(sample(1:n)[1:n.obs])

        Y <- sample_supergauss_v2(kappa = kappa, sigma = sigma, sigma.e = sigma.e, obs.ind = obs.ind, nu = nu, acf = acf, Sigma = Sigma)
        Sigma.hat <- Sigma[obs.ind,obs.ind] + sigma.e^2*diag(n.obs)
        mu <- Sigma[,obs.ind]%*%solve(Sigma.hat,Y)
        
        for(j in 1:length(m.vec)) { 
            
            m <- m.vec[j]
            
            mn <- m_pca_fun(m, alpha, n, n.obs)
            
            K <- eigen_cov$vec[,1:mn]    
            D <- Diagonal(mn,eigen_cov$val[1:mn]) 
            
            Bo <- K[obs.ind,] 
            Q.hat <- solve(D) + t(Bo)%*%Bo/sigma.e^2 
            mu.pca <- K%*%solve(Q.hat, t(Bo)%*%Y/sigma.e^2) 
            
            err.pca[i,j] <- err.pca[i,j] + sqrt((loc[2]-loc[1])*sum((mu-mu.pca)^2))/n.rep

                print("PCA")
                print(err.pca[1,j])            
            
            #######################
            # Fourier prediction
            #####################
            cat(kk, j, "Fourier\n")

            err.tmp <- 0
            for(jj in 1:samples.fourier){
                K <-  ff.comp(m = mn, kappa = kappa, alpha = nu + 0.5, loc = loc) * sigma^2
                D <- Matrix::Diagonal(x = 1, n = ncol(K))
                Bo <- K[obs.ind,]
                Q.hat <- solve(D) + t(Bo)%*%Bo/sigma.e^2
                mu.fourier <- K%*%solve(Q.hat, t(Bo)%*%Y/sigma.e^2)
                err.tmp <-  err.tmp + sqrt((loc[2]-loc[1])*sum((mu-mu.fourier)^2))/(n.rep*samples.fourier)
            }
            err.fourier[i,j] <- err.tmp

                print("fourier")
                print(err.fourier[1,j])

            ########################
            # Statespace prediction
            #######################
            cat(kk, j, "Statespace\n")
            ind = 1 + 100*(0:(n-1))
            h2 = seq(from=0,to=max(loc),length.out=100*(n-1)+1)
            
            mn <- max(c(1,m - floor(alpha)))
            coeff <- spec.coeff(kappa = kappa,alpha = nu + 0.5,mn)
            S1 <- ab2spec(coeff$a,coeff$b,h2, flim = 2)
            r1 <- S2cov(S1,h2,flim = 2)
            acf <- r1[ind]
            acf <- acf * sigma^2
            cov_mat <- toeplitz(acf, symmetric=TRUE)
            acf2 <- acf
            acf2[1] <- acf2[1] + sigma.e^2
            cov_mat_nugget <-  toeplitz(as.vector(acf2), symmetric=TRUE)
            cov_mat_nugget <- cov_mat_nugget[obs.ind,obs.ind]
            d <- solve(cov_mat_nugget, Y)
            cov_mat <- cov_mat[, obs.ind]
            mu.ss <- cov_mat%*%d
            
            err.ss[i,j] <- err.ss[i,j] + sqrt((loc[2]-loc[1])*sum((mu-mu.ss)^2))/n.rep   

                print("ss")
                print(err.ss[1,j])

        }
        }

    dir.create(file.path(folder_to_save, "pred_tables"), showWarnings = FALSE)
    dir.create(file.path(paste0(folder_to_save, "/pred_tables/"), as.character(n)), showWarnings = FALSE)
    dir.create(file.path(paste0(folder_to_save, "/pred_tables/", as.character(n)),paste0("range_",as.character(range))), showWarnings = FALSE)
    saveRDS(res, paste0(folder_to_save,"/pred_tables/",as.character(n),"/range_",as.character(range),"/res_",as.character(nu),"_",as.character(n),"_range_",as.character(range),"_ss_pca_fourier.RDS"))

    return(list(err.nn = err.nn, 
                err.rat = err.rat, 
                err.pca = err.pca, 
                err.fourier = err.fourier, 
                err.ss = err.ss,
                nu = nu))    
}


error.computations_norat_nonngp <- function(range, sigma, sigma.e, n, n.obs, samples.fourier, loc, nu, m.vec, n.rep, folder_to_save) {

    set.seed(123)
    m.vec <- 1:6
   err.nn <- err.rat <- err.pca <- err.fourier <- err.ss <- matrix(0,nrow= 1, ncol = length(m.vec))
    
    alpha <- nu + 1/2
    kappa = sqrt(8*nu)/range
    loc2 <- loc - loc[1]
    acf = rSPDE::matern.covariance(h=loc2,kappa=kappa,nu=nu,sigma=sigma)
    acf <- as.vector(acf)
    acf[1] = acf[1]+sigma.e^2
    acf2 =rSPDE::matern.covariance(h=loc2,kappa=kappa,nu=nu,sigma=sigma)
    acf2 =as.vector(acf2)
    Sigma <- toeplitz(acf2)
    Sigma.hat <- Sigma + sigma.e^2*diag(n)

    eigen_cov <- eigen(Sigma)   
    for(kk in 1:n.rep) {
        cat(kk, "True pred\n")
        obs.ind <- sort(sample(1:n)[1:n.obs])

        Y <- sample_supergauss_v2(kappa = kappa, sigma = sigma, sigma.e = sigma.e, obs.ind = obs.ind, nu = nu, acf = acf, Sigma = Sigma)
        Sigma.hat <- Sigma[obs.ind,obs.ind] + sigma.e^2*diag(n.obs)
        mu <- Sigma[,obs.ind]%*%solve(Sigma.hat,Y)
        
        for(j in 1:length(m.vec)) { 
            
            m <- m.vec[j]
            
            mn <- m_pca_fun(m, alpha, n, n.obs)
            
            K <- eigen_cov$vec[,1:mn]    
            D <- Diagonal(mn,eigen_cov$val[1:mn]) 
            
            Bo <- K[obs.ind,] 
            Q.hat <- solve(D) + t(Bo)%*%Bo/sigma.e^2 
            mu.pca <- K%*%solve(Q.hat, t(Bo)%*%Y/sigma.e^2) 
            
            err.pca[i,j] <- err.pca[i,j] + sqrt((loc[2]-loc[1])*sum((mu-mu.pca)^2))/n.rep

                print("PCA")
                print(err.pca[1,j])            
            
            #######################
            # Fourier prediction
            #####################
            cat(kk, j, "Fourier\n")

            err.tmp <- 0
            for(jj in 1:samples.fourier){
                K <-  ff.comp(m = mn, kappa = kappa, alpha = nu + 0.5, loc = loc) * sigma^2
                D <- Matrix::Diagonal(x = 1, n = ncol(K))
                Bo <- K[obs.ind,]
                Q.hat <- solve(D) + t(Bo)%*%Bo/sigma.e^2
                mu.fourier <- K%*%solve(Q.hat, t(Bo)%*%Y/sigma.e^2)
                err.tmp <-  err.tmp + sqrt((loc[2]-loc[1])*sum((mu-mu.fourier)^2))/(n.rep*samples.fourier)
            }
            err.fourier[i,j] <- err.tmp

                print("fourier")
                print(err.fourier[1,j])

            ########################
            # Statespace prediction
            #######################
            cat(kk, j, "Statespace\n")
            ind = 1 + 100*(0:(n-1))
            h2 = seq(from=0,to=max(loc),length.out=100*(n-1)+1)
            
            mn <- max(c(1,m - floor(alpha)))
            coeff <- spec.coeff(kappa = kappa,alpha = nu + 0.5,mn)
            S1 <- ab2spec(coeff$a,coeff$b,h2, flim = 2)
            r1 <- S2cov(S1,h2,flim = 2)
            acf <- r1[ind]
            acf <- acf * sigma^2
            cov_mat <- toeplitz(acf, symmetric=TRUE)
            acf2 <- acf
            acf2[1] <- acf2[1] + sigma.e^2
            cov_mat_nugget <-  toeplitz(as.vector(acf2), symmetric=TRUE)
            cov_mat_nugget <- cov_mat_nugget[obs.ind,obs.ind]
            d <- solve(cov_mat_nugget, Y)
            cov_mat <- cov_mat[, obs.ind]
            mu.ss <- cov_mat%*%d
            
            err.ss[i,j] <- err.ss[i,j] + sqrt((loc[2]-loc[1])*sum((mu-mu.ss)^2))/n.rep   

                print("ss")
                print(err.ss[1,j])

        }
        }

    dir.create(file.path(folder_to_save, "pred_tables"), showWarnings = FALSE)
    dir.create(file.path(paste0(folder_to_save, "/pred_tables/"), as.character(n)), showWarnings = FALSE)
    dir.create(file.path(paste0(folder_to_save, "/pred_tables/", as.character(n)),paste0("range_",as.character(range))), showWarnings = FALSE)
    saveRDS(res, paste0(folder_to_save,"/pred_tables/",as.character(n),"/range_",as.character(range),"/res_",as.character(nu),"_",as.character(n),"_range_",as.character(range),"_ss_pca_fourier.RDS"))

    return(list(err.nn = err.nn, 
                err.rat = err.rat, 
                err.pca = err.pca, 
                err.fourier = err.fourier, 
                err.ss = err.ss,
                nu = nu))    
}

error.computations_general <- function(method, range, sigma, sigma.e, n, n.obs, samples.fourier, loc, nu, m.vec, n.rep, folder_to_save) {

    sim_data_name <- "sim_data_result"
    true_mean_name <- "true_mean_result"
    true_sigma_name <- "true_sigma_result"
    nu_vec_python <- "nu_vec"
    obs_ind_python <- "obs_ind_result"

    file_path <- sprintf("python_codes/results/combined_results_n%s_nobs%s_range%s_sigmae%.2f.h5", 
                 n, n.obs, range, sigma.e)

    full_sim_data <- rhdf5::h5read(file_path, sim_data_name)
    full_true_pred <- rhdf5::h5read(file_path, true_mean_name)
    full_true_sigma <- rhdf5::h5read(file_path, true_sigma_name)

    nu_vec_python <- rhdf5::h5read(file_path, nu_vec_python)

    ind_nu <- 249-nu*100+1

    full_sim_data <- full_sim_data[,,ind_nu]
    full_true_pred <- full_true_pred[,,ind_nu]
    full_true_sigma <- full_true_sigma[,,ind_nu]
    
    if (n == 10000 && n.obs == 5000) {
        obs.ind_full <- rhdf5::h5read(file_path, obs_ind_python)
        obs.ind_full <- obs.ind_full[,,ind_nu]
    } else {
        obs.ind <- 1:n
    }

    if (n == 5000) {
        full_sim_data <- full_sim_data[obs.ind,]
        full_true_pred <- full_true_pred[obs.ind,]
        full_true_sigma <- full_true_sigma[obs.ind,]
    }

    set.seed(123)
    m.vec <- 1:6
    err_mu <- matrix(0, nrow = 1, ncol = length(m.vec))
    err_sigma <- matrix(0, nrow = 1, ncol = length(m.vec))

    alpha <- nu + 1 / 2
    kappa <- sqrt(8 * nu) / range

    for (kk in 1:n.rep) {
        cat(kk, "Processing replicate\n")

        if (n == 10000 && n.obs == 5000) {
            obs.ind <- obs.ind_full[, kk] + 1
        }

        Y <- full_sim_data[obs.ind, kk]
        mu_true <- full_true_pred[, kk]
        sigma_true <- full_true_sigma[, kk]

        for (j in 1:length(m.vec)) {
            cat("Processing m =", m.vec[j], "\n")
            m <- m.vec[j]

            if (method == "rational") {
                Qrat <- rSPDE:::matern.rational.precision(loc, order = m, nu, kappa, sigma, 
                                                          cumsum = FALSE, ordering = "location",
                                                          type_rational = "brasil", type_interp = "spline")
                A_obs <- Qrat$A[obs.ind,]
                Q.hat <- t(A_obs) %*% A_obs / sigma.e^2 + Qrat$Q
                mu_est <- Qrat$A %*% Matrix::solve(Matrix::Cholesky(Q.hat, perm = FALSE), 
                                                   t(A_obs) %*% Y / sigma.e^2, system = "A")
                Sigma_post <- Qrat$A %*% Matrix::solve(Q.hat, t(Qrat$A))

            } else if (method == "nngp") {
                mn <- m_nngp_fun(m, alpha, n, n.obs)
                Qnn <- get.nnQ(loc = loc[obs.ind], kappa, nu, sigma, mn)
                Q.hat <- Qnn + Diagonal(n.obs) / sigma.e^2
                mu_obs <- solve(Q.hat, Y / sigma.e^2)
                Bp <- get.nn.pred(loc, kappa, nu, sigma, mn, obs.ind)$B
                mu_est <- Bp %*% mu_obs
                Sigma_post <- Bp %*% solve(Q.hat, t(Bp))

            } else if (method == "taper") {
                mn <- m_taper_fun(m, alpha, n, n.obs)
                cov_mat <- taper_matern_efficient(m = mn, loc = loc, nu = nu, kappa = kappa, sigma = sigma)
                cov_mat_nugget <- cov_mat[obs.ind, obs.ind] + Diagonal(x = sigma.e^2, n = length(obs.ind))
                mu_obs <- solve(cov_mat_nugget, Y)
                mu_est <- cov_mat %*% mu_obs
                Sigma_post <- cov_mat - cov_mat[, obs.ind] %*% solve(cov_mat_nugget, t(cov_mat[, obs.ind]))

            } else if (method == "fem") {
                mn <- m_fem_fun(m, alpha, n, n.obs)
                pred <- fem_pred(Y, loc, obs.ind, nu, kappa, sigma, sigma.e, m, mn)
                mu_est <- pred$mean
                Sigma_post <- Matrix::Diagonal(x = pred$variance)

            } else if (method == "pca") {
                mn <- m_pca_fun(m, alpha, n, n.obs)
                eigen_cov <- eigen_matern_covariance(loc, kappa, sigma, nu)
                K <- eigen_cov$vectors[, 1:mn]
                D <- diag(1 / eigen_cov$values[1:mn])
                Bo <- K[obs.ind,]
                Q.hat <- solve(D + t(Bo) %*% Bo / sigma.e^2)
                mu_est <- K %*% (Q.hat %*% (t(Bo) %*% Y / sigma.e^2))
                Sigma_post <- K %*% solve(Q.hat, t(K))

            } else if (method == "fourier") {
                mn <- m_pca_fun(m, alpha, n, n.obs)
                K <- ff.comp(m = mn, kappa = kappa, alpha = nu + 0.5, loc = loc) * sigma^2
                D <- diag(1 / rep(1, ncol(K)))
                Bo <- K[obs.ind,]
                Q.hat <- solve(D + t(Bo) %*% Bo / sigma.e^2)
                mu_est <- K %*% (Q.hat %*% (t(Bo) %*% Y / sigma.e^2))
                Sigma_post <- K %*% solve(Q.hat, t(K))

            } else if (method == "statespace") {
                ind <- 1 + 100 * (0:(n - 1))
                h2 <- seq(0, max(loc), length.out = 100 * (n - 1) + 1)
                coeff <- spec.coeff(kappa, alpha = nu + 0.5, mn)
                S1 <- ab2spec(coeff$a, coeff$b, h2, flim = 2)
                r1 <- S2cov(S1, h2, flim = 2)
                acf <- r1[ind] * sigma^2
                cov_mat <- toeplitz(acf, symmetric = TRUE)
                cov_mat_nugget <- toeplitz(acf + c(sigma.e^2, rep(0, length(acf) - 1)), symmetric = TRUE)
                d <- solve(cov_mat_nugget[obs.ind, obs.ind], Y)
                mu_est <- cov_mat[, obs.ind] %*% d
                Sigma_post <- cov_mat - cov_mat[, obs.ind] %*% solve(cov_mat_nugget[obs.ind, obs.ind], t(cov_mat[, obs.ind]))
            }

            # Compute errors
            err_mu[1, j] <- err_mu[1, j] + sqrt(mean((mu_true - mu_est)^2)) / n.rep
            err_sigma[1, j] <- err_sigma[1, j] + sqrt(mean((sigma_true - sqrt(diag(Sigma_post)))^2)) / n.rep
            print(method)
            print("mean errors")
            print(err_mu[1,j])
            print("sd errors")
            print(err_sigma[1,j])
        }
    }

    # Save results
    results <- list(
        mu_errors = as.data.frame(err_mu),
        sigma_errors = as.data.frame(err_sigma),
        nu = nu
    )
    colnames(results$mu_errors) <- as.character(m.vec)
    colnames(results$sigma_errors) <- as.character(m.vec)

    save_dir <- file.path(
        folder_to_save, "pred_tables",
        paste0(n, "_", n.obs),
        paste0("range_", range),
        method
    )
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

    saveRDS(
        results,
        file.path(
            save_dir,
            sprintf("res_%.2f_%s_%s_range_%s_sigmae_%.2f_%s.RDS",
                    nu, n, n.obs, range, sigma.e, method)
        )
    )
    return(results)
}



sample_predict_save <- function(n, range, sigma, sigma.e, nu_vec, n_rep){
    sim_data_result <- list()
    true_mu_result <- list()
    loc <- seq(0,n/100,length.out=n)
    for(nu in nu_vec){
        sim_data_result[[as.character(nu)]] <- matrix(ncol=n_rep, nrow=n)
        true_mu_result[[as.character(nu)]] <- matrix(ncol=n_rep, nrow=n)

        time1 <- Sys.time()

        kappa = sqrt(8*nu)/range
        acf = rSPDE::matern.covariance(h=loc,kappa=kappa,nu=nu,sigma=sigma)
        acf <- as.vector(acf)
        Sigma <- toeplitz(acf)
        Sigma.hat <- Sigma + sigma.e^2*diag(n)
        R <- chol(Sigma)
        for(i in 1:n_rep){
            X <- t(R)%*%rnorm(length(loc))
            sim <- as.vector(X + sigma.e*rnorm(length(loc)))
            mu <- Sigma%*%solve(Sigma.hat,sim)
            sim_data_result[[as.character(nu)]][,i] <- sim
            true_mu_result[[as.character(nu)]][,i] <- mu
        }
        tmp_list <- list(sim_data = sim_data_result[[as.character(nu)]], true_mu = true_mu_result[[as.character(nu)]], nu = nu)
        saveRDS(tmp_list, paste0("python_codes/partial_results_R/simulation_results_n",n,"_range",range,"_nu",nu,".RDS"))

        time2 <- Sys.time()
        print(paste("Time for nu = ",nu))
        print(time2-time1)
    }
    result <- list(sim_data_result = sim_data_result, true_mu_result = true_mu_result, nu_vec = nu_vec)
    saveRDS(result, paste0("python_codes/simulation_results_n",n,"_range",range,".RDS"))
}





error.computations_nopca_nofourier_noss <- function(range, sigma, sigma.e, n, n.obs, loc, nu, m.vec, Dists, n.rep) {

    set.seed(123)
    m.vec <- 1:6
    err.nn <- err.rat <- err.ss <- matrix(0,nrow=1, ncol = length(m.vec))
    range <- range * max(loc)

    alpha <- nu + 1/2
    kappa = sqrt(8*nu)/range
    Sigma <- rSPDE::matern.covariance(h=Dists,kappa=kappa,nu=nu,sigma=sigma)

    for(kk in 1:n.rep) {
            cat(kk, "True pred\n")
            obs.ind <- sort(sample(1:n)[1:n.obs])
            Y <- sample_supergauss(kappa = kappa, sigma = sigma, sigma.e = sigma.e, obs.ind = obs.ind, nu = nu, loc = loc)

            Sigma.hat <- Sigma[obs.ind,obs.ind] + sigma.e^2*diag(n.obs)
            mu <- Sigma[,obs.ind]%*%solve(Sigma.hat,Y)            

            for(j in 1:length(m.vec)) { 
                
                m <- m.vec[j]
                
                #########################
                ## Rational prediction
                #########################
                cat(kk, j, "Rational\n")
                if((nu + 0.5)%%1 > 1e-10){
                    Qrat <-rSPDE:::matern.rational.precision(loc = loc, order = m, nu = nu, kappa = kappa, sigma = sigma, 
                                                             cumsum = FALSE, ordering = "location",
                                                             type_rational = "brasil", type_interp =  "spline")    

                    A_obs <- Qrat$A[obs.ind,]
                    A_mat = Matrix::t(A_obs)
                    Q_xgiveny <-(A_mat%*% (A_obs))/sigma.e^2 + Qrat$Q
                    post_y <- (A_mat%*% Y)/sigma.e^2
                    R <- Matrix::Cholesky(Q_xgiveny, perm = FALSE)         
                    mu_xgiveny <- Matrix::solve(R, post_y, system = "A")
                    mu.rat <-  Qrat$A %*% mu_xgiveny

                    err.rat[1,j] <- err.rat[1,j] + sqrt((loc[2]-loc[1])*sum((mu-mu.rat)^2))/n.rep
                }

                #########################
                ## nngp prediction
                #########################

                mn <- m_nngp_fun(m, alpha, n, n.obs)

                Qnn <- get.nnQ(loc = loc[obs.ind],kappa = kappa,nu = nu,sigma = sigma, n.nbr = mn)
                Qhat <- Qnn + Diagonal(n.obs)/sigma.e^2        
                mu.nn <- solve(Qhat, Y/sigma.e^2)
                Bp <- get.nn.pred(loc = loc, kappa = kappa, nu = nu, sigma = sigma, n.nbr = mn, S = obs.ind)$B
                mu.nn <- Bp%*%mu.nn

                err.nn[1,j] <- err.nn[1,j] + sqrt((loc[2]-loc[1])*sum((mu-mu.nn)^2))/n.rep
                

                # ########################
                # # Statespace prediction
                # #######################
                # cat(kk, j, "Statespace\n")
                # t1 <- Sys.time()
                # ind = 1 + 100*(0:(n-1))
                # h2 = seq(from=0,to=max(loc),length.out=100*(n-1)+1)
                
                # mn <- max(c(1,m - floor(alpha)))
                # coeff <- spec.coeff(kappa = kappa,alpha = nu + 0.5,mn)
                # S1 <- ab2spec(coeff$a,coeff$b,h2, flim = 2)
                # r1 <- S2cov(S1,h2,flim = 2)
                # acf <- r1[ind]
                # acf <- acf * sigma^2
                # cov_mat <- toeplitz(acf, symmetric=TRUE)
                # acf2 <- acf
                # acf2[1] <- acf2[1] + sigma.e^2
                # cov_mat_nugget <-  toeplitz(as.vector(acf2), symmetric=TRUE)
                # cov_mat_nugget <- cov_mat_nugget[obs.ind,obs.ind]
                # d <- solve(cov_mat_nugget, Y)
                # cov_mat <- cov_mat[, obs.ind]
                # mu.ss <- cov_mat%*%d
                # t2 <- Sys.time()
                # print("Statespace time")
                # print(t2 - t1)
                
                # err.ss[1,j] <- err.ss[1,j] + sqrt((loc[2]-loc[1])*sum((mu-mu.ss)^2))/n.rep   
            }
    }
    return(list(err.nn = err.nn, 
                err.rat = err.rat,  
                err.ss = err.ss,
                nu = nu.vec))    
}



