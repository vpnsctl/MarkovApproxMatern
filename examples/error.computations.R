
m_nngp_fun <- function(m, alpha, n, n.obs){
            if(alpha<1) {
                mn <- m - 1
                if(mn < 1){
                    mn <- 1
                }
            } else if (alpha < 2) {
                if(n == 5000){
                    m_vec <- c(1, 2, 15, 23, 28, 33)
                } else if(n.obs == 5000){
                    m_vec <- c(1, 2, 3, 12, 18, 23)
                } else if(n.obs == 10000){
                    m_vec <- c(17, 31, 39, 46, 51, 56)
                } else{
                    stop("not implemented")
                }
                mn <- m_vec[m]
            } else {
                if(n == 5000){
                    m_vec <- c(17, 31, 39, 46, 51, 56)
                } else if(n.obs == 5000){
                    m_vec <- c(15, 23, 33, 43, 50, 56)
                } else if(n.obs == 10000){
                    m_vec <- c(15, 31, 40, 48, 55, 60)
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


error.computations <- function(range, sigma, sigma.e, n, n.obs, samples.fourier, loc, nu.vec, m.vec, Dists) {

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
            R <- chol(Sigma[obs.ind,obs.ind])
            X <- t(R)%*%rnorm(n.obs)
            Y <- X + sigma.e*rnorm(n.obs)
            Sigma.hat <- Sigma[obs.ind,obs.ind] + sigma.e^2*diag(n.obs)
            mu <- Sigma[,obs.ind]%*%solve(Sigma.hat,Y)
            
            for(j in 1:length(m.vec)) { 
                
                m <- m.vec[j]
                
                #########################
                ## Rational prediction
                #########################
                cat(i/length(nu.vec), kk, j, "Rational\n")
                Qrat <-rSPDE:::matern.rational.precision(loc = loc, order = m, nu = nu, kappa = kappa, sigma = sigma, 
                                                         cumsum = TRUE, ordering = "location",
                                                         type_rational = "brasil", type_interp =  "spline")    
                
                A_obs <- Qrat$A[obs.ind,]
                A_mat = Matrix::t(A_obs)
                Q_xgiveny <-(A_mat%*% (A_obs))/sigma.e^2 + Qrat$Q
                post_y <- (A_mat%*% Y)/sigma.e^2
                R <- Matrix::Cholesky(Q_xgiveny, perm = FALSE)         
                mu_xgiveny <- Matrix::solve(R, post_y, system = "A")
                mu.rat <-  Qrat$A %*% mu_xgiveny
                
                err.rat[i,j] <- err.rat[i,j] + sqrt((loc[2]-loc[1])*sum((mu-mu.rat)^2))/n.rep
                
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



