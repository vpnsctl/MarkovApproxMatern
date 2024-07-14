# install.packages("mvtnorm")
library(mvtnorm)

# ordered locations
get_cov_mat <- function(loc, m, method, nu, kappa, sigma, samples = NULL, L=NULL){
    if(method == "rat_markov"){
        if(nu < 0.5) {
            Qrat <- matern.rational.ldl(loc = loc, order = m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp =  "spline", equally_spaced = FALSE)    
        } else {
            Qrat <- matern.rational.ldl(loc = loc, order = m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp =  "spline", equally_spaced = FALSE)    
        }
        Q <- t(Qrat$L)%*%Qrat$D%*%Qrat$L
        return(list(A = Qrat$A, Q = Q))
    } else if(method == "nngp"){
            return(get.nnQ(loc = loc,kappa = kappa, nu = nu,sigma = sigma, n.nbr = m))
    } else if(method == "pca"){
        D_loc <- dist2matR(dist(loc))
        cov_mat <- rSPDE::matern.covariance(h=D_loc,kappa=kappa,nu=nu,sigma=sigma)
        eigen_cov <- eigen(cov_mat)       
        K <- eigen_cov$vec[,1:m]    
        D <- diag(eigen_cov$val[1:m])    
        cov_KL <- K%*%D%*%t(K)
        return(cov_KL)
    } else if(method == "fourier"){
        D_loc <- dist2matR(dist(loc_full))
        K <- ff.comp(m = m, kappa = kappa, alpha = nu + 0.5, loc = loc) * sigma^2
        if(samples > 1){
            for(jj in 2:samples){
                K <- K + ff.comp(m = m, kappa = kappa, alpha = nu + 0.5, loc = loc) * sigma^2
            }
        }
        return(K/samples)
    } else if(method == "statespace"){
        flim <- 2
        fact <- 100
        N <- length(loc)
        ind = 1 + fact*(0:(N-1))
        h2 = seq(from=0,to=L,length.out=fact*(N-1)+1)
        coeff <- spec.coeff(kappa = kappa,alpha = nu + 0.5,m)
        S1 <- ab2spec(coeff$a,coeff$b,h2, flim = flim)
        r1 <- S2cov(S1,h2,flim = flim)
        acf <- r1[ind]
        acf <- acf * sigma^2
        cov_mat <- toeplitz(acf, symmetric=TRUE)
        return(cov_mat)
    } else if(method == "true"){
        D_loc <- dist2matR(dist(loc))
        cov_mat <- rSPDE::matern.covariance(h=D_loc,kappa=kappa,nu=nu,sigma=sigma)
        return(cov_mat)
    } else{
        stop("Method not implemented!")
    }
}

posterior_constructor_cov <- function(cov_mat, y, sigma_e){
    cov_mat_nugget <-  cov_mat + Matrix::Diagonal(x=sigma_e^2, n=nrow(cov_mat))
    d <- solve(cov_mat_nugget, y)    
    post_mean <- cov_mat%*%d
    post_cov <- cov_mat - cov_mat %*% solve(cov_mat_nugget, cov_mat)
    return(list(post_mean = post_mean, post_cov = post_cov))
}

posterior_constructor_prec <- function(prec_mat, y, sigma_e, A = NULL){
    if(is.null(A)){        
        I <- Matrix::Diagonal(n = ncol(prec_mat), x = 1)
        Q_xgiveny <- 1/sigma_e^2 * I + prec_mat
        post_y <- y/sigma_e^2
        R <- Matrix::Cholesky(Q_xgiveny, perm=FALSE)
        mu_xgiveny <- solve(R, post_y, system = "A")
    } else{
            A_mat = t(A)
            Q_xgiveny <-(A_mat%*%A)/sigma_e^2 + prec_mat
            post_y <- (A_mat%*% y)/sigma_e^2
            R <- Matrix::Cholesky(Q_xgiveny)            
            mu_xgiveny <- solve(R, post_y, system = "A")
            mu_xgiveny <-  A %*% mu_xgiveny
    }

    return(list(post_mean = mu_xgiveny, post_cov = A%*%solve(Q_xgiveny, A_mat)))
}

compute_prob_ests <- function(loc_full, loc_prob_idx, L, range, sigma, nu_vec, sigma_e, m_rat, 
m_nngp_fun, m_pca_fun, m_fourier_fun, m_statespace_fun, samples_fourier = 10, seed = 1){
    set.seed(seed)
    prob_ests <- list()
    N <- length(loc_full)
    for(nu in nu_vec){
        print(paste0("nu = ", nu))
        alpha <- nu + 0.5  
        kappa <- sqrt(8*nu)/range
        true_cov <- get_cov_mat(loc = loc_full, m = NULL, method = "true", nu = nu, kappa = kappa, sigma = sigma, samples = NULL, L=NULL)
        z <- matrix(rnorm(N), ncol = N, nrow = 1)
        L <- chol(true_cov)
        y <- t(L)%*%t(z)
        post_true <- posterior_constructor_cov(true_cov, y, sigma_e)
        post_mean_true <- post_true$post_mean
        post_cov_true <- post_true$post_cov
        post_mean_true_obs <- post_mean_true[loc_prob_idx]
        post_cov_true_obs <- post_cov_true[loc_prob_idx, loc_prob_idx]
        prob_true <- pmvnorm(lower=-Inf,upper=y[loc_prob_idx],mean=post_mean_true_obs,sigma = as.matrix(post_cov_true_obs))
        prob_true <- c(prob_true)
        prob_ests[[as.character(nu)]] <- list()
        prob_ests[[as.character(nu)]][["true"]] <- prob_true
        m_nngp <- m_nngp_fun(m_rat, alpha)
        m_pca <- m_pca_fun(m_rat, alpha)
        m_fourier <- m_fourier_fun(m_rat, alpha)
        m_statespace <- m_statespace_fun(m_rat, alpha)
        print("true")
        print(prob_true)
        prob_ests[[as.character(nu)]][["rat"]] <- list()
        for(i_m in m_rat){
            prec_rat <- get_cov_mat(loc = loc, m = i_m, method = "rat_markov", nu = nu, kappa = kappa, sigma = sigma, samples = NULL, L=NULL)
            post_rat <- posterior_constructor_prec(prec_rat$Q, y, sigma_e, A = prec_rat$A)
            post_mean_rat <- post_rat$post_mean
            post_cov_rat <- post_rat$post_cov
            post_mean_rat_obs <- post_mean_rat[loc_prob_idx]
            post_cov_rat_obs <- post_cov_rat[loc_prob_idx, loc_prob_idx]
            prob_rat <- pmvnorm(lower=-Inf,upper=y[loc_prob_idx],mean=post_mean_rat_obs,sigma = as.matrix(post_cov_rat_obs))
            prob_rat <- c(prob_rat)
            prob_ests[[as.character(nu)]] <- list()
            prob_ests[[as.character(nu)]][["rat"]][[as.character(i_m)]] <- prob_rat
            print("rat")
            print(prob_rat)

            i_nngp <- m_nngp_fun(i_m, nu + 0.5)

            prec_nngp <- get_cov_mat(loc = loc, m = i_nngp, method = "nngp", nu = nu, kappa = kappa, sigma = sigma, samples = NULL, L=NULL)
            post_nngp <- posterior_constructor_prec(prec_nngp, y, sigma_e)
            post_mean_nngp <- post_nngp$post_mean
            post_cov_nngp <- post_nngp$post_cov
            post_mean_nngp_obs <- post_mean_nngp[loc_prob_idx]
            post_cov_nngp_obs <- post_cov_nngp[loc_prob_idx, loc_prob_idx]
            prob_nngp <- pmvnorm(lower=-Inf,upper=y[loc_prob_idx],mean=post_mean_nngp_obs,sigma = as.matrix(post_cov_nngp_obs))
            prob_nngp <- c(prob_nngp)
            prob_ests[[as.character(nu)]] <- list()
            prob_ests[[as.character(nu)]][["nngp"]][[as.character(i_m)]] <- prob_nngp
            print("nngp")
            print(prob_nngp)

        }

    }
}