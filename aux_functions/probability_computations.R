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
        Sigma_fou <- matrix(0, ncol=length(loc), nrow=length(loc)) 
        if(samples > 1){
            for(jj in 2:samples){
                Sigma_fou <- Sigma_fou + ff.approx(m = m, kappa = kappa, alpha = nu + 0.5, loc = loc) * sigma^2
            }
        }
        return(Sigma_fou/samples)
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

posterior_constructor_cov <- function(cov_mat, y, sigma_e, obs_ind, pred_ind, type){
    cov_mat_nugget <-  cov_mat + Matrix::Diagonal(x=sigma_e^2, n=nrow(cov_mat))
    cov_mat_nugget <- cov_mat_nugget[obs_ind, obs_ind]
    d <- solve(cov_mat_nugget, y)    
    if(type == "unobserved"){
        post_mean <- cov_mat[pred_ind, obs_ind]%*%d
        post_cov <- cov_mat[pred_ind,pred_ind] - cov_mat[pred_ind,obs_ind] %*% solve(cov_mat_nugget, t(cov_mat[pred_ind,obs_ind]))
    } else{
        post_mean <- cov_mat[, obs_ind]%*%d
        post_cov <- cov_mat - cov_mat[,obs_ind] %*% solve(cov_mat_nugget, t(cov_mat[,obs_ind]))
    }

    return(list(post_mean = post_mean, post_cov = post_cov))
}

posterior_constructor_prec <- function(prec_mat, y, sigma_e, A = NULL, idx_obs, idx_pred, type){
    if(is.null(A)){        
        A <-  Matrix::Diagonal(n = ncol(prec_mat), x = 1)
    }
        if(type == "unobserved"){
            A_pred <- A[idx_pred,]
        } else{
            A_pred <- A
        }

        A_obs <- A[idx_obs,]
        A_mat <- t(A_obs)
        Q_xgiveny <- (A_mat%*%A_obs)/sigma_e^2 + prec_mat
        post_y <- (A_mat%*%y)/sigma_e^2
        R <- Matrix::Cholesky(Q_xgiveny)            
        mu_xgiveny <- solve(R, post_y, system = "A")
        mu_xgiveny <-  A_pred%*%mu_xgiveny
        return(list(post_mean = mu_xgiveny, post_cov = A_pred%*%solve(Q_xgiveny, t(A_pred))))         
}

posterior_constructor_nngp <- function(prec_mat, y, sigma_e, idx_pred, obs.ind, loc, i_m, nu, kappa, sigma, type){ 
    n.obs <- length(obs.ind)
    Qhat <- prec_mat + Diagonal(n.obs)/sigma_e^2        
    mu.nn <- solve(Qhat, y/sigma_e^2)
    tmp <- get.nn.pred(loc = loc, kappa = kappa, nu = nu, sigma = sigma, n.nbr = i_m, S = obs.ind)
    Bp <- tmp$B
    mu.nn <- Bp%*%mu.nn
    post_cov <- Bp %*% solve(Qhat, t(Bp))+ tmp$F
    if(type == "unobserved"){
        mu.nn <- mu.nn[idx_pred]
        post_cov <- post_cov[idx_pred, idx_pred] 
    }
    return(list(post_mean = mu.nn, post_cov = post_cov))
}


# Type indicates which locations to obtain the posterior probabilities
# all indicates all locations
# unobserved indicates only the unobserved locations
compute_prob_ests <- function(type = c("unobserved", "all"), loc_full, obs_ind, loc_prob_idx, L_statespace, range, sigma, nu_vec, sigma_e, m_rat, 
m_nngp_fun, m_pca_fun, m_fourier_fun, m_statespace_fun, samples_fourier = 10, seed = 1){
    type <- type[[1]]
    set.seed(seed)
    prob_ests <- list()
    N <- length(loc_full)
    N_pred <- length(loc_prob_idx)
    for(nu in nu_vec){
        print(paste0("nu = ", nu))
        alpha <- nu + 0.5  
        kappa <- sqrt(8*nu)/range
        true_cov <- get_cov_mat(loc = loc_full, m = NULL, method = "true", nu = nu, kappa = kappa, sigma = sigma, samples = NULL, L=NULL)
        z <- matrix(rnorm(N-N_pred), ncol = N-N_pred, nrow = 1)
        L <- chol(true_cov[obs_ind,obs_ind])
        y <- t(L)%*%t(z) + sigma_e * rnorm(N-N_pred)
        post_true <- posterior_constructor_cov(true_cov, y, sigma_e, obs_ind, loc_prob_idx, type = type)
        post_mean_true <- as.vector(post_true$post_mean)
        post_cov_true <- as.matrix(post_true$post_cov)
        
        lb_prob <- min(y) 
        ub_prob <- max(y)

        prob_true <- pmvnorm(lower=lb_prob,upper=ub_prob,mean=post_mean_true,sigma = post_cov_true)
        prob_true <- c(prob_true)
        prob_ests[[as.character(nu)]] <- list()
        prob_ests[[as.character(nu)]][["true"]] <- prob_true

        prob_ests[[as.character(nu)]][["rat"]] <- list()
        for(i_m in m_rat){
            prec_rat <- get_cov_mat(loc = loc_full, m = i_m, method = "rat_markov", nu = nu, kappa = kappa, sigma = sigma, samples = NULL, L=NULL)
            post_rat <- posterior_constructor_prec(prec_rat$Q, y, sigma_e, A = prec_rat$A, obs_ind, loc_prob_idx, type = type)
            post_mean_rat <- as.vector(post_rat$post_mean)
            post_cov_rat <- as.matrix(post_rat$post_cov)
            prob_rat <- pmvnorm(lower=lb_prob,upper=ub_prob,mean=post_mean_rat,sigma = as.matrix(post_cov_rat))
            prob_rat <- c(prob_rat)

            prob_ests[[as.character(nu)]][["rat"]][[as.character(i_m)]] <- prob_rat

            i_nngp <- m_nngp_fun(i_m, nu + 0.5)

            prec_nngp <- get_cov_mat(loc = loc_full[obs_ind], m = i_nngp, method = "nngp", nu = nu, kappa = kappa, sigma = sigma, samples = NULL, L=NULL)
            post_nngp <- posterior_constructor_nngp(prec_mat = prec_nngp, y=y, sigma_e = sigma_e, loc_prob_idx, obs_ind, loc_full, i_nngp, nu, kappa, sigma, type = type)
            post_mean_nngp <- as.vector(post_nngp$post_mean)
            post_cov_nngp <- as.matrix(post_nngp$post_cov)
            prob_nngp <- pmvnorm(lower=lb_prob,upper=ub_prob,mean=post_mean_nngp,sigma = post_cov_nngp)
            prob_nngp <- c(prob_nngp)
            prob_ests[[as.character(nu)]][["nngp"]][[as.character(i_m)]] <- prob_nngp
            
            i_pca <- m_pca_fun(i_m, nu + 0.5)
            i_pca <- min(i_pca, length(obs_ind))
            pca_cov <- get_cov_mat(loc = loc_full, m = i_pca, method = "pca", nu = nu, kappa = kappa, sigma = sigma, samples = NULL, L=NULL)
            post_pca <- posterior_constructor_cov(pca_cov, y, sigma_e, obs_ind, loc_prob_idx, type = type)
            post_mean_pca <- as.vector(post_pca$post_mean)
            post_cov_pca <- as.matrix(post_pca$post_cov)
            prob_pca <- pmvnorm(lower=lb_prob,upper=ub_prob,mean=post_mean_pca,sigma = post_cov_pca)
            prob_pca <- c(prob_pca)
            prob_ests[[as.character(nu)]][["pca"]][[as.character(i_m)]] <- prob_pca

            i_fourier <- m_fourier_fun(i_m, nu + 0.5)
            i_fourier <- min(i_fourier, length(obs_ind))
            fourier_cov <- get_cov_mat(loc = loc_full, m = i_fourier, method = "fourier", nu = nu, kappa = kappa, sigma = sigma, samples = samples, L=NULL)
            post_fourier <- posterior_constructor_cov(fourier_cov, y, sigma_e, obs_ind, loc_prob_idx, type = type)
            post_mean_fourier <- as.vector(post_fourier$post_mean)
            post_cov_fourier <- as.matrix(post_fourier$post_cov)
            prob_fourier <- pmvnorm(lower=lb_prob,upper=ub_prob,mean=post_mean_fourier, sigma = post_cov_fourier)
            prob_fourier <- c(prob_fourier)
            prob_ests[[as.character(nu)]][["fourier"]][[as.character(i_m)]] <- prob_fourier

            i_statespace <- m_statespace_fun(i_m, nu + 0.5)
            statespace_cov <- get_cov_mat(loc = loc_full, m = i_statespace, method = "statespace", nu = nu, kappa = kappa, sigma = sigma, samples = NULL, L=L_statespace)
            post_statespace <- posterior_constructor_cov(statespace_cov, y, sigma_e, obs_ind, loc_prob_idx, type = type)
            post_mean_statespace <- as.vector(post_statespace$post_mean)
            post_cov_statespace <- as.matrix(post_statespace$post_cov)
            prob_statespace <- pmvnorm(lower=lb_prob,upper=ub_prob,mean=post_mean_statespace, sigma = post_cov_statespace)
            prob_statespace <- c(prob_statespace)
            prob_ests[[as.character(nu)]][["statespace"]][[as.character(i_m)]] <- prob_statespace            
        }

    }
    return(prob_ests)
}