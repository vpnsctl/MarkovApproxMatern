source("aux_functions/aux_functions_cov.R")

sample_matern <- function(loc, nu, kappa, sigma, nsim = 1){
    Dloc <- dist(loc)
    Dloc <- as.matrix(Dloc)
    N <- length(loc)
    cov_mat <- rSPDE::matern.covariance(h=Dloc,kappa=kappa,nu=nu,sigma=sigma)
    z <- matrix(rnorm(N * nsim), ncol = N, nrow = nsim)
    L <- chol(cov_mat)
    return(t(L)%*%t(z))
}

# # Example:
# loc <- seq(0,1,by = 0.001)
# nu <- 0.6
# kappa <- 10
# sigma <- 1
# sim <- sample_matern(loc = loc, nu = nu, kappa = kappa, sigma = sigma, nsim = 10000)
# library(rSPDE)
# c.true <- matern.covariance(0.5-loc, kappa=kappa, nu=nu, sigma=sigma)
# plot(loc, c.true,
#      type = "l", ylab = "C(|s-0.5|)", xlab = "s", ylim = c(0, 5),
#      cex.main = 0.8, cex.axis = 0.8, cex.lab = 0.8
#    )
# lines(loc, cov(t(sim))[(length(loc)-1)/2+1,], col = 2)


sample_y <- function(loc, nu, kappa, sigma, sigma_e, seed=123){
    set.seed(seed)
    z <- sample_matern(loc = loc, nu = nu, kappa = kappa, sigma = sigma, nsim = 1)
    return(z + sigma_e * rnorm(length(z)))
}

# # Example:
# nu <- 0.8
# kappa <- 10
# sigma <- 1
# N_full <- 100
# n_loc <- 50
# loc_full <- seq(0,10,length.out = N_full)
# idx_pred <- sample(1:N_full, n_loc)
# idx_obs <- sample(1:N_full, n_loc)
# idx_pred <- setdiff(idx_pred, idx_obs)
# y <- sample_y(loc_full[idx_obs],nu,kappa,sigma,0.1, 1)
# loc_pred <- loc_full[idx_pred]

# loc - observation locations
# loc_pred - locations to obtain predictions, if NULL, the observations locations will be used.
# loc_full - optional, a vector containing locations to extract both obs and pred locations.
# idx_obs - if loc_full is not NULL, the indices of the observation locations.
# idx_pred - if loc_full is not NULL, the indices of the locations to obtain predictions.

true_pred <- function(y, loc = NULL, loc_pred = NULL, loc_full = NULL, idx_obs = NULL, idx_pred = NULL, nu, kappa, sigma, sigma_e){
    if(is.null(loc_full)){
        if(is.null(loc)){
            stop("either loc or loc_full needs to be non-NULL")
        }
        loc_full <- c(loc_pred,loc)
        idx_obs <- (length(loc_pred)+1):length(loc_full)
        if(!is.null(loc_pred)){
            idx_pred <- 1:length(loc_pred)
        } else{
            idx_pred <- idx_obs
        }
    }

    Dloc <- dist(loc_full)
    Dloc <- as.matrix(Dloc)
    cov_mat = rSPDE::matern.covariance(h=Dloc,kappa=kappa,nu=nu,sigma=sigma)
    cov_mat_nugget <-  cov_mat[idx_obs, idx_obs] + Matrix::Diagonal(x=sigma_e^2, n=length(idx_obs))
    cov_mat <- cov_mat[idx_pred,idx_obs]
    d <- solve(cov_mat_nugget, y)
    return(cov_mat%*%d)
}

# # Example:
# post_mean_true <- true_pred(y, loc=loc_full[idx_obs], loc_pred = loc_full[idx_pred], nu = nu, kappa=kappa, sigma=sigma, sigma_e = 0.1)
# post_mean_true2 <- true_pred(y, loc_full = loc_full, idx_obs = idx_obs, idx_pred = idx_pred, nu = nu, kappa=kappa, sigma=sigma, sigma_e = 0.1)

# post_mean_true_obs <- true_pred(y, loc=loc_full[idx_obs], loc_pred = loc_full[idx_obs], nu = nu, kappa=kappa, sigma=sigma, sigma_e = 0.1)



#### Predict error computation

## Predict rational markov (our method) - Compute the predictions from a full vector of locations (typically equally spaced)

pred_rat_markov <- function(y, loc_full, idx_obs, idx_pred, m, nu, kappa, sigma, sigma_e, equally_spaced = FALSE){
pred <- list()
    for(i_m in m){
            r <- rSPDE:::matern.rational.precision(loc = loc_full, order = i_m, nu = nu, kappa = kappa, cumsum = TRUE, ordering = "location", sigma = sigma, type_rational = "brasil", type_interp = "spline")
            
            A_obs <- r$A[idx_obs,]
            A_pred <- r$A[idx_pred,]
            A_mat = t(A_obs)
            Q_xgiveny <-(A_mat%*% (A_obs))/sigma_e^2 + r$Q
            post_y <- (A_mat%*% y)/sigma_e^2
            R <- Matrix::Cholesky(Q_xgiveny, perm = FALSE)         
            mu_xgiveny <- solve(R, post_y, system = "A")
            approx_mean1 <-  A_pred %*% mu_xgiveny
            pred[[as.character(i_m)]] <- approx_mean1
        }
    return(pred)
}



# # Example:
# post_mean_rat <- pred_rat_markov(y, loc_full = loc_full, idx_obs = idx_obs, idx_pred = idx_pred, m =1:6, nu = nu, kappa=kappa, sigma=sigma, sigma_e = 0.1, equally_spaced = TRUE)
# ord_lf <- order(c(loc_full[idx_obs], loc_full[idx_pred]))

# start <- Sys.time()
# post_mean_rat_obs <- pred_rat_markov(y, loc_full = loc_full, idx_obs = idx_obs, idx_pred = idx_obs, m = 1:6, nu = nu, kappa=kappa, sigma=sigma, sigma_e = 0.1, equally_spaced = TRUE)
# end <- Sys.time()
# end - start

## Predict rational markov (our method) - Compute the predictions using LDL decomposition from a full vector of locations (typically equally spaced)

pred_rat_markov_ldl <- function(y, loc_full, idx_obs,idx_pred, m, nu, kappa, sigma, sigma_e, equally_spaced = FALSE){
pred <- list()
    for(i_m in m){
        if(nu < 0.5){
            r <- rSPDE:::matern.rational.ldl(loc_full, order = i_m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp = "spline", equally_spaced = equally_spaced, ordering = "field")
        } else {
            r <- rSPDE:::matern.rational.ldl(loc_full, order = i_m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp = "spline", equally_spaced = equally_spaced, ordering = "field")
        } 
            r$Q <- t(r$L) %*% r$D %*% r$L
                A_obs <- r$A[idx_obs,]
                A_pred <- r$A[idx_pred,]     
            A_mat = t(A_obs)
            Q_xgiveny <-(A_mat%*% (A_obs))/sigma_e^2 + r$Q
            post_y <- (A_mat%*% y)/sigma_e^2
            R <- Matrix::Cholesky(Q_xgiveny, perm=FALSE)
            mu_xgiveny <- solve(R, post_y, system = "A")
            approx_mean1 <-  A_pred %*% mu_xgiveny
            pred[[as.character(i_m)]] <- approx_mean1
        }
    return(pred)
}


# # Example:
# start <- Sys.time()
# post_mean_rat_ldl <- pred_rat_markov_ldl(y, loc_full = loc_full, idx_obs = idx_obs, idx_pred = idx_pred, m = 1:6, nu = nu, kappa=kappa, sigma=sigma, sigma_e = 0.1, equally_spaced = TRUE)
# end <- Sys.time()
# end - start


## Predict rational markov (our method) - Compute the predictions at the observations

pred_rat_markov_obs <- function(y, loc_obs, m, nu, kappa, sigma, sigma_e, sorted = FALSE, test_equally_spaced = FALSE, return_full = FALSE){   
    if(!sorted){
        ord <- order(loc_obs)
        y <- y[ord]
        loc_obs <- loc_obs[ord]
    }

    equally_spaced <- FALSE

    if(test_equally_spaced){
        dif_loc <- diff(loc_obs)
        if(min(dif_loc) == max(dif_loc)){
            equally_spaced <- TRUE
        }
    }
    
    pred <- list()
    if(return_timings){
        pred[["timings"]] <- list()
    }    
    for(i_m in m){
        if(return_timings){
            t1 <- Sys.time()
        }        
            if(nu < 0.5){
                r <- rSPDE:::matern.rational.precision(loc_obs, order = i_m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp = "spline", equally_spaced = equally_spaced, ordering = "field")            
            } else {
                r <- rSPDE:::matern.rational.precision(loc_obs, order = i_m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp = "spline", equally_spaced = equally_spaced, ordering = "field")            
            } 
         if(return_timings){
            t2 <- Sys.time()
        }            
                A_mat = t(r$A)
                Q_xgiveny <-(A_mat%*%r$A)/sigma_e^2 + r$Q
                post_y <- (A_mat%*% y)/sigma_e^2
                R <- Matrix::Cholesky(Q_xgiveny)            
                mu_xgiveny <- solve(R, post_y, system = "A")
                if(!return_full){
                    approx_mean1 <-  r$A %*% mu_xgiveny
                    if(!sorted){
                        approx_mean1[ord] <- approx_mean1
                    }                    
                } else{
                    approx_mean1 <-  mu_xgiveny
                    pred[["A"]][[as.character(i_m)]] <- r$A
                }
                pred[[as.character(i_m)]] <- approx_mean1
                if(return_timings){
                    t3 <- Sys.time()
                    pred[["timings"]][[as.character(i_m)]][["build_Q"]] <- t2 - t1
                    pred[["timings"]][[as.character(i_m)]][["get_pred"]] <- t3 - t2
                }                     
            }
    return(pred)
}

# # Example:
# loc_obs <- loc_full[idx_obs]
# start <- Sys.time()
# post_mean_rat_obs <- pred_rat_markov_obs(y, loc_obs = loc_obs, m = 1:6, nu = nu, kappa=kappa, sigma=sigma, sigma_e = 0.1)
# end <- Sys.time()
# end - start


## Predict rational markov obs (our method) - Compute the predictions at the observations using LDL

pred_rat_markov_ldl_obs <- function(y, loc_obs, m, nu, kappa, sigma, sigma_e, sorted = FALSE, test_equally_spaced = FALSE, equally_spaced = FALSE, return_full = FALSE, return_timings = FALSE){
    if(!sorted){
        ord <- order(loc_obs)
        y <- y[ord]
        loc_obs <- loc_obs[ord]
    }
   

    if(test_equally_spaced){
        dif_loc <- diff(loc_obs)
        if(min(dif_loc) == max(dif_loc)){
            equally_spaced <- TRUE
        }
    }
    
    pred <- list()
    if(return_timings){
        pred[["timings"]] <- list()
    }
    for(i_m in m){
        if(return_timings){
            t1 <- Sys.time()
        }
            if(nu < 0.5){
                r <- rSPDE:::matern.rational.ldl(loc_obs, order = i_m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp = "spline", equally_spaced = equally_spaced, ordering = "field")          
            } else {
                r <- rSPDE:::matern.rational.ldl(loc_obs, order = i_m, nu = nu, kappa = kappa, sigma = sigma, type_rational = "brasil", type_interp = "spline", equally_spaced = equally_spaced, ordering = "field")      
            } 
                r$Q <- t(r$L) %*% r$D %*% r$L
         if(return_timings){
            t2 <- Sys.time()
        }
                A_mat = t(r$A)
                Q_xgiveny <-(A_mat%*%r$A)/sigma_e^2 + r$Q
                post_y <- (A_mat%*% y)/sigma_e^2
                R <- Matrix::Cholesky(Q_xgiveny)            
                mu_xgiveny <- solve(R, post_y, system = "A")
                if(!return_full){
                    approx_mean1 <-  r$A %*% mu_xgiveny
                    if(!sorted){
                        approx_mean1[ord] <- approx_mean1
                    }                    
                } else{
                    approx_mean1 <-  mu_xgiveny
                    pred[["A"]][[as.character(i_m)]] <- r$A                    
                }
                pred[[as.character(i_m)]] <- approx_mean1
                if(return_timings){
                    t3 <- Sys.time()
                    pred[["timings"]][[as.character(i_m)]][["build_Q"]] <- t2 - t1
                    pred[["timings"]][[as.character(i_m)]][["get_pred"]] <- t3 - t2
                }                
            }
    return(pred)
}

# # Example:
# start <- Sys.time()
# post_mean_rat_ldl_obs <- pred_rat_markov_ldl_obs(y, loc_obs = loc_full[idx_obs], m = 1:6, nu = nu, kappa=kappa, sigma=sigma, sigma_e = 0.1)
# end <- Sys.time()
# end - start

## Predict rational markov at any location based on observation predictions (our method) - Compute the predictions at any location using observation predictions

pred_rat_markov_pred <- function(y, loc_obs, loc_pred = NULL, m, nu, kappa, sigma, sigma_e, use_LDL = TRUE, sorted = FALSE, test_equally_spaced = FALSE, equally_spaced = FALSE, return_timings = FALSE){
    ord_obs <- order(loc_obs)
    sorted_loc_obs <- loc_obs[ord_obs]
    y_sorted <- y[ord_obs]

    if(test_equally_spaced){
        dif_loc <- diff(loc_obs)
        if(min(dif_loc) == max(dif_loc)){
            equally_spaced <- TRUE
        }
    }

    if(use_LDL){
        pred_obs <- pred_rat_markov_ldl_obs(y_sorted, loc_obs = sorted_loc_obs, m = m, nu = nu, kappa=kappa, sigma=sigma, sigma_e = sigma_e, sorted = sorted, test_equally_spaced = test_equally_spaced, equally_spaced = equally_spaced, return_full = TRUE, return_timings = return_timings)
    } else{
        pred_obs <- pred_rat_markov_obs(y_sorted, loc_obs = sorted_loc_obs, m = m, nu = nu, kappa=kappa, sigma=sigma, sigma_e = sigma_e, sorted = sorted, test_equally_spaced = test_equally_spaced,  return_full = TRUE, return_timings = return_timings)
    }

    # full_loc_pred <- loc_pred
    # full_pred_obs <- pred_obs

    # idx_pred_obs <- (loc_pred %in% loc_obs)

    # loc_pred <- loc_pred[!idx_pred_obs]

    # print("getting obs predictions")
    # print(end-start)

    # print("creating nbrs matrix")
    start <- Sys.time()
    nbrs_mat <- matrix(nrow = length(loc_pred), ncol=2)
    loc_pred_obs <- rep(FALSE, length(loc_pred))
    loc_pred_obs_loc <- rep(NA, length(loc_pred))
    for(i in 1:length(loc_pred)){
        dists_full <- (loc_pred[i] - sorted_loc_obs)
        dists_tmp <- dists_tmp2 <- dists_full
        dists_tmp[dists_full>0] <- -Inf
        dists_tmp2[dists_full < 0] <- Inf
        if(all(dists_tmp == -Inf)){
            dists_tmp <- dists_tmp2
        }
        nbrs_mat[i,] <- sort(c(which.min(abs(dists_tmp)), which.min(dists_tmp2)))
        if(min(abs(dists_tmp)) == 0 || min(dists_tmp2) == 0){
            loc_pred_obs[i] <- TRUE
            loc_pred_obs_loc[i] <- nbrs_mat[i,1]
        }
    }

    loc_pred <- loc_pred[!loc_pred_obs]
    nbrs_mat <- nbrs_mat[!loc_pred_obs,]

    # different neighbors
    d.nbrs <- (nbrs_mat[,1] - nbrs_mat[,2])
    # Total neighbors
    T.nbrs <- sum(d.nbrs == 0) + 2*sum(d.nbrs != 0)

    cov_coeff <- list()
    A <- list()
    alpha <- nu + 0.5
    exp_1 <- max(floor(alpha)-1,0)

    k <- length(loc_pred)
    n_obs <- length(sorted_loc_obs)

    end <- Sys.time()

    if(return_timings){
        pred_obs[["timings"]][["setup"]] <- end - start
    }
    
    for(i_m in m){
        if(return_timings){
            start <- Sys.time()    
        }

        N <- T.nbrs *((exp_1+1) + i_m * (floor(alpha)+1))^2
        ii <- numeric(N)
        jj <- numeric(N)
        val <- numeric(N)
        ii_A <- numeric(k * (i_m+1))
        jj_A <- numeric(k * (i_m+1))
        val_A <- rep(1,(k * (i_m+1)))

        counter <- 0     
        counter_A <- 0

        coeff <- rSPDE:::interp_rational_coefficients(order = i_m, 
                                              type_rational_approx = "brasil", 
                                              type_interp = "spline", 
                                              alpha = alpha)
        r_coeff <- coeff$r
        p_coeff <- coeff$p
        k_coeff <- coeff$k  

        for(i in 1:k){

            nbrs <- unique(nbrs_mat[i,])
            n.nbr <- length(nbrs)
            po_loc <- c(loc_pred[i],sorted_loc_obs[nbrs])
            ord_po <- order(po_loc)

                if(length(po_loc) == 2){
                    tmp_1 <- rSPDE:::matern.k.joint(po_loc[ord_po][1],po_loc[ord_po][1], kappa = kappa, alpha = alpha)
                    tmp_2 <- rSPDE:::matern.k.joint(po_loc[ord_po][1], po_loc[ord_po][2], kappa = kappa, alpha = alpha)

                    Sigma_po <- k_coeff*sigma^2 * rbind(cbind(tmp_1, 
                                     tmp_2),cbind(
                               t(tmp_2), 
                                     tmp_1))             
                } else{
                    tmp_1 <- rSPDE:::matern.k.joint(po_loc[ord_po][1],po_loc[ord_po][1],kappa,alpha)
                    tmp_2 <- rSPDE:::matern.k.joint(po_loc[ord_po][1],po_loc[ord_po][2],kappa,alpha)
                    tmp_3 <- rSPDE:::matern.k.joint(po_loc[ord_po][1],po_loc[ord_po][3],kappa,alpha)
                    tmp_4 <- rSPDE:::matern.k.joint(po_loc[ord_po][2],po_loc[ord_po][3],kappa,alpha)

                    Sigma_po <- k_coeff*sigma^2 * rbind(cbind(tmp_1, 
                                     tmp_2,
                                     tmp_3),
                               cbind(t(tmp_2),
                                     tmp_1,
                                     tmp_4),
                               cbind(t(tmp_3),
                                     t(tmp_4),
                                     tmp_1))                    
                }
            # }

            for(ii_p in 1:length(p_coeff)){
                if(length(po_loc) == 2){
                    tmp_1 <- rSPDE:::matern.p.joint(po_loc[ord_po][1],po_loc[ord_po][1], kappa = kappa,
                                      p = p_coeff[ii_p], alpha = alpha)
                    tmp_2 <- rSPDE:::matern.p.joint(po_loc[ord_po][1], po_loc[ord_po][2], kappa = kappa,
                                      p = p_coeff[ii_p], alpha = alpha)

                    Sigma_tmp <- r_coeff[ii_p]*sigma^2* rbind(cbind(tmp_1, 
                                     tmp_2),
                               cbind(t(tmp_2), 
                                     tmp_1))                    
                } else{
                    tmp_1 <- rSPDE:::matern.p.joint(po_loc[ord_po][1],po_loc[ord_po][1],kappa,p_coeff[ii_p],alpha)
                    tmp_2 <- rSPDE:::matern.p.joint(po_loc[ord_po][1],po_loc[ord_po][2],kappa,p_coeff[ii_p],alpha)
                    tmp_3 <- rSPDE:::matern.p.joint(po_loc[ord_po][1],po_loc[ord_po][3],kappa,p_coeff[ii_p],alpha)
                    tmp_4 <- rSPDE:::matern.p.joint(po_loc[ord_po][2],po_loc[ord_po][3],kappa,p_coeff[ii_p],alpha)

                    Sigma_tmp <- r_coeff[ii_p]*sigma^2*rbind(cbind(tmp_1, 
                                     tmp_2,
                                     tmp_3),
                               cbind(t(tmp_2),
                                     tmp_1,
                                     tmp_4),
                               cbind(t(tmp_3),
                                     t(tmp_4),
                                     tmp_1))                    
                }
                Sigma_po <- bdiag(Sigma_po, Sigma_tmp)
            }           

            idx_p <- c()
            if(length(nbrs)==1){
                if(ord_po[1] == 1){
                    idx_p <- 1:(1+exp_1)
                    tmp_num <- 1+exp_1 + exp_1+2
                    idx_p <- c(idx_p, tmp_num:(tmp_num + floor(alpha)))               
                } else{
                    idx_p <- (2+exp_1):(2+2*exp_1)
                    tmp_num <- 2+2*exp_1 + floor(alpha)+1 + 1
                    idx_p <- c(idx_p, tmp_num:(tmp_num + floor(alpha)))                    
                }
                if(i_m > 1){
                        for(tmp in 2:i_m){
                            tmp_num <- max(idx_p)+(floor(alpha)+1)+1
                            idx_p <- c(idx_p, tmp_num:(tmp_num + floor(alpha)))
                        }
                }                       
            } else{
                # in this case, ord_po[1] is always 2
                idx_p <- c((2+exp_1):(2+2*exp_1))
                idx_p <- unique(idx_p)
                tmp_num <- 2 + 2*exp_1 + (exp_1+2) + (floor(alpha)+1)
                idx_p <- c(idx_p, tmp_num:(tmp_num + floor(alpha)))
                if(i_m > 1){
                    for(tmp in 2:i_m){
                        tmp_num <- max(idx_p)+2*(floor(alpha)+1)+1
                        idx_p <- c(idx_p, tmp_num:(tmp_num + floor(alpha)))
                    }
                }
            }
            
            Sigma_oo <- Sigma_po[-idx_p, -idx_p]
            Sigma_po <- Sigma_po[idx_p, -idx_p]

            n_terms_i <-  (exp_1 + 1 + i_m * (1+floor(alpha)))*n.nbr

            n.terms <- length(Sigma_po)

            val[counter + 1:n.terms] <- (solve(Sigma_oo,t(Sigma_po)))   

            if(i == 1){
                tmp_idx <- 1
            } else{
                tmp_idx <- (i-1)*(exp_1+1) + 1
            }
            ii[counter +  1:n_terms_i] <- rep(tmp_idx, n_terms_i)
            if(exp_1>0){
                for(tmp in 1:exp_1){
                    ii[counter +  (n_terms_i*tmp+1):((tmp+1)*n_terms_i)] <- rep(tmp_idx+tmp,n_terms_i)
                }
            }
            
            ii[counter + n_terms_i*(exp_1+1) +  1:(n_terms_i*(floor(alpha)+1))] <- rep(i + (exp_1+1)*k + (i-1)*(floor(alpha)) + 0:floor(alpha),each=n_terms_i)            
            if(i_m>1){
                for(tmp in 2:i_m){
                    ii[counter + n_terms_i*(exp_1+1) +  ((tmp-1)*(floor(alpha)+1)*n_terms_i+1):(tmp*n_terms_i*(floor(alpha)+1))] <- rep(i - floor(alpha) + (exp_1+1)*k + k*(tmp-1)*(floor(alpha)+1) + i*(floor(alpha)) + 0:floor(alpha),each=n_terms_i)
                }
            }

            jj_A[counter_A + 1] <- tmp_idx
            ii_A[counter_A + 1] <- i
            for(tmp in 1:i_m){
                if(floor(alpha) == 0){
                    jj_A[counter_A + 1 + tmp] <- i + tmp*k
                } else{
                    jj_A[counter_A + 1 + tmp] <-  i + (exp_1+1)*tmp*k + (tmp-1)*k + floor(alpha)*(i-1)
                }
                    ii_A[counter_A + 1 + tmp] <- i
            }
            
            tmp_jj <- c()
            if(floor(alpha)== 0){
                for(tmp_nb in nbrs){
                    tmp_jj <- c(tmp_jj, (tmp_nb - 1)*(exp_1+1) + 1:(exp_1+1))
                }            
                for(tmp_nb in nbrs){
                    tmp_jj <- c(tmp_jj, n_obs + (tmp_nb - 1)*(floor(alpha)+1) + 1:(floor(alpha)+1))
                }      

                if(i_m > 1){
                    for(tmp_coeff in 2:i_m){
                        for(tmp_nb in nbrs){
                            tmp_jj <- c(tmp_jj, n_obs * tmp_coeff  + (tmp_nb - 1)*(floor(alpha)+1) + 1:(floor(alpha)+1))
                        }     
                    }
                }
            } else{
                for(tmp_nb in nbrs){
                    tmp_jj <- c(tmp_jj, (tmp_nb - 1)*(exp_1+1) + 1:(exp_1+1))
                }            
                for(tmp_nb in nbrs){
                    tmp_jj <- c(tmp_jj, (exp_1+1)*n_obs + (tmp_nb - 1)*(floor(alpha)+1) + 1:(floor(alpha)+1))
                }      

                if(i_m > 1){
                    for(tmp_coeff in 2:i_m){
                        for(tmp_nb in nbrs){
                            tmp_jj <- c(tmp_jj, n_obs * (floor(alpha)+1) * (tmp_coeff-1) + n_obs * (floor(alpha)) + (tmp_nb - 1)*(floor(alpha)+1) + 1:(floor(alpha)+1))
                        }     
                    }
                }
            }

            jj[counter +  1:n.terms] <- tmp_jj
            counter <- counter + n.terms
            counter_A <- counter_A + i_m + 1         
        }      

    cov_coeff[[as.character(i_m)]] <- Matrix::sparseMatrix(i   = ii,
                                j    = jj,
                                x    = val,
                                dims = c( k*(exp_1+1 + i_m * (floor(alpha)+1)), n_obs*(exp_1+1 + i_m * (floor(alpha)+1))))                         
    
    A[[as.character(i_m)]] <- Matrix::sparseMatrix(i   = ii_A,
                                j    = jj_A,
                                x    = val_A,
                                dims = c( k, k*(exp_1+1 + i_m * (floor(alpha)+1))))
        if(return_timings){
            end <- Sys.time()
            pred_obs[["timings"]][[as.character(i_m)]][["build_pred_matrices"]] <- end-start
        }
    }
    for(i_m in m){
        if(return_timings){
            start <- Sys.time()
        }
        pred_obs_full <- rep(0, length(loc_pred_obs))
        pred_obs_full[loc_pred_obs] <- (pred_obs[["A"]][[as.character(i_m)]] %*% pred_obs[[as.character(i_m)]])[loc_pred_obs_loc[loc_pred_obs]]
        pred_obs_full[!loc_pred_obs] <- A[[as.character(i_m)]] %*% cov_coeff[[as.character(i_m)]] %*% pred_obs[[as.character(i_m)]]
        pred_obs[[as.character(i_m)]] <- pred_obs_full
        if(return_timings){
            end <- Sys.time()
            pred_obs[["timings"]][[as.character(i_m)]][["assemble_preds"]] <- end-start
        }
    }

    return(pred_obs)
}

# # Example 

# start <- Sys.time()
#     post_mean_rat_pred <- pred_rat_markov_pred(y, loc_obs = loc_full[idx_obs], loc_pred = loc_full[idx_pred] ,m = 1:2, nu = 0.3, kappa=kappa, sigma=sigma, sigma_e = 0.1)
# end <- Sys.time()

# start <- Sys.time()
#     post_mean_rat <- pred_rat_markov(y, loc_full = loc_full, idx_obs = idx_obs, idx_pred = idx_pred, m = 1:6, nu = 3.2,
#      kappa=kappa, sigma=sigma, sigma_e = 0.1, equally_spaced = FALSE) 
# end <- Sys.time()

# end - start

# end - start


# for(nu in c(0.3, 0.8, 1.2, 1.8, 2.1, 2.8, 3.2)){
#     start <- Sys.time()
#     post_mean_rat_pred <- pred_rat_markov_pred(y, loc_obs = loc_full[idx_obs], loc_pred = loc_full[idx_pred] ,m = 1:6, nu = nu, kappa=kappa, sigma=sigma, sigma_e = 0.1)
#     end <- Sys.time()
#     print("new")
#     print(end - start)
    
#     start <- Sys.time()
#     post_mean_rat <- pred_rat_markov(y, loc_full = loc_full, idx_obs = idx_obs, idx_pred = idx_pred, m = 1:6, nu = nu,
#      kappa=kappa, sigma=sigma, sigma_e = 0.1, equally_spaced = FALSE)
#     end <- Sys.time()
#     print("old")
#     print(end-start)


#      print("nu")
#      print(nu)
#      for(i_m in 1:6){
#         print(paste0("m=",i_m))
#         print(sum(abs(post_mean_rat_pred[[i_m]] - post_mean_rat[[i_m]])))
#      }
# }
# end <- Sys.time()
# end - start



# start <- Sys.time()
# post_mean_rat_pred <- pred_rat_markov_pred(y, loc_obs = loc_full[idx_obs], loc_pred = loc_full[idx_pred],m = 1:6, nu = nu, kappa=kappa, sigma=sigma, sigma_e = 0.1)
# end <- Sys.time()
# end - start


## Predict PCA

pred_PCA <- function(y, loc, loc_pred = NULL, m, nu, kappa, sigma, sigma_e, method = c("standard", "woodbury", "svd"), m_pca_fun){
loc_full <- c(loc_pred,loc)
N <- length(loc_full)
method <- method[[1]]
pred <- list()
D_loc <- as.matrix(dist(loc_full))
cov_mat <- rSPDE::matern.covariance(h=D_loc,kappa=kappa,nu=nu,sigma=sigma)
eigen_cov <- eigen(cov_mat)
rat_m <- m
m <- m_pca_fun(m, nu + 0.5)
count <- 1
for(i_m in m){
    K <- eigen_cov$vec[,1:i_m]    
    D <- Diagonal(i_m,eigen_cov$val[1:i_m]) 
    if(method == "woodbury"){
        K <- K%*%sqrt(D)
        if(!is.null(loc_pred)){
            K_pred <- K[1:length(loc_pred),]
            K_obs <- K[(length(loc_pred)+1):length(loc_full),]
        } else{
            K_pred <- K_obs <- K
        }
        tK_obs <- t(K_obs)
        I_mat_low <- Matrix::Diagonal(x = 1, n = ncol(K_obs))
        I_mat_high <- Matrix::Diagonal(x = 1, n = nrow(K_obs))
        diag_eps_inv <- Matrix::Diagonal(x = sigma_e^(-2), n = nrow(K_obs))
        tKSig <- tK_obs %*% diag_eps_inv
        SigK <- t(tKSig)
        inv_Part <- solve(I_mat_low + tK_obs%*%SigK)
        nugget_part <- inv_Part %*% tKSig
        cov_mat_nugget_inv <- SigK %*% nugget_part
        solve_nugget <- diag_eps_inv - cov_mat_nugget_inv
        post_mean <- tK_obs %*% solve_nugget %*% y
        post_mean <- K_pred%*%post_mean
    } else if(method == "svd"){
        obs_ind <- (length(loc_pred)+1):length(loc_full)       
        cov_KL <- K%*%D%*%t(K)
        cov_KL <- cov_KL[, obs_ind]        
        K <- K[obs_ind, ]
        svd_K <- svd(K%*%sqrt(D))
        cov_KL_svd_U <- cov_KL %*% svd_K$u
        y_new <- t(svd_K$u) %*% y
        prec_nugget <- cov_KL_svd_U %*% Matrix::Diagonal(x = 1/(svd_K$d^2 + sigma_e^2)) 
        post_mean = prec_nugget%*%y_new
        if(!is.null(loc_pred)){
            post_mean <- post_mean[1:length(loc_pred)]
        }        
    } else {
        obs_ind <- (length(loc_pred)+1):length(loc_full)
        Bo <- K[obs_ind,] 
        Q.hat <- solve(D) + t(Bo)%*%Bo/sigma_e^2 
        post_mean <- K%*%solve(Q.hat, t(Bo)%*%y/sigma_e^2) 
        if(!is.null(loc_pred)){
            post_mean <- post_mean[1:length(loc_pred)]
        }      
    }
    pred[[as.character(rat_m[count])]] <- post_mean
    count <- count + 1
}
return(pred)
}

# # Example:
# post_mean_pca <- pred_PCA(y, loc=s, loc_pred, m = 1:6, nu = nu, kappa=kappa, sigma=sigma, sigma_e = 0.1)


# predict NN 

# pred_rat_NN <- function(y, loc, loc_pred = NULL, m, nu, kappa, sigma, sigma_e){
# loc_full <- c(loc_pred, loc)

# pred <- list()     
# rat_m <- m
# m <- get_m(nu = nu, m = m, method = "nngp", type = "prediction")
# count <- 1
# for(i_m in m){
#         # prec_mat <- get.nnQ(Sigma, i_m)
#         prec_mat <- get.nnQ(loc=loc_full,kappa=kappa,nu=nu,sigma=sigma, n.nbr = i_m)
#         I <- Matrix::Diagonal(n = ncol(prec_mat), x = 1)
#         if(!is.null(loc_pred)){
#             A_pred <- I[1:length(loc_pred),]
#             A_obs <- I[(length(loc_pred)+1):length(loc_full),]
#         } else{
#             A <- I
#         }
#         A_mat <- t(A_obs)
#         Q_xgiveny <- (A_mat%*%A_obs)/sigma_e^2 + prec_mat
#         post_y <- (A_mat%*%y)/sigma_e^2
#         R <- Matrix::Cholesky(Q_xgiveny, perm=FALSE)
#         mu_xgiveny <- solve(R, post_y, system = "A")
#         pred[[as.character(rat_m[count])]] <- A_pred%*%mu_xgiveny
#         count <- count + 1
#     }
#     return(pred)
# }

pred_rat_NN <- function(y, loc_full, idx_pred, idx_obs, m, nu, kappa, sigma, sigma_e, m_nngp_fun){
pred <- list()     
rat_m <- m
m <- m_nngp_fun(m, nu + 0.5)
count <- 1
n.obs <- length(idx_obs)
for(i_m in m){
        # prec_mat <- get.nnQ(Sigma, i_m)
        prec_mat <- get.nnQ(loc=loc_full[idx_obs],kappa=kappa,nu=nu,sigma=sigma, n.nbr = i_m)
        Qhat <- prec_mat + Diagonal(n.obs)/sigma_e^2   
        mu.nn <- solve(Qhat, y/sigma_e^2)
        Bp <- get.nn.pred(loc = loc_full, kappa = kappa, nu = nu, sigma = sigma, n.nbr = i_m, S = idx_obs)$B
        mu.nn <- Bp%*%mu.nn        

        pred[[as.character(rat_m[count])]] <- mu.nn[idx_pred]
        count <- count + 1
    }
    return(pred)
}



# # # Example:
# start <- Sys.time()
# post_mean_nn <- pred_rat_NN(y, loc=s, m = 1:6, nu = nu, kappa=kappa, sigma=sigma, sigma_e = 0.1)
# end <- Sys.time()
# end-start


# Predict Fourier

pred_Fourier <- function(y, loc, loc_pred = NULL, m, nu, kappa, sigma, sigma_e,samples = 100, m_fourier_fun){
loc_full <- c(loc_pred, loc)
N <- length(loc_full)
pred <- list()
D_loc <- as.matrix(dist(loc_full))
rat_m <- m
m <- m_fourier_fun(m, nu + 0.5)
count <- 1
for(i_m in m){
    post_mean <- rep(0, length(loc_pred))
    for(jj in 1:samples){
        K <-  ff.comp(m = i_m, kappa = kappa, alpha = nu + 0.5, loc = loc_full) * sigma^2
        if(!is.null(loc_pred)){
            K_pred <- K[1:length(loc_pred),]
            K_obs <- K[(length(loc_pred)+1):length(loc_full),]
        } else{
            K_pred <- K_obs <- K
        }
        tK_obs <- t(K_obs)
        I_mat_low <- Matrix::Diagonal(x = 1, n = ncol(K_obs))
        I_mat_high <- Matrix::Diagonal(x = 1, n = nrow(K_obs))
        diag_eps_inv <- Matrix::Diagonal(x = sigma_e^(-2), n = nrow(K_obs))
        tKSig <- tK_obs %*% diag_eps_inv
        SigK <- t(tKSig)
        inv_Part <- solve(I_mat_low + tK_obs%*%SigK)
        nugget_part <- inv_Part %*% tKSig
        cov_mat_nugget_inv <- SigK %*% nugget_part
        solve_nugget <- diag_eps_inv - cov_mat_nugget_inv
        post_mean_tmp <- tK_obs %*% solve_nugget %*% y
        post_mean <- post_mean + K_pred%*%post_mean_tmp
    }
    pred[[as.character(rat_m[count])]] <- post_mean/samples
    count <- count + 1
}
return(pred)
}

# # Example:
# post_mean_fourier <- pred_Fourier(y, loc=s, m = 1:6, nu = nu, kappa=kappa, sigma_e = 0.1)


# # Predict SS
# # loc in [0,1]

# pred_statespace <- function(y, loc, loc_pred=NULL, loc_full, m, nu, kappa, sigma, sigma_e, L=1, flim = 2, fact = 100){
# loc_full <- c(loc_pred,loc)
# N <- length(loc_full)
# ind = 1 + fact*(0:(N-1))
# h2 = seq(from=0,to=L,length.out=fact*(N-1)+1)
# pred <- list()     
# rat_m <- m
# m <- get_m(nu = nu, m = m, method = "statespace", type = "prediction")
# count <- 1
# for(i_m in m){
#         coeff <- spec.coeff(kappa = kappa,alpha = nu + 0.5,i_m)
#         S1 <- ab2spec(coeff$a,coeff$b,h2, flim = flim)
#         r1 <- S2cov(S1,h2,flim = flim)
#         acf <- r1[ind]
#         acf <- acf * sigma^2
#         cov_mat <- toeplitz(acf, symmetric=TRUE)
#         acf2 <- acf
#         acf2[1] <- acf2[1] + sigma_e^2
#         cov_mat_nugget <-  toeplitz(as.vector(acf2), symmetric=TRUE)
#         d <- solve(cov_mat_nugget, y)
#         if(!is.null(loc_pred)){
#             cov_mat <- cov_mat[1:length(loc_pred),]
#         }
#         pred[[as.character(rat_m[count])]] <- cov_mat%*%d
#         count <- count + 1
#     }
#     return(pred)
# }

# # Example:
# post_mean_ss <- pred_statespace(y, loc=s, m = 1:6, nu = nu, kappa=kappa, sigma_e = 0.1)


# Predict SS
# loc in [0,1]

pred_statespace_idx <- function(y, idx_obs, idx_pred, loc_full, m, nu, kappa, sigma, sigma_e, L=1, flim = 2, fact = 100, m_statespace_fun){
N <- length(loc_full)
ind = 1 + fact*(0:(N-1))
h2 = seq(from=0,to=L,length.out=fact*(N-1)+1)
pred <- list()     
rat_m <- m
m <- m_statespace_fun(m, nu + 0.5)
count <- 1
for(i_m in m){
        coeff <- spec.coeff(kappa = kappa,alpha = nu + 0.5,i_m)
        S1 <- ab2spec(coeff$a,coeff$b,h2, flim = flim)
        r1 <- S2cov(S1,h2,flim = flim)
        acf <- r1[ind]
        acf <- acf * sigma^2
        cov_mat <- toeplitz(acf, symmetric=TRUE)
        acf2 <- acf
        acf2[1] <- acf2[1] + sigma_e^2
        cov_mat_nugget <-  toeplitz(as.vector(acf2), symmetric=TRUE)
        cov_mat_nugget <- cov_mat_nugget[idx_obs,idx_obs]
        d <- solve(cov_mat_nugget, y)
        cov_mat <- cov_mat[idx_pred, idx_obs]
        pred[[as.character(rat_m[count])]] <- cov_mat%*%d
        count <- count + 1
    }
    return(pred)
}



# Generate complete table
# L is the length of the interval

compute_pred_errors <- function(N, n_obs, range, sigma, nu.vec, m.vec, sigma_e, L = 1, seed = 123, 
m_nngp_fun, m_pca_fun, m_fourier_fun, m_statespace_fun, method_pca = "standard",
print = FALSE){
    post_mean_true <- list()
    post_mean_rat <- list()
    post_mean_rat_ldl <- list()    
    post_mean_PCA <- list()
    post_mean_nnGP <- list()
    post_mean_fourier <- list()
    post_mean_statespace <- list()

    #Rational

    for(ii in 1:length(N)){
        # loc_full <- seq(0,L, length.out = N_full)
        # idx_pred <- sample(1:N_full, n_loc)
        # idx_obs <- sample(1:N_full, n_loc)
        # loc_pred <- loc_full[idx_pred]
        # loc <- loc_full[idx_obs]

        n_loc <- N[ii]
        n_loc_obs <- n_obs[ii]
        loc <- seq(0, L, length.out = n_loc)
        loc_full <- loc_pred <- loc
        obs.ind <- idx_obs <- sort(sample(1:n_loc)[1:n_loc_obs])
        idx_pred <- 1:n_loc
        loc <- loc[idx_obs]

        post_mean_true[[as.character(n_loc)]] <- list()
        post_mean_rat[[as.character(n_loc)]] <- list()
        post_mean_rat_ldl[[as.character(n_loc)]] <- list()        
        post_mean_PCA[[as.character(n_loc)]] <- list()
        post_mean_nnGP[[as.character(n_loc)]] <- list()
        post_mean_fourier[[as.character(n_loc)]] <- list()
        post_mean_statespace[[as.character(n_loc)]] <- list()

        for(i in 1:length(nu.vec)) {
            cat(i/length(nu.vec)," ")
            nu <- nu.vec[i]
            alpha <- nu + 0.5  
            kappa <- sqrt(8*nu)/range
            y <- sample_y(loc = loc,nu = nu,kappa = kappa ,sigma = sigma, sigma_e = sigma_e, seed = seed)
            
            if(print){
                message("Starting true posterior")
            }
            post_mean_true[[as.character(n_loc)]][[as.character(nu)]] <- true_pred(y=y, loc=loc, loc_pred = loc_pred, nu=nu, kappa=kappa, sigma=sigma, sigma_e=sigma_e)
            if(print){
                message("Starting rational posterior")
            }            

            # post_mean_rat[[as.character(n_loc)]][[as.character(nu)]] <- pred_rat_markov(y=y, loc=loc, loc_pred = loc_pred, m=m.vec, nu=nu, kappa=kappa, sigma=sigma, sigma_e=sigma_e, equally_spaced = TRUE)
            post_mean_rat[[as.character(n_loc)]][[as.character(nu)]] <- pred_rat_markov(y=y, loc_full=loc_full, idx_obs=idx_obs, idx_pred = idx_pred, m=m.vec, nu=nu, kappa=kappa, sigma=sigma, sigma_e=sigma_e, equally_spaced = TRUE)            
            
            # if(print){
            #     message("Starting rational ldl posterior")
            # }            
            # post_mean_rat_ldl[[as.character(n_loc)]][[as.character(nu)]] <- pred_rat_markov(y=y, loc=loc, loc_pred = loc_pred, m=m.vec, nu=nu, kappa=kappa, sigma=sigma, sigma_e=sigma_e, equally_spaced = TRUE) 
            # post_mean_rat_ldl[[as.character(n_loc)]][[as.character(nu)]] <- pred_rat_markov_ldl(y=y, loc=loc_full, idx_pred = idx_pred,idx_obs=idx_obs, m=m.vec, nu=nu, kappa=kappa, sigma=sigma, sigma_e=sigma_e, equally_spaced = TRUE)      
            # post_mean_rat_ldl[[as.character(n_loc)]][[as.character(nu)]] <- post_mean_rat_ldl[[as.character(n_loc)]][[as.character(nu)]]           
            if(print){
                message("Starting PCA posterior")
            }
            post_mean_PCA[[as.character(n_loc)]][[as.character(nu)]] <- pred_PCA(y=y, loc=loc, loc_pred = loc_pred, m=m.vec, nu=nu, kappa=kappa, sigma=sigma, sigma_e = sigma_e, m_pca_fun = m_pca_fun, method = method_pca)  
            if(print){
                message("Starting nnGP posterior")
            }
            # post_mean_nnGP[[as.character(n_loc)]][[as.character(nu)]] <- tryCatch(pred_rat_NN(y = y, loc_pred = loc_pred, loc=loc, m=m.vec, nu=nu, kappa=kappa, sigma=sigma, sigma_e=sigma_e), error = function(e){NULL})
            post_mean_nnGP[[as.character(n_loc)]][[as.character(nu)]] <- tryCatch(pred_rat_NN(y = y, loc_full = loc_full, idx_obs=idx_obs, idx_pred = idx_pred, m=m.vec, nu=nu, kappa=kappa, sigma=sigma, sigma_e=sigma_e, m_nngp_fun = m_nngp_fun), error = function(e){NULL})            
            if(print){
                message("Starting Fourier posterior")
            }
            post_mean_fourier[[as.character(n_loc)]][[as.character(nu)]] <- pred_Fourier(y=y, loc=loc, loc_pred = loc_pred, m=m.vec, nu=nu, kappa=kappa, sigma = sigma, sigma_e=sigma_e,samples = 100, m_fourier_fun = m_fourier_fun)
            if(print){
                message("Starting state-space posterior")
            }
            post_mean_statespace[[as.character(n_loc)]][[as.character(nu)]] <- pred_statespace_idx(y=y, idx_obs=idx_obs, idx_pred = idx_pred, loc_full=loc_full, m=m.vec, nu=nu, kappa=kappa, sigma = sigma, sigma_e=sigma_e, flim = 2, L = L, fact = 100, m_statespace_fun = m_statespace_fun)
        }
    }

    df_pred <- data.frame(Method = "Rational", nu = nu.vec[1], Norm = "l2", N = N[1], m = m.vec)
    tmp_true <- post_mean_true[[as.character(N[[1]])]][[as.character(nu[[1]])]]
    error_l2_tmp <- unlist(lapply(post_mean_rat[[as.character(N[1])]][[as.character(nu.vec[1])]], function(y_hat){base::norm(drop(y_hat)-drop(tmp_true),"2")}))
    error_max_tmp <- unlist(lapply(post_mean_rat[[as.character(N[1])]][[as.character(nu.vec[1])]], function(y_hat){max(abs(as.vector(y_hat)-tmp_true))}))    
    df_pred[["Error"]] <- error_l2_tmp
    df_tmp <- data.frame(Method = "Rational", nu = nu.vec[1], Norm = "max", N = N[1], m = m.vec)
    df_tmp[["Error"]] <- error_max_tmp
    df_pred <- rbind(df_pred, df_tmp)
    if(length(nu.vec)>1){
        for(j in 2:length(nu.vec)){
                df_tmp <- data.frame(Method = "Rational", nu = nu.vec[j], Norm = "l2", N = N[1], m = m.vec)
                tmp_true <- post_mean_true[[as.character(N[1])]][[as.character(nu.vec[j])]]
                error_l2_tmp <- unlist(lapply(post_mean_rat[[as.character(N[1])]][[as.character(nu.vec[j])]], function(y_hat){base::norm(drop(y_hat)-drop(tmp_true),"2")}))
                error_max_tmp <- unlist(lapply(post_mean_rat[[as.character(N[1])]][[as.character(nu.vec[j])]], function(y_hat){max(abs(as.vector(y_hat)-tmp_true))}))    
                df_tmp[["Error"]] <- error_l2_tmp
                df_pred <- rbind(df_pred, df_tmp)
                df_tmp <- data.frame(Method = "Rational", nu = nu.vec[j], Norm = "max", N = N[1], m = m.vec)
                df_tmp[["Error"]] <- error_max_tmp
                df_pred <- rbind(df_pred, df_tmp)            
        }
    }
    if(length(N)>1){
        for(i in 2:length(N)){
            for(j in 1:length(nu.vec)){
                    df_tmp <- data.frame(Method = "Rational", nu = nu.vec[j], Norm = "l2", N = N[i], m = m.vec)
                    tmp_true <- post_mean_true[[as.character(N[i])]][[as.character(nu.vec[j])]]
                    error_l2_tmp <- unlist(lapply(post_mean_rat[[as.character(N[i])]][[as.character(nu.vec[j])]], function(y_hat){base::norm(drop(y_hat)-drop(tmp_true),"2")}))
                    error_max_tmp <- unlist(lapply(post_mean_rat[[as.character(N[i])]][[as.character(nu.vec[j])]], function(y_hat){max(abs(as.vector(y_hat)-tmp_true))}))    
                    df_tmp[["Error"]] <- error_l2_tmp
                    df_pred <- rbind(df_pred, df_tmp)
                    df_tmp <- data.frame(Method = "Rational", nu = nu.vec[j], Norm = "max", N = N[i], m = m.vec)
                    df_tmp[["Error"]] <- error_max_tmp
                    df_pred <- rbind(df_pred, df_tmp)            
            }
        }
    }

    # ## LDL
    # for(i in 1:length(N)){
    #     for(j in 1:length(nu.vec)){
    #             df_tmp <- data.frame(Method = "Rational_LDL", nu = nu.vec[j], Norm = "l2", N = N[i], m = m.vec)
    #             tmp_true <- post_mean_true[[as.character(N[i])]][[as.character(nu.vec[j])]]
    #             error_l2_tmp <- unlist(lapply(post_mean_rat_ldl[[as.character(N[i])]][[as.character(nu.vec[j])]], function(y_hat){base::norm(drop(y_hat)-drop(tmp_true),"2")}))
    #             error_max_tmp <- unlist(lapply(post_mean_rat_ldl[[as.character(N[i])]][[as.character(nu.vec[j])]], function(y_hat){max(abs(as.vector(y_hat)-tmp_true))}))    
    #             df_tmp[["Error"]] <- error_l2_tmp
    #             df_pred <- rbind(df_pred, df_tmp)
    #             df_tmp <- data.frame(Method = "Rational_LDL", nu = nu.vec[j], Norm = "max", N = N[i], m = m.vec)
    #             df_tmp[["Error"]] <- error_max_tmp
    #             df_pred <- rbind(df_pred, df_tmp)            
    #     }
    # }    

    ## PCA
    for(i in 1:length(N)){
        for(j in 1:length(nu.vec)){
                df_tmp <- data.frame(Method = "PCA", nu = nu.vec[j], Norm = "l2", N = N[i], m = m.vec)
                tmp_true <- post_mean_true[[as.character(N[i])]][[as.character(nu.vec[j])]]
                error_l2_tmp <- unlist(lapply(post_mean_PCA[[as.character(N[i])]][[as.character(nu.vec[j])]], function(y_hat){base::norm(drop(y_hat)-drop(tmp_true),"2")}))
                error_max_tmp <- unlist(lapply(post_mean_PCA[[as.character(N[i])]][[as.character(nu.vec[j])]], function(y_hat){max(abs(as.vector(y_hat)-tmp_true))}))    
                df_tmp[["Error"]] <- error_l2_tmp
                df_pred <- rbind(df_pred, df_tmp)
                df_tmp <- data.frame(Method = "PCA", nu = nu.vec[j], Norm = "max", N = N[i], m = m.vec)
                df_tmp[["Error"]] <- error_max_tmp
                df_pred <- rbind(df_pred, df_tmp)            
        }
    }

    ## nnGP
    for(i in 1:length(N)){
        for(j in 1:length(nu.vec)){
                df_tmp <- data.frame(Method = "nnGP", nu = nu.vec[j], Norm = "l2", N = N[i], m = m.vec)
                tmp_true <- post_mean_true[[as.character(N[i])]][[as.character(nu.vec[j])]]
                error_l2_tmp <- unlist(lapply(post_mean_nnGP[[as.character(N[i])]][[as.character(nu.vec[j])]], function(y_hat){base::norm(drop(y_hat)-drop(tmp_true),"2")}))
                error_max_tmp <- unlist(lapply(post_mean_nnGP[[as.character(N[i])]][[as.character(nu.vec[j])]], function(y_hat){max(abs(as.vector(y_hat)-tmp_true))}))    
                df_tmp[["Error"]] <- error_l2_tmp
                df_pred <- rbind(df_pred, df_tmp)
                df_tmp <- data.frame(Method = "nnGP", nu = nu.vec[j], Norm = "max", N = N[i], m = m.vec)
                df_tmp[["Error"]] <- error_max_tmp
                df_pred <- rbind(df_pred, df_tmp)            
        }
    }

    ## State-space
    for(i in 1:length(N)){
        for(j in 1:length(nu.vec)){
                df_tmp <- data.frame(Method = "Fourier", nu = nu.vec[j], Norm = "l2", N = N[i], m = m.vec)
                tmp_true <- post_mean_true[[as.character(N[i])]][[as.character(nu.vec[j])]]
                error_l2_tmp <- unlist(lapply(post_mean_statespace[[as.character(N[i])]][[as.character(nu.vec[j])]], function(y_hat){base::norm(drop(y_hat)-drop(tmp_true),"2")}))
                error_max_tmp <- unlist(lapply(post_mean_statespace[[as.character(N[i])]][[as.character(nu.vec[j])]], function(y_hat){max(abs(as.vector(y_hat)-tmp_true))}))    
                df_tmp[["Error"]] <- error_l2_tmp
                df_pred <- rbind(df_pred, df_tmp)
                df_tmp <- data.frame(Method = "Fourier", nu = nu.vec[j], Norm = "max", N = N[i], m = m.vec)
                df_tmp[["Error"]] <- error_max_tmp
                df_pred <- rbind(df_pred, df_tmp)            
        }
    }

    ## Fourier
    for(i in 1:length(N)){
        for(j in 1:length(nu.vec)){
                df_tmp <- data.frame(Method = "State-Space", nu = nu.vec[j], Norm = "l2", N = N[i], m = m.vec)
                tmp_true <- post_mean_true[[as.character(N[i])]][[as.character(nu.vec[j])]]
                error_l2_tmp <- unlist(lapply(post_mean_fourier[[as.character(N[i])]][[as.character(nu.vec[j])]], function(y_hat){base::norm(drop(y_hat)-drop(tmp_true),"2")}))
                error_max_tmp <- unlist(lapply(post_mean_fourier[[as.character(N[i])]][[as.character(nu.vec[j])]], function(y_hat){max(abs(as.vector(y_hat)-tmp_true))}))    
                df_tmp[["Error"]] <- error_l2_tmp
                df_pred <- rbind(df_pred, df_tmp)
                df_tmp <- data.frame(Method = "State-Space", nu = nu.vec[j], Norm = "max", N = N[i], m = m.vec)
                df_tmp[["Error"]] <- error_max_tmp
                df_pred <- rbind(df_pred, df_tmp)            
        }
    }
    df_pred[["m"]] <- as.character(df_pred[["m"]])
    return(df_pred)
}




# Function to plot the L2 and Linfinity distances between the covariance functions
# If methods is null, all methods will be considered
# If n_loc is null, the first value of n_loc will be considered. n_loc must be a character
# n_loc is the number of locations

library(ggplot2)

plot_pred <- function(df_pred, n_loc = NULL, norm = c("l2", "max"), methods = NULL, logscale = TRUE){
    if(is.null(n_loc)){
        n_loc <- df_pred[["N"]][1]
    }
    n_loc <- as.character(n_loc)
    df_pred <- df_pred |> dplyr::filter(N == n_loc)
    if(is.null(methods)){
        methods <- unique(df_pred[["Method"]])
    }
    norm <- norm[[1]]
    if(logscale){
        return(df_pred |> dplyr::filter(Norm == norm, Method %in% methods) |> ggplot2::ggplot() + 
            geom_line(aes(x = nu, y = Error, col = m,linetype=Method)) + scale_y_log10())
    } else{
        return(df_pred |> dplyr::filter(Norm == norm, Method %in% methods) |> ggplot2::ggplot() + 
            geom_line(aes(x = nu, y = Error, col = m,linetype=Method)))
    }
}
