rm(list=ls())
source("aux_functions/calibration_functions.R")
source("aux_functions/aux_functions_cov.R")
source("aux_functions/probability_computations.R")

library(rSPDE)
library(excursions)
library(mvtnorm)

# ------------------------------------------------------------------------------
# 1) Define paths and create folders if needed
# ------------------------------------------------------------------------------
calibration_file_ratnngp <- "calibration_results/calibrated_m_list.RDS"
calibration_file_fem     <- "calibration_results/calibrated_m_list_fem.RDS"
partial_results_folder   <- "prob_tables/partial_results"

if (!dir.exists("calibration_results")) {
  dir.create("calibration_results")
}
if (!dir.exists("prob_tables")) {
  dir.create("prob_tables")
}
if (!dir.exists(partial_results_folder)) {
  dir.create(partial_results_folder)
}

# ------------------------------------------------------------------------------
# 2) Set parameters (adjust as needed)
# ------------------------------------------------------------------------------
range        <- 0.5             # or 0.5, or any other
sigma        <- 1
sigma.e      <- sqrt(0.1) # or 0.1
n            <- c(0, 25, 50, 100, 150, 250, 300, 400, 500,
                  750, 1000, 1250, 1500, 1750, 2000, 2500, 3000)
n            <- rev(n)        # Reverse the vector
n.obs        <- 250
n.rep        <- 10
samples_calibration <- 50
max_it_per_m <- 20
nu           <- 1
alpha        <- nu + 1/2
m_vec        <- 1:6
domain_upper_limit <- 10
use.excursions <- TRUE
coverage     <- 0.9

kappa <- sqrt(8*nu)/range

# ------------------------------------------------------------------------------
# 3) Load or initialize calibration for RAT/NNGP
# ------------------------------------------------------------------------------
calibrated_m_ratnngp <- list()
if(file.exists(calibration_file_ratnngp)) {
  calibrated_m_ratnngp <- readRDS(calibration_file_ratnngp)
  cat("Loaded RAT/NNGP calibrated m from", calibration_file_ratnngp, "\n")
}

# ------------------------------------------------------------------------------
# 4) Load or initialize calibration for FEM
# ------------------------------------------------------------------------------
calibrated_m_fem <- list()
if(file.exists(calibration_file_fem)) {
  calibrated_m_fem <- readRDS(calibration_file_fem)
  cat("Loaded FEM calibrated m from", calibration_file_fem, "\n")
}

# ------------------------------------------------------------------------------
# 5) Prepare location sets for each n
# ------------------------------------------------------------------------------
locs <- lapply(n, function(x) {
  if(x > 0) {
    seq(from = 1, to = domain_upper_limit, length.out = x)
  } else {
    NULL
  }
})

# ------------------------------------------------------------------------------
# 6) Allocate error matrices: RAT, NNGP, FEM
# ------------------------------------------------------------------------------
err_rat  <- matrix(0, nrow = length(m_vec), ncol = length(n))
err_nngp <- matrix(0, nrow = length(m_vec), ncol = length(n))
err_fem  <- matrix(0, nrow = length(m_vec), ncol = length(n))

# ------------------------------------------------------------------------------
# 7) Main replicate loop
# ------------------------------------------------------------------------------
previous_calibration_ratnngp <- NULL
previous_calibration_fem     <- NULL

for (j in seq_len(n.rep)) {
  cat("\n\n============================\n")
  cat("Starting replicate j =", j, "\n")
  cat("============================\n")
  
  # 7.1 Generate observation locations and data
  obs_loc <- sort(runif(n.obs, 0, domain_upper_limit))
  while (min(diff(obs_loc)) < 1e-3) {
    obs_loc <- sort(runif(n.obs, 0, domain_upper_limit))
  }
  
  # Simulate data from the true covariance
  z <- matrix(rnorm(n.obs), ncol = n.obs)
  L <- chol(get_cov_mat(obs_loc, NULL, "true", nu, kappa, sigma, NULL, NULL))
  y <- t(L) %*% t(z) + sigma.e * rnorm(n.obs)
  
  # 7.2 Loop over each n
  for (i in seq_along(n)) {
    # Partial filename includes sigma.e
    partial_file <- file.path(
      partial_results_folder,
      paste0("partial_result_rep", j,
             "_n", n[i],
             "_range", range,
             "_sigmaE", sigma.e,
             "_nu", nu, ".RDS")
    )
    
    # If the file already exists, we can load existing partial results
    if(file.exists(partial_file)) {
      cat("Skipping replicate =", j, ", n =", n[i],
          "- partial file exists.\n")
      partial_saved <- readRDS(partial_file)
      # If the partial file has them, load them
      if ("err_rat"  %in% names(partial_saved)) {
        err_rat[, i]  <- partial_saved$err_rat
      }
      if ("err_nngp" %in% names(partial_saved)) {
        err_nngp[, i] <- partial_saved$err_nngp
      }
      if ("err_fem"  %in% names(partial_saved)) {
        err_fem[, i]  <- partial_saved$err_fem
      }
      # Move on
      next
    }
    
    # 7.3 Build combined location vector and fix ordering
    if(n[i] > 0) {
      loc <- c(obs_loc, locs[[i]])
      obs.ind <- seq_along(loc)
      
      # Sort everything so loc is in ascending order
      tmp <- sort(loc, index.return = TRUE)
      loc <- tmp$x
      obs.ind[tmp$ix] <- seq_along(loc)
      obs.ind <- obs.ind[1:n.obs]
      
      # If any duplicates (or extremely close points), remove them carefully
      if(any(diff(loc) < 1e-3)) {
        ind_remove <- which(diff(loc) < 1e-3)
        ind_remove <- unique(c(ind_remove, ind_remove + 1))
        # Don't remove actual observation points:
        ind_remove <- setdiff(ind_remove, obs.ind)
        loc <- loc[-ind_remove]
        # Re-match obs.ind by closeness
        tolerance <- 1e-6
        obs.ind <- sapply(obs_loc, function(x) {
          which(abs(loc - x) < tolerance)[1]
        })
      }
    } else {
      # If n = 0, only observation locations
      loc     <- obs_loc
      obs.ind <- seq_len(n.obs)
    }
    
    # 7.4 Possibly calibrate RAT/NNGP if needed
    if(is.null(calibrated_m_ratnngp[[as.character(n[i])]])) {
      if(n[i] > 999) {
        cat("Calibrating RAT/NNGP for n =", n[i], "...\n")
        previous_calibration_ratnngp <-
          auto_calibration_nngp_rat(
            n                  = n[i] + n.obs,
            n_obs              = n.obs,
            nu                 = nu,
            range              = range,
            sigma              = sigma,
            sigma_e            = sigma.e,
            samples            = samples_calibration,
            m_rat              = m_vec,
            previous_calibration = previous_calibration_ratnngp,
            max_it_per_m       = max_it_per_m,
            print              = FALSE
          )
        calibrated_m_ratnngp[[as.character(n[i])]] <-
          previous_calibration_ratnngp
        cat("RAT/NNGP calibration:", previous_calibration_ratnngp, "\n")
      } else {
        # If not large n, use calibration from n=1000 (or nearest big)
        calibrated_m_ratnngp[[as.character(n[i])]] <-
          calibrated_m_ratnngp[["1000"]]
      }
    }
    
    # 7.5 Possibly calibrate FEM if needed
    if(is.null(calibrated_m_fem[[as.character(n[i])]])) {
      if(n[i] > 999) {
        cat("Calibrating FEM for n =", n[i], "...\n")
        previous_calibration_fem <-
          auto_calibration_fem_rat(
            n                  = n[i] + n.obs,
            n_obs              = n.obs,
            nu                 = nu,
            range              = range,
            sigma              = sigma,
            sigma_e            = sigma.e,
            samples            = samples_calibration,
            m_rat              = m_vec,
            max_it_per_m       = max_it_per_m,
            print              = FALSE
          )
        calibrated_m_fem[[as.character(n[i])]] <- previous_calibration_fem
        cat("FEM calibration:", previous_calibration_fem, "\n")
      } else {
        # If not large n, use calibration from n=1000 (or nearest big)
        calibrated_m_fem[[as.character(n[i])]] <- calibrated_m_fem[["1000"]]
      }
    }
    
    # 7.6 Compute the "truth"
    cat("Computing true posterior ...\n")
    true_cov <- get_cov_mat(loc = loc, m = NULL, method = "true",
                            nu = nu, kappa = kappa, sigma = sigma,
                            samples = NULL, L = NULL)
    post_true <- posterior_constructor_cov(
        cov_mat  = true_cov, 
        y        = y, 
        sigma_e  = sigma.e, 
        obs_ind  = obs.ind,   
        pred_ind = setdiff(seq_along(loc), obs.ind),  
        type     = "all"
    )
    post_mean_true <- as.vector(post_true$post_mean)
    post_cov_true  <- as.matrix(post_true$post_cov)
    
    # Confidence region bounds
    Q.true <- solve(post_cov_true)
    conf   <- simconf(
      alpha = 1 - coverage,
      mu    = post_mean_true,
      Q     = Q.true,
      vars  = diag(post_cov_true),
      n.iter= 1e5
    )
    lb_prob <- conf$a
    ub_prob <- conf$b
    
    if(use.excursions && n[i] > 999) {
      prob_true <- gaussint(
        a     = lb_prob,
        b     = ub_prob,
        mu    = post_mean_true,
        Q     = Q.true,
        n.iter= 1e5
      )$P
    } else {
      prob_true <- pmvnorm(
        lower = lb_prob,
        upper = ub_prob,
        mean  = post_mean_true,
        sigma = post_cov_true
      )
    }
    cat("Truth probability:", prob_true, "\n")
    
    # 7.7 Evaluate each method for each m
    for(k in seq_along(m_vec)) {
      m <- m_vec[k]
      # ------------------------------------------------------------------------
      # (a) Rational (RAT)
      # ------------------------------------------------------------------------
      cat("Rep =", j, ", n =", n[i], ", RAT, m =", m, "\n")
      Qrat <- rSPDE:::matern.rational.precision(
        loc           = loc,
        order         = m,
        nu            = nu,
        kappa         = kappa,
        sigma         = sigma,
        type_rational = "brasil",
        type_interp   = "spline"
      )
      A_obs   <- Qrat$A[obs.ind, , drop=FALSE]
      A_mat   <- t(A_obs)
      Q_xgiveny <- (A_mat %*% A_obs)/sigma.e^2 + Qrat$Q
      post_y    <- (A_mat %*% y)/sigma.e^2
      R         <- Matrix::Cholesky(Q_xgiveny, perm = FALSE)
      mu_xgiveny<- solve(R, post_y, system = "A")
      mu_rat    <- as.vector(Qrat$A %*% mu_xgiveny)
      Sigma_rat <- as.matrix(Qrat$A %*% solve(Q_xgiveny, t(Qrat$A)))
      
      if(use.excursions && n[i] > 999) {
        prob_rat <- gaussint(
          a     = lb_prob,
          b     = ub_prob,
          mu    = mu_rat,
          Q     = solve(Sigma_rat),
          n.iter= 1e5
        )$P
      } else {
        prob_rat <- pmvnorm(
          lower = lb_prob,
          upper = ub_prob,
          mean  = mu_rat,
          sigma = Sigma_rat
        )
      }
      err_rat[k, i] <- err_rat[k, i] + (prob_rat - prob_true)
      
      # ------------------------------------------------------------------------
      # (b) NNGP
      # ------------------------------------------------------------------------
      # Use calibrated m for NNGP
      mn_nngp <- calibrated_m_ratnngp[[as.character(n[i])]][k]
      cat("Rep =", j, ", n =", n[i], ", NNGP, m =", mn_nngp, "\n")
      prec_nngp <- get_cov_mat(
        loc    = loc[obs.ind],
        m      = mn_nngp,
        method = "nngp",
        nu     = nu,
        kappa  = kappa,
        sigma  = sigma
      )
      post_nngp_obj <- posterior_constructor_nngp(
        prec_mat = prec_nngp,
        y        = y,
        sigma_e  = sigma.e,
        idx_pred = setdiff(seq_along(loc), obs.ind),
        obs.ind  = obs.ind,        # match function argument
        loc      = loc,
        i_m      = mn_nngp,        # match 'i_m' instead of 'm'
        nu       = nu,
        kappa    = kappa,
        sigma    = sigma,
        type     = "full"
      )
      mu_nngp    <- as.vector(post_nngp_obj$post_mean)
      Sigma_nngp <- as.matrix(post_nngp_obj$post_cov)
      
      if(use.excursions && n[i] > 999) {
        prob_nngp <- gaussint(
          a     = lb_prob,
          b     = ub_prob,
          mu    = mu_nngp,
          Q     = solve(Sigma_nngp),
          n.iter= 1e5
        )$P
      } else {
        prob_nngp <- pmvnorm(
          lower = lb_prob,
          upper = ub_prob,
          mean  = mu_nngp,
          sigma = Sigma_nngp
        )
      }
      err_nngp[k, i] <- err_nngp[k, i] + (prob_nngp - prob_true)
      
      # ------------------------------------------------------------------------
      # (c) FEM
      # ------------------------------------------------------------------------
      cat("Rep =", j, ", n =", n[i], ", FEM, m =", m, "\n")
      mesh_fem <- calibrated_m_fem[[as.character(n[i])]][k]
      post_fem_mat <- fem_post_calculations(
        y       = y,
        loc     = loc,
        idx_obs = obs.ind,
        nu      = nu,
        kappa   = kappa,
        sigma   = sigma,
        sigma_e = sigma.e,
        m       = m,         # polynomial order
        mesh_fem= mesh_fem   # mesh size or dof from calibration
      )
      mu_fem    <- post_fem_mat$mu.fem
      Sigma_fem <- post_fem_mat$Sigma.fem
      
      if(use.excursions && n[i] > 999) {
        prob_fem <- gaussint(
          a     = lb_prob,
          b     = ub_prob,
          mu    = mu_fem,
          Q     = solve(Sigma_fem),
          n.iter= 1e5
        )$P
      } else {
        prob_fem <- pmvnorm(
          lower = lb_prob,
          upper = ub_prob,
          mean  = mu_fem,
          sigma = Sigma_fem
        )
      }
      err_fem[k, i] <- err_fem[k, i] + (prob_fem - prob_true)
    }  # end for k in m_vec
    
    # 7.8 Save partial results (RAT, NNGP, FEM) for replicate j, n[i]
    partial_res_list <- list(
      err_rat  = err_rat[, i],
      err_nngp = err_nngp[, i],
      err_fem  = err_fem[, i],
      n        = n[i],
      rep      = j,
      m        = m_vec,
      sigma_e  = sigma.e  # <--- Store the sigma.e used
    )
    saveRDS(partial_res_list, file = partial_file)
  } # end for i in seq_along(n)
} # end for j in seq_len(n.rep)

# ------------------------------------------------------------------------------
# 8) Average the errors across replicates
# ------------------------------------------------------------------------------
err_rat  <- err_rat  / n.rep
err_nngp <- err_nngp / n.rep
err_fem  <- err_fem  / n.rep

# ------------------------------------------------------------------------------
# 9) Save final results (name includes sigma.e)
# ------------------------------------------------------------------------------
res_list <- list(
  err_rat   = err_rat,
  err_nngp  = err_nngp,
  err_fem   = err_fem,
  nu        = nu,
  m_vec     = m_vec,
  n_vec     = n,
  range     = range,
  sigma_e   = sigma.e  # store it in final results, too
)

final_filename <- paste0(
  "prob_tables/rat_nngp_fem_range", range,
  "_sigmaE", sigma.e,
  "_nu", nu, ".RDS"
)
saveRDS(res_list, file = final_filename)

cat("\nDONE! All results saved to:\n", final_filename, "\n")