source("aux_functions/aux_functions_cov.R")
source("aux_functions/predict_methods.R")
source("aux_functions/prediction_computations.R")


N <- 500
n_obs <- 400
L <- 10
range <- 1
sigma <- 1
# nu_vec=seq(0.1, 1.5,  by = 0.1)
nu_vec = c(0.3, 1.2, 2.2)
idx <- (nu_vec + 0.5)%%1 > 1e-10
nu_vec <- nu_vec[idx]
samples = 10
m <- 1:6
# Std. dev. of the measurement error
sigma_e <- 0.01

# m_nngp_fun <- function(m, alpha){
#               if(alpha<1) {
#                 mn <- m
#             } else if (alpha < 2) {
#                 mn <- round(m*(ceil(alpha)+1)^1.5)    
#             } else {
#                 mn <- round(m*(ceil(alpha)+1)^1.5)
#             }
#             return(mn)
# } 



m_nngp_fun <- function(m, alpha){
  mn <- numeric(length(m)) 
    for (i in 1:length(m))
    {  

if(alpha<1) {
                mn[i] <- m[i] - 1
            } else if (alpha < 2) {
                 if(m[i] == 2){
                     mn[i] <-  7  
                 } else if(m[i] == 3) {
                    mn[i] <- 16
                } else if(m[i] == 4) {
                    mn[i] <- 23
                 } else if (m[i] == 5)
                mn[i] = 31
                else{
                  mn[i] = 36
                }
                }

             else {

              if (m[i]== 2)
              {
                mn[i] = 34
                }

      else if (m[i] == 3) {
         mn[i] <- 43
       } else if (m[i] == 4) {
         mn[i] <- 50
       } else if (m[i] == 5){
         mn[i] = 57
       }
       else{
         mn[i] = 63
       }

 
}
           
}
            return(mn)
} 


m_pca_fun <- function(m, alpha){
      return(pmax(9*round((m+1)*sqrt(m) * ceil(alpha)^(3/2)),2))
}
m_fourier_fun <- function(m, alpha){
    return(pmax(9*round((m+1)*sqrt(m) * ceil(alpha)^(3/2)),2))
}
m_statespace_fun <- function(m, alpha){
    return(m+floor(alpha) + 1)
}

t1 = Sys.time()
df_pred_error <- compute_pred_errors(N = N, n_obs = n_obs, range = range, sigma=sigma, nu.vec = nu_vec, 
m.vec = m, m_nngp_fun = m_nngp_fun, m_pca_fun = m_pca_fun, m_fourier_fun = m_fourier_fun, m_statespace_fun = m_statespace_fun,
sigma_e = sigma_e, L = L, seed = 123, print=FALSE)
t2 = Sys.time()

time = time = as.numeric(t2-t1, units = "secs")