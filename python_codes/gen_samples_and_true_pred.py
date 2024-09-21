import tensorflow as tf

tf.random.set_seed(123) # Setting the seed for reproducibility

import time

import tensorflow_probability as tfp   # For special functions like Bessel K and gamma

import math

import numpy as np

n = 10000 # 10000, 5000
n_obs = 10000 # 10000, 5000

n_rep = 100

range = 0.5 # range = 0.5 , 1 and 2
sigma = 1
sigma_e = 0.1

nu_vec = tf.range(0.01, 2.50, 0.01)

import numpy as np

def compute_kappa(nu, range_value):
    kappa = np.sqrt(8 * nu) / range_value
    return kappa

def generate_obs_indices(n, n_obs):
    indices = tf.range(0, n, dtype = tf.int32)
    indices = tf.random.shuffle(indices)
    indices = tf.sort(indices[:n_obs])
    
    return indices


def matern_covariance(h, kappa, nu, sigma):
    h = tf.cast(h, tf.float64)
    kappa = tf.cast(kappa, tf.float64)
    nu = tf.cast(nu, tf.float64)
    sigma = tf.cast(sigma, tf.float64)

    if nu == 0.5:
        C = sigma**2 * tf.exp(-kappa * tf.abs(h))
    else:
        # Using exp(lgamma(nu)) to calculate gamma(nu)
        C = (sigma**2 / (2**(nu - 1) * tf.exp(tf.math.lgamma(nu)))) * \
            (kappa * tf.abs(h))**nu * tfp.math.bessel_kve(nu, kappa * tf.abs(h))/tf.exp(kappa * tf.abs(h))
        C = tf.where(h == 0, sigma**2, C)  # Ensuring the value at h == 0 is sigma^2
    return C

def compute_matern_covariance_toeplitz(n_points, kappa, sigma, nu, sigma_e, ret_operator = False):
    loc = tf.linspace(0.0, n_points / 100.0, n_points)
    Sigma_row = matern_covariance(loc, kappa=kappa, nu=nu, sigma=sigma)
    Sigma_row = tf.tensor_scatter_nd_update(Sigma_row, [[0]], [Sigma_row[0] + sigma_e**2])
    
    if(ret_operator):
        return tf.linalg.LinearOperatorToeplitz(row=Sigma_row, col=Sigma_row)
    else:
        return Sigma_row

def simulate_from_covariance(Sigma_row, obs_ind, sigma_e, ret_numpy = False):
    try:
        obs_ind = tf.convert_to_tensor(obs_ind, dtype=tf.int32)
        
        Sigma = tf.linalg.LinearOperatorToeplitz(row=Sigma_row, col=Sigma_row).to_dense()
        Sigma_sub = tf.gather(Sigma, obs_ind) 
        Sigma_sub = tf.gather(Sigma_sub, obs_ind, axis=1)
        R = tf.linalg.cholesky(Sigma_sub) # This is the transpose of the Cholesky we obtain in R
        
        random_normal = tf.random.normal(shape=(len(obs_ind),1), dtype=tf.float64)
                                
        X = tf.matmul(R, random_normal)
        
        noise = sigma_e * tf.random.normal(shape=(len(obs_ind),1), dtype=tf.float64)
        sim = tf.add(X, noise)
        
        if(ret_numpy):
            return sim.numpy()
        else:
            return sim

    except tf.errors.InvalidArgumentError:
        return None


import tensorflow as tf

def compute_true_mu(Sigma_row, obs_ind, sigma_e, Y):
    Sigma = tf.linalg.LinearOperatorToeplitz(row=Sigma_row, col=Sigma_row).to_dense()
    Sigma = tf.gather(Sigma, obs_ind, axis=1)  # Make sure Sigma is a dense matrix
    Sigma_row_sigma_e = tf.tensor_scatter_nd_update(Sigma_row, [[0]], [Sigma_row[0] + sigma_e**2])
    Sigma_hat = tf.linalg.LinearOperatorToeplitz(row=Sigma_row_sigma_e, col=Sigma_row_sigma_e).to_dense()
    Sigma_hat = tf.gather(tf.gather(Sigma_hat, obs_ind), obs_ind, axis = 1)
    
    # Solve the linear system: Sigma_hat * mu = Y
    Sigma_hat_inv_Y = tf.linalg.solve(Sigma_hat, Y)  # Equivalent to solve(Sigma.hat, Y)
    
    # Compute mu: Sigma[,obs.ind] %*% solve(Sigma.hat, Y)
    mu = tf.matmul(Sigma, Sigma_hat_inv_Y)
    
    return mu


sim_data_result = np.zeros((tf.size(nu_vec), n_rep, n_obs))
true_mean_result = np.zeros((tf.size(nu_vec), n_rep, n))

# Sigma_row = compute_matern_covariance_toeplitz(n_points = n, kappa = kappa, sigma = sigma, nu = nu, sigma_e = 0, ret_operator = False)

# obs_ind = generate_obs_indices(n=n, n_obs=n_obs)

# sim_data = simulate_from_covariance(Sigma_row, obs_ind, sigma_e = sigma_e)

# mu = compute_true_mu(Sigma_row = Sigma_row, obs_ind = obs_ind, sigma_e = sigma_e, Y = sim_data)
