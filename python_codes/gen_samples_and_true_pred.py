import tensorflow as tf

tf.random.set_seed(123) # Setting the seed for reproducibility

import time

import tensorflow_probability as tfp   # For special functions like Bessel K and gamma

import math

n = 10000
n_obs = 10000

n_rep = 100

range = 0.5
sigma = 1
sigma_e = 0.1

nu = 0.6

kappa = math.sqrt(8*nu)/range


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

def compute_matern_covariance_toeplitz(n_points, kappa, sigma, nu, sigma_e, ret_dense):
    loc = tf.linspace(0.0, n_points / 100.0, n_points)
    C_row = matern_covariance(loc, kappa=kappa, nu=nu, sigma=sigma)
    C_row = tf.tensor_scatter_nd_update(C_row, [[0]], [C_row[0] + sigma_e**2])
    
    toeplitz_operator = tf.linalg.LinearOperatorToeplitz(row=C_row, col=C_row)
    
    C_full = toeplitz_operator.to_dense()
    if(ret_dense):
        return C_full
    else:
        return toeplitz_operator

def simulate_from_covariance(Sigma, obs_ind, sigma_e, ret_numpy = False):
    try:
        obs_ind = tf.convert_to_tensor(obs_ind, dtype=tf.int32)
        
        Sigma_dense = Sigma.to_dense()
        Sigma_sub = tf.gather(Sigma_dense, obs_ind) 
        Sigma_sub = tf.gather(Sigma_sub, obs_ind, axis=1)
        R = tf.linalg.cholesky(Sigma_sub) # This is the transpose of the Cholesky we obtain in R
        
        random_normal = tf.random.normal(shape=(len(obs_ind),1), dtype=tf.float64)
                                
        X = tf.matmul(R, random_normal)
        
        noise = sigma_e * tf.random.normal(shape=(len(obs_ind),), dtype=tf.float64)
        sim = tf.add(X, noise)
        
        sim = X
        
        if(ret_numpy):
            return sim.numpy()
        else:
            return sim

    except tf.errors.InvalidArgumentError:
        return None

Sigma = compute_matern_covariance_toeplitz(n_points = n, kappa = kappa, sigma = sigma, nu = nu, sigma_e = 0, ret_dense = False)

obs_ind = generate_obs_indices(n=n, n_obs=n_obs)

sim_data = simulate_from_covariance(Sigma, obs_ind, sigma_e = sigma_e)