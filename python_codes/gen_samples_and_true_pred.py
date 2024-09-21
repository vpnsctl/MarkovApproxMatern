import tensorflow as tf
import numpy as np
import h5py
import time
import tensorflow_probability as tfp

tf.random.set_seed(123)  # Setting the seed for reproducibility

n = 5000  # 10000, 5000
n_obs = 5000  # 10000, 5000
n_rep = 100
range_value = 0.5  # range = 0.5, 1 and 2
sigma = 1
sigma_e = 0.1

# Reversed order for nu_vec
nu_vec = tf.range(2.50, 0.01, -0.01)

def compute_kappa(nu, range_value):
    kappa = np.sqrt(8 * nu) / range_value
    return kappa

def generate_obs_indices(n, n_obs):
    indices = tf.range(0, n, dtype=tf.int32)
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
        C = (sigma**2 / (2**(nu - 1) * tf.exp(tf.math.lgamma(nu)))) * \
            (kappa * tf.abs(h))**nu * tfp.math.bessel_kve(nu, kappa * tf.abs(h)) / tf.exp(kappa * tf.abs(h))
        C = tf.where(h == 0, sigma**2, C)  # Ensuring the value at h == 0 is sigma^2
    return C

def compute_matern_covariance_toeplitz(n_points, kappa, sigma, nu, sigma_e, ret_operator=False):
    loc = tf.linspace(0.0, n_points / 100.0, n_points)
    Sigma_row = matern_covariance(loc, kappa=kappa, nu=nu, sigma=sigma)
    Sigma_row = tf.tensor_scatter_nd_update(Sigma_row, [[0]], [Sigma_row[0] + sigma_e**2])
    
    if ret_operator:
        return tf.linalg.LinearOperatorToeplitz(row=Sigma_row, col=Sigma_row)
    else:
        return Sigma_row

def simulate_from_covariance(Sigma_row, obs_ind, sigma_e, ret_numpy=False):
    obs_ind = tf.convert_to_tensor(obs_ind, dtype=tf.int32)
    
    Sigma = tf.linalg.LinearOperatorToeplitz(row=Sigma_row, col=Sigma_row).to_dense()
    Sigma_sub = tf.gather(Sigma, obs_ind)
    Sigma_sub = tf.gather(Sigma_sub, obs_ind, axis=1)
    
    R = tfp.experimental.linalg.simple_robustified_cholesky(Sigma_sub)

    random_normal = tf.random.normal(shape=(len(obs_ind), 1), dtype=tf.float64)
    X = tf.matmul(R, random_normal)

    noise = sigma_e * tf.random.normal(shape=(len(obs_ind), 1), dtype=tf.float64)
    sim = tf.add(X, noise)
    
    if ret_numpy:
        return sim.numpy()
    else:
        return sim

def compute_true_mu(Sigma_row, obs_ind, sigma_e, Y):
    Sigma = tf.linalg.LinearOperatorToeplitz(row=Sigma_row, col=Sigma_row).to_dense()
    Sigma = tf.gather(Sigma, obs_ind, axis=1)  # Make sure Sigma is a dense matrix
    Sigma_row_sigma_e = tf.tensor_scatter_nd_update(Sigma_row, [[0]], [Sigma_row[0] + sigma_e**2])
    Sigma_hat = tf.linalg.LinearOperatorToeplitz(row=Sigma_row_sigma_e, col=Sigma_row_sigma_e).to_dense()
    Sigma_hat = tf.gather(tf.gather(Sigma_hat, obs_ind), obs_ind, axis=1)
    
    Sigma_hat_inv_Y = tf.linalg.solve(Sigma_hat, Y)  # Equivalent to solve(Sigma.hat, Y)
    mu = tf.matmul(Sigma, Sigma_hat_inv_Y)
    
    return mu

sim_data_result = np.zeros((tf.size(nu_vec), n_rep, n_obs))
true_mean_result = np.zeros((tf.size(nu_vec), n_rep, n))

start_time = time.time()

for idx, nu in enumerate(nu_vec):
    nu_start_time = time.time()

    kappa = compute_kappa(nu, range_value=range_value)
    Sigma_row = compute_matern_covariance_toeplitz(n_points=n, kappa=kappa, sigma=sigma, nu=nu, sigma_e=0, ret_operator=False)

    for i in range(n_rep):
        obs_ind = generate_obs_indices(n=n, n_obs=n_obs)
        sim_data = simulate_from_covariance(Sigma_row, obs_ind, sigma_e=sigma_e)
        mu = compute_true_mu(Sigma_row=Sigma_row, obs_ind=obs_ind, sigma_e=sigma_e, Y=sim_data)

        sim_data_result[idx, i, :] = sim_data.numpy().flatten()
        true_mean_result[idx, i, :] = mu.numpy().flatten()

    nu_elapsed_time = time.time() - nu_start_time
    print(f"Time for nu = {nu:.2f}: {nu_elapsed_time:.2f} seconds")

total_time = time.time() - start_time
print(f"Total computation time: {total_time:.2f} seconds")

filename = f'simulation_results_n{n}_nobs{n_obs}_range{range_value}.h5'

with h5py.File(filename, 'w') as f:
    f.create_dataset('sim_data_result', data=sim_data_result)
    f.create_dataset('true_mean_result', data=true_mean_result)
    f.create_dataset('nu_vec', data=nu_vec.numpy())  # Save the nu_vec in the HDF5 file

print(f"Results saved to {filename}")
