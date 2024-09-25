import tensorflow as tf
import numpy as np
import h5py
import time
import tensorflow_probability as tfp
import os

# Set the seed for reproducibility
tf.random.set_seed(123)

# Define your parameters
n = 10000
n_obs = 10000
n_rep = 100
range_value = 2
sigma = 1
sigma_e = 0.1

if(n == n_obs and range_value == 2):
    jitter = 0 # To stabilize Cholesky (needed for high nu, but we will add for all of them for simplicity, as they can be considered measurement error)
else:
    jitter = 0

# Reversed order for nu_vec
nu_vec = tf.range(2.49, 0.01, -0.01, dtype = tf.float64)

# Function to compute kappa
def compute_kappa(nu, range_value):
    kappa = np.sqrt(8 * nu) / range_value
    return kappa

# Function to generate observation indices
def generate_obs_indices(n, n_obs):
    indices = tf.range(0, n, dtype=tf.int32)
    indices = tf.random.shuffle(indices)
    indices = tf.sort(indices[:n_obs])
    return indices

# Matern covariance function
def matern_covariance(h, kappa, nu, sigma):
    h = tf.convert_to_tensor(h, tf.float64)
    kappa = tf.convert_to_tensor(kappa, tf.float64)
    nu = tf.convert_to_tensor(nu, tf.float64)
    sigma = tf.convert_to_tensor(sigma, tf.float64)

    if nu == 0.5:
        C = sigma**2 * tf.exp(-kappa * tf.abs(h))
    else:
        C = (sigma**2 / (2**(nu - 1) * tf.exp(tf.math.lgamma(nu)))) * \
            (kappa * tf.abs(h))**nu * tfp.math.bessel_kve(nu, kappa * tf.abs(h)) / tf.exp(kappa * tf.abs(h))
        C = tf.where(h == 0, sigma**2, C)
    return C

# Function to compute the covariance matrix
def compute_matern_covariance_toeplitz(n_points, kappa, sigma, nu, sigma_e, ret_operator=False):
    loc = np.linspace(0.0, n_points / 100.0, n_points, dtype='float64')
    loc = tf.convert_to_tensor(loc)
    Sigma_row = matern_covariance(loc, kappa=kappa, nu=nu, sigma=sigma)
    Sigma_row = tf.tensor_scatter_nd_update(Sigma_row, [[0]], [Sigma_row[0] + sigma_e**2])
    
    if ret_operator:
        return tf.linalg.LinearOperatorToeplitz(row=Sigma_row, col=Sigma_row)
    else:
        return Sigma_row

# Simulate data from covariance
def simulate_from_covariance(Sigma_row, obs_ind, sigma_e, ret_numpy=False, jitter = 1e-5, R = None):
    obs_ind = tf.convert_to_tensor(obs_ind, dtype=tf.int32)
    
    if R is None:
        Sigma = tf.linalg.LinearOperatorToeplitz(row=Sigma_row, col=Sigma_row).to_dense()
        Sigma_sub = tf.gather(Sigma, obs_ind)
        Sigma_sub = tf.gather(Sigma_sub, obs_ind, axis=1)
        Sigma_sub =  Sigma_sub + tf.eye(Sigma_sub.shape[0], dtype=Sigma_sub.dtype) * jitter
        R = tf.linalg.cholesky(Sigma_sub)

    random_normal = tf.random.normal(shape=(len(obs_ind), 1), dtype=tf.float64)
    X = tf.matmul(R, random_normal)

    noise = sigma_e * tf.random.normal(shape=(len(obs_ind), 1), dtype=tf.float64)
    sim = tf.add(X, noise)
    
    if ret_numpy:
        return sim.numpy()
    else:
        return sim

# Compute the true mean
def compute_true_mu(Sigma_row, obs_ind, sigma_e, Y):
    Sigma = tf.linalg.LinearOperatorToeplitz(row=Sigma_row, col=Sigma_row).to_dense()
    Sigma = tf.gather(Sigma, obs_ind, axis=1)
    Sigma_row_sigma_e = tf.tensor_scatter_nd_update(Sigma_row, [[0]], [Sigma_row[0] + sigma_e**2])
    Sigma_hat = tf.linalg.LinearOperatorToeplitz(row=Sigma_row_sigma_e, col=Sigma_row_sigma_e).to_dense()
    Sigma_hat = tf.gather(tf.gather(Sigma_hat, obs_ind), obs_ind, axis=1)
    
    Sigma_hat_inv_Y = tf.linalg.solve(Sigma_hat, Y)
    mu = tf.matmul(Sigma, Sigma_hat_inv_Y)
    
    return mu

# Initialize result arrays
sim_data_result = np.zeros((tf.size(nu_vec), n_rep, n_obs))
true_mean_result = np.zeros((tf.size(nu_vec), n_rep, n))

sim_data_nu = np.zeros((n_rep, n_obs))
true_mean_nu = np.zeros((n_rep, n))

if(n != n_obs):
    obs_ind_result = np.zeros((tf.size(nu_vec), n_rep, n_obs))
    obs_ind_nu = np.zeros((n_rep, n_obs))

# Create directories for saving results
output_dir = 'python_codes/partial_results'
os.makedirs(output_dir, exist_ok=True)

# Start timing the computations
start_time = time.time()

# Loop over nu values
for idx, nu in enumerate(nu_vec):
    nu_start_time = time.time()

    # Define filename for current nu
    filename = f"{output_dir}/simulation_results_n{n}_nobs{n_obs}_range{range_value}_nu{nu.numpy():.2f}.h5"

    # Check if the file exists
    if os.path.exists(filename):
        with h5py.File(filename, 'r') as f:
            sim_data_result[idx, :, :] = f['sim_data_nu'][:]
            true_mean_result[idx, :, :] = f['true_mean_nu'][:]
            if n != n_obs:
                obs_ind_result[idx, :, :] = f['obs_ind_nu'][:]
        print(f"Loaded data from {filename}")
    else:
        # Compute kappa and Sigma_row since file does not exist
        kappa = compute_kappa(nu, range_value=range_value)
        Sigma_row = compute_matern_covariance_toeplitz(n_points=n, kappa=kappa, sigma=sigma, nu=nu, sigma_e=0, ret_operator=False)
        
        if n == n_obs:
            Sigma = tf.linalg.LinearOperatorToeplitz(row=Sigma_row, col=Sigma_row).to_dense()
            Sigma = Sigma + tf.eye(Sigma.shape[0], dtype=Sigma.dtype) * jitter
            R = tf.linalg.cholesky(Sigma)
        else:
            R = None

        for i in range(n_rep):
            obs_ind = generate_obs_indices(n=n, n_obs=n_obs)
            sim_data = simulate_from_covariance(Sigma_row, obs_ind, sigma_e=sigma_e, R=R)
            mu = compute_true_mu(Sigma_row=Sigma_row, obs_ind=obs_ind, sigma_e=sigma_e, Y=sim_data)

            sim_data_nu[i, :] = sim_data_result[idx, i, :] = sim_data.numpy().flatten()
            true_mean_nu[i, :] = true_mean_result[idx, i, :] = mu.numpy().flatten()
            if n != n_obs:
                obs_ind_nu[i, :] = obs_ind_result[idx, i, :] = obs_ind.numpy().flatten()

        # Save computed results to HDF5 file
        with h5py.File(filename, 'w') as f:
            f.create_dataset('sim_data_nu', data=sim_data_nu)
            f.create_dataset('true_mean_nu', data=true_mean_nu)
            if n != n_obs:
                f.create_dataset('obs_ind_nu', data=obs_ind_nu)

        print(f"Computed and saved data to {filename}")

    nu_elapsed_time = time.time() - nu_start_time
    print(f"Time for nu = {nu:.2f}: {nu_elapsed_time:.2f} seconds")

total_time = time.time() - start_time
print(f"Total computation time: {total_time:.2f} seconds")

# Save overall results to HDF5 file
filename = f'simulation_results_n{n}_nobs{n_obs}_range{range_value}.h5'

with h5py.File(filename, 'w') as f:
    f.create_dataset('sim_data_result', data=sim_data_result)
    f.create_dataset('true_mean_result', data=true_mean_result)
    f.create_dataset('nu_vec', data=nu_vec.numpy())
    if n != n_obs:
        f.create_dataset('obs_ind_result', data=obs_ind_result)

print(f"Results saved to {filename}")
