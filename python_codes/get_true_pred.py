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
n_obs = 5000
n_rep = 100
range_value = 1
sigma = 1
sigma_e = 0.1
# Reversed order for nu_vec
nu_vec = tf.range(2.49, 0.00, -0.01)

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
    h = tf.cast(h, tf.float64)
    kappa = tf.cast(kappa, tf.float64)
    nu = tf.cast(nu, tf.float64)
    sigma = tf.cast(sigma, tf.float64)

    if nu == 0.5:
        C = sigma**2 * tf.exp(-kappa * tf.abs(h))
    else:
        C = (sigma**2 / (2**(nu - 1) * tf.exp(tf.math.lgamma(nu)))) * \
            (kappa * tf.abs(h))**nu * tfp.math.bessel_kve(nu, kappa * tf.abs(h)) / tf.exp(kappa * tf.abs(h))
        C = tf.where(h == 0, sigma**2, C)
    return C

# Function to compute the covariance matrix
def compute_matern_covariance_toeplitz(n_points, kappa, sigma, nu, sigma_e, ret_operator=False):
    loc = tf.linspace(0.0, n_points / 100.0, n_points)
    Sigma_row = matern_covariance(loc, kappa=kappa, nu=nu, sigma=sigma)
    Sigma_row = tf.tensor_scatter_nd_update(Sigma_row, [[0]], [Sigma_row[0] + sigma_e**2])
    
    if ret_operator:
        return tf.linalg.LinearOperatorToeplitz(row=Sigma_row, col=Sigma_row)
    else:
        return Sigma_row

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
true_mean_result = np.zeros((tf.size(nu_vec), n_rep, n))

true_mean_nu = np.zeros((n_rep, n))

if(n != n_obs):
    obs_ind_result = np.zeros((tf.size(nu_vec), n_rep, n_obs))
    obs_ind_nu = np.zeros((n_rep, n_obs))
    
# Create directories for saving results
output_dir = 'python_codes/partial_results'
os.makedirs(output_dir, exist_ok=True)


directory = 'python_codes/results'
filename = f'simulation_results_n10000_nobs10000_range{range_value}.h5'
file_path = os.path.join(directory, filename)

# Ensure the file exists before attempting to read
if not os.path.exists(file_path):
    raise FileNotFoundError(f"The file {file_path} does not exist.")

# Read the data from the HDF5 file
with h5py.File(file_path, 'r') as f:
    # Assuming 'sim_data_result' is the dataset in the HDF5 file
    sim_data_result = f['sim_data_result'][:]

# Start timing the computations
start_time = time.time()

# Loop over nu values
for idx, nu in enumerate(nu_vec):
    nu_start_time = time.time()

    kappa = compute_kappa(nu, range_value=range_value)
    Sigma_row = compute_matern_covariance_toeplitz(n_points=10000, kappa=kappa, sigma=sigma, nu=nu, sigma_e=0, ret_operator=False)

    for i in range(n_rep):
        if(n == 10000 and n_obs == 5000):
            obs_ind = generate_obs_indices(n=n, n_obs=n_obs)
        elif (n == 5000):
            obs_ind = tf.range(0, n, dtype=tf.int32)
        if(n_obs == 5000):
            sim_data = sim_data_result[idx, i, obs_ind] 
        else:
            sim_data = sim_data_result[idx, i, :]
             
        sim_data = tf.convert_to_tensor(sim_data.flatten(), dtype=tf.float64)
        sim_data = tf.reshape(sim_data, (-1, 1))
        mu = compute_true_mu(Sigma_row=Sigma_row, obs_ind=obs_ind, sigma_e=sigma_e, Y=sim_data)

        true_mean_nu[i, :] = true_mean_result[idx, i, :] = mu.numpy().flatten()
        if(n != n_obs):
            obs_ind_nu[i, :] = obs_ind_result[idx, i, :] = obs_ind.numpy().flatten()        
        
    filename = f"{output_dir}/simulation_results_n{n}_nobs{n_obs}_range{range_value}_nu{nu.numpy():.2f}.h5"
    
    with h5py.File(filename, 'w') as f:
        f.create_dataset('true_mean_nu', data=true_mean_nu)
        if(n != n_obs):
            f.create_dataset('obs_ind_nu', data=obs_ind_nu)

    nu_elapsed_time = time.time() - nu_start_time
    print(f"Time for nu = {nu:.2f}: {nu_elapsed_time:.2f} seconds")

total_time = time.time() - start_time
print(f"Total computation time: {total_time:.2f} seconds")

# Save results to HDF5 file
filename = f'simulation_results_n{n}_nobs{n_obs}_range{range_value}.h5'

with h5py.File(filename, 'w') as f:
    f.create_dataset('true_mean_result', data=true_mean_result)
    f.create_dataset('nu_vec', data=nu_vec.numpy())
    if(n != n_obs):
        f.create_dataset('obs_ind_result', data=obs_ind_result)

print(f"Results saved to {filename}")

