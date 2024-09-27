import h5py
import numpy as np
import tensorflow as tf
import os
import time
import pandas as pd
import math
from scipy.special import gamma
from scipy.stats import t

def m_pca_fun(m, alpha, n, n_obs):
    if alpha < 1:
        if n == 5000:
            m_vec = [268, 308, 355, 406, 433, 478]
        elif n_obs == 5000:
            m_vec = [381, 448, 533, 588, 654, 727]
        elif n_obs == 10000:
            m_vec = [271, 311, 367, 412, 452, 493]
        else:
            raise ValueError("Not implemented")
    elif alpha < 2:
        if n == 5000:
            m_vec = [380, 473, 561, 651, 708, 776]
        elif n_obs == 5000:
            m_vec = [532, 704, 844, 953, 1065, 1162]
        elif n_obs == 10000:
            m_vec = [372, 493, 591, 672, 751, 821]
        else:
            raise ValueError("Not implemented")
    else:
        if n == 5000:
            m_vec = [611, 810, 945, 1082, 1205, 1325]
        elif n_obs == 5000:
            m_vec = [904, 1202, 1431, 1622, 1808, 1965]
        elif n_obs == 10000:
            m_vec = [640, 848, 1016, 1168, 1299, 1420]
        else:
            raise ValueError("Not implemented")
    return m_vec[m - 1]

def load_hdf5_data(file_path, dataset_name):
    with h5py.File(file_path, 'r') as f:
        data = f[dataset_name][:]
    return data

# def ff_comp(m, kappa, alpha, loc):
#     w = t.rvs(df = 2*alpha - 1, scale = kappa/np.sqrt(2*alpha-1), size = m)
#     b = tf.random.uniform([m], 0, 2 * np.pi, dtype=tf.float64)
#     ZX = np.sqrt(2) * np.cos(w[:, np.newaxis] * loc + b[:, np.newaxis]) / np.sqrt(m)
#     ZX_matrix = tf.convert_to_tensor(ZX.T, dtype=tf.float64)
#     return ZX_matrix

@tf.function
def ff_comp(m, loc, w):
    b = tf.random.uniform([m], 0, 2 * np.pi, dtype=tf.float64)
    loc = tf.convert_to_tensor(loc, dtype=tf.float64)  
    ZX = tf.sqrt(tf.convert_to_tensor(2.0, dtype=tf.float64)) * tf.cos(tf.expand_dims(w, axis=1) * tf.expand_dims(loc, axis=0) + tf.expand_dims(b, axis=1)) / tf.sqrt(tf.cast(m, tf.float64))
    return tf.transpose(ZX)  


# def fourier_prediction(Y, obs_ind, mn, kappa, nu, loc, sigma, sigma_e):
#     sigma = tf.convert_to_tensor(sigma, dtype=tf.float64)
#     sigma_e = tf.convert_to_tensor(sigma_e, dtype=tf.float64)
#     K = ff_comp(m=mn, kappa=kappa, alpha=nu + 0.5, loc=loc) * sigma**2
#     D = tf.linalg.diag(tf.fill([tf.shape(K)[1]], tf.convert_to_tensor(1, dtype=tf.float64)))
#     Bo = tf.gather(K, obs_ind, axis=0)
#     # D is identity, so no need to invert
#     Q_hat = D + tf.matmul(tf.transpose(Bo), Bo) / sigma_e**2
#     mu_fourier = tf.matmul(K, tf.linalg.solve(Q_hat, tf.matmul(tf.transpose(Bo), Y) / sigma_e**2))
#     return mu_fourier

def sample_t_mat(kappa, alpha, m):
    return tf.convert_to_tensor(t.rvs(df=2 * alpha - 1, scale=kappa / tf.sqrt(2 * alpha - 1), size=m), dtype=tf.float64)
    
def fourier_prediction(Y, obs_ind, mn, kappa, nu, loc, sigma, sigma_e):
    sigma = tf.convert_to_tensor(sigma, dtype=tf.float64)
    sigma_e = tf.convert_to_tensor(sigma_e, dtype=tf.float64)
    w = sample_t_mat(kappa = kappa, alpha = nu + 0.5, m = mn)
    K = ff_comp(m=mn, loc=loc, w=w) * tf.square(sigma)
    Bo = tf.gather(K, obs_ind, axis=0)
    Q_hat = tf.linalg.eye(tf.shape(K)[1], dtype=tf.float64) + tf.matmul(tf.transpose(Bo), Bo) / tf.square(sigma_e)
    mu_fourier = tf.matmul(K, tf.linalg.solve(Q_hat, tf.matmul(tf.transpose(Bo), Y) / tf.square(sigma_e)))
    return mu_fourier


def get_fourier_errors(n, n_obs, range_val, n_rep, sigma, sigma_e, samples_fourier, folder_to_save):
    sim_data_name = "sim_data_result"
    true_mean_name = "true_mean_result"
    nu_vec_python = "nu_vec"
    obs_ind_python = "obs_ind_result"

    sim_file_path = f"python_codes/results/simulation_results_n10000_nobs10000_range{range_val}.h5"
    true_pred_file_path = f"python_codes/results/simulation_results_n{n}_nobs{n_obs}_range{range_val}.h5"

    full_sim_data = load_hdf5_data(sim_file_path, sim_data_name)
    full_true_pred = load_hdf5_data(true_pred_file_path, true_mean_name)
    nu_vec_loaded = load_hdf5_data(sim_file_path, nu_vec_python)

    np.random.seed(123)
    m_vec = np.arange(1, 7)
    all_err_fourier = []

    nu_vec = tf.range(2.49, 0.01, -0.01, dtype = tf.float64)

    loc = np.linspace(0.0, n / 100.0, n, dtype='float64')
    loc = tf.convert_to_tensor(loc)

    for nu in nu_vec:
        ind_nu = np.argmin((nu - nu_vec_loaded)**2)

        full_sim_data_nu = full_sim_data[ind_nu,:, :]
        full_true_pred_nu = full_true_pred[ind_nu, :, :]
        
        full_sim_data_nu = tf.convert_to_tensor(full_sim_data_nu)
        full_true_pred_nu = tf.convert_to_tensor(full_true_pred_nu)

        if n == 10000 and n_obs == 5000:
            obs_ind_full = load_hdf5_data(true_pred_file_path, obs_ind_python)[ind_nu, :, :]
            obs_ind_full = tf.convert_to_tensor(obs_ind_full, dtype=tf.int32)
        else:
            obs_ind = tf.range(0, n, dtype=tf.int32)

        if n == 5000:
            full_sim_data_nu = tf.gather(full_sim_data_nu, obs_ind)
            full_true_pred_nu = tf.gather(full_true_pred_nu, obs_ind)

        m_vec = np.arange(1, 7)
        err_fourier = np.zeros((1, len(m_vec)))

        alpha = nu + 0.5
        kappa_val = np.sqrt(8 * nu) / range_val        
        
        print(f"Getting True pred for nu = {nu:.2f}")

        time1 = time.time()
        for kk in range(n_rep):

            if n == 10000 and n_obs == 5000:
                obs_ind = obs_ind_full[kk]

            Y = tf.gather(full_sim_data_nu[kk], obs_ind)
            Y = tf.reshape(Y, [-1, 1])
            mu = full_true_pred_nu[kk]
            mu = tf.reshape(mu, [-1, 1])
            max_attempts_fourier = 10  # number of times to try to sample t distribution (needed for nu = 0.01, when it sometimes gets numerically singular)

            for j, m in enumerate(m_vec):
                mn = m_pca_fun(m, alpha, n, n_obs)
                err_tmp = 0 
                loc_diff = loc[1] - loc[0]

                # Check if partial results exist
                partial_save_path = os.path.join(folder_to_save, "pred_tables", f"{n}_{n_obs}", f"range_{range_val}", "fourier")
                os.makedirs(partial_save_path, exist_ok=True)
                partial_result_file = os.path.join(partial_save_path, f"partial_results_nu_{nu:.2f}_m_{m}.npz")

                if os.path.exists(partial_result_file):
                    print(f"Loading partial result for nu={nu:.2f}, m={m}")
                    partial_result = np.load(partial_result_file)
                    err_tmp = partial_result['errors']
                else:
                    # Compute if file doesn't exist
                    for jj in range(samples_fourier): 
                        attempts = 0
                        while attempts < max_attempts_fourier:
                            try:
                                mu_fourier = fourier_prediction(Y, obs_ind, mn, kappa_val, nu, loc, sigma, sigma_e)
                                break  # Exit the loop if successful
                            except Exception as e:
                                attempts += 1
                                print(f"Attempt {attempts} failed with error: {e}")
                                if attempts == max_attempts_fourier:
                                    print("Max retry attempts reached. Exiting.")
                                    raise  # Re-raise the exception after maximum attempts
                        err_tmp += tf.sqrt(loc_diff * tf.reduce_sum(tf.square(mu - mu_fourier))) / (n_rep * samples_fourier)

                    # Save partial result
                    np.savez(partial_result_file, errors=err_tmp.numpy())
                
                err_fourier[0, j] += err_tmp 

        print(f"Time: {time.time() - time1}")
        print("Fourier Error :", err_fourier[0])
        all_err_fourier.append((nu, err_fourier[0])) 

    # Save final results
    final_save_path = os.path.join(folder_to_save, "pred_tables")
    os.makedirs(final_save_path, exist_ok=True)
    
    final_results = {"nu": [], "errors": []}
    for nu, errors in all_err_fourier:  
        final_results["nu"].append(nu)
        final_results["errors"].append(errors)

    # Convert to DataFrame for saving
    final_df = pd.DataFrame(final_results)
    final_df.to_pickle(os.path.join(final_save_path, f"res_{n}_{n_obs}_range{range_val}_fourier.pkl"))

    return final_df

n = 10000
n_obs = 10000
range_val = 0.5
sigma = 1.0
sigma_e = 0.1
folder_to_save = os.getcwd()
n_rep = 100
samples_fourier = 10

get_fourier_errors(n, n_obs, range_val, n_rep, sigma, sigma_e, samples_fourier, folder_to_save)
