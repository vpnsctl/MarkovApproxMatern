import h5py
import numpy as np
import tensorflow as tf
import tensorflow_probability as tfp
import os
import time
import pandas as pd

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

def compute_matern_covariance_toeplitz(loc, kappa, sigma, nu, sigma_e, ret_operator=False):
    Sigma_row = matern_covariance(loc, kappa=kappa, nu=nu, sigma=sigma)
    Sigma_row = tf.tensor_scatter_nd_update(Sigma_row, [[0]], [Sigma_row[0] + sigma_e**2])
    
    if ret_operator:
        return tf.linalg.LinearOperatorToeplitz(row=Sigma_row, col=Sigma_row)
    else:
        return Sigma_row

def compute_eigen_cov(Sigma):
    eigenvalues, eigenvectors = tf.linalg.eigh(Sigma)
    eigenvectors = tf.reverse(eigenvectors, axis=[1])
    eigenvalues = tf.reverse(eigenvalues, axis=[0])
    
    eigen_cov = {
        "val": eigenvalues,
        "vec": eigenvectors
    }
    return eigen_cov

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

def pca_prediction(Y, obs_ind, eigen_cov, sigma_e, m):
    K = eigen_cov['vec'][:, :m]
    reciprocals_eigenval = tf.math.reciprocal(eigen_cov['val'][:m])
    Di = tf.linalg.diag(reciprocals_eigenval)
    Bo = tf.gather(K, obs_ind , axis=0)
    Q_hat = Di + tf.matmul(tf.transpose(Bo), Bo) / sigma_e**2
    mu_pca = tf.matmul(K, tf.linalg.solve(Q_hat, tf.matmul(tf.transpose(Bo), Y) / sigma_e**2))
    return mu_pca

def get_pca_errors(n, n_obs, range_val, n_rep, sigma, sigma_e, folder_to_save):
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
    all_err_pca = []

    nu_vec = np.arange(2.49, 0.01, -0.01)

    loc = tf.linspace(0.0, n / 100.0, n)

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
        err_pca = np.zeros((1, len(m_vec)))

        alpha = nu + 0.5
        kappa_val = np.sqrt(8 * nu) / range_val

        loc = tf.linspace(0.0, n / 100.0, n)
        
        Sigma = compute_matern_covariance_toeplitz(loc=loc, kappa=kappa_val, sigma=sigma, nu=nu, sigma_e=0, ret_operator=True).to_dense()
        
        eigen_cov = compute_eigen_cov(Sigma)
        
        print(f"Getting True pred for nu = {nu:.2f}")

        time1 = time.time()
        for kk in range(n_rep):

            if n == 10000 and n_obs == 5000:
                obs_ind = obs_ind_full[kk]

            Y = tf.gather(full_sim_data_nu[kk], obs_ind)
            Y = tf.reshape(Y, [-1, 1])
            mu = full_true_pred_nu[kk]
            mu = tf.reshape(mu, [-1, 1])

            for j, m in enumerate(m_vec):
                mn = m_pca_fun(m, alpha, n, n_obs)
                mu_pca = pca_prediction(Y, obs_ind, eigen_cov, sigma_e, mn)

                loc_diff = loc[1] - loc[0]
                err_pca[0, j] += np.sqrt(loc_diff * np.sum((mu - mu_pca)**2)) / n_rep

            
        print(f"Time: {time.time() - time1}")
        print("PCA Error :", err_pca[0])
        all_err_pca.append((nu, err_pca[0]))

        # Save partial results
        partial_save_path = os.path.join(folder_to_save, "pred_tables", f"{n}_{n_obs}", f"range_{range_val}", "pca")
        os.makedirs(partial_save_path, exist_ok=True)
        
        np.savez(os.path.join(partial_save_path, f"partial_results_nu_{nu:.2f}.npz"), errors=err_pca[0])

    # Save final results
    final_save_path = os.path.join(folder_to_save, "pred_tables")
    os.makedirs(final_save_path, exist_ok=True)
    
    final_results = {"nu": [], "errors": []}
    for nu, errors in all_err_pca:
        final_results["nu"].append(nu)
        final_results["errors"].append(errors)

    # Convert to DataFrame for saving
    final_df = pd.DataFrame(final_results)
    final_df.to_pickle(os.path.join(final_save_path, f"res_{n}_{n_obs}_range{range_val}_pca.RDS"))

    return final_df

n = 10000
n_obs = 5000
range_val = 1
sigma = 1.0
sigma_e = 0.1
folder_to_save = os.getcwd()
n_rep = 100

get_pca_errors(n, n_obs, range_val, n_rep, sigma, sigma_e, folder_to_save)
