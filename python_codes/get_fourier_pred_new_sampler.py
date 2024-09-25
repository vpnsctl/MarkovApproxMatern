import h5py
import numpy as np
import tensorflow as tf
import math
import os
import time
import pandas as pd
from scipy.special import gamma
from scipy.stats import cauchy
import concurrent.futures
from scipy.interpolate import UnivariateSpline
from scipy.special import gamma, hyp2f1
from scipy.optimize import brentq
import numpy as np
import tensorflow as tf
from scipy.stats import cauchy
max_workers = 24

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

# CDF obtained from mathematica
# def mat_spec_cdf(x, kappa, alpha):
#     A = (gamma(alpha) * np.sqrt(4 * np.pi) * kappa**(2 * (alpha - 0.5))) / (2 * np.pi * gamma(alpha - 0.5))
#     term1 = x * (1 + (x**2 / kappa**2))**alpha
#     term2 = (x**2 + kappa**2)**(-alpha)
#     hypergeo = hyp2f1(0.5, alpha, 1.5, -x**2 / kappa**2)
#     result = A * term1 * term2 * hypergeo + 0.5
#     return result

def mat_spec(x, kappa, alpha):
    A = gamma(alpha) * np.sqrt(4 * np.pi) * kappa**(2 * (alpha - 0.5)) / (2 * np.pi * gamma(alpha - 0.5))
    return A / ((kappa**2 + x**2)**alpha)

def sample_single_rejection(kappa, alpha):
    c = mat_spec(0.0, kappa=kappa, alpha=alpha) / cauchy.pdf(0.0, loc=0, scale=kappa)
    while True:
        X = cauchy.rvs(loc=0, scale=kappa)
        U1 = np.random.uniform(0, 1)
        fx = cauchy.pdf(X, loc=0, scale=kappa)
        fy = mat_spec(X, kappa=kappa, alpha=alpha)
        Y = c * fx * U1
        
        if Y < fy:
            return X

def sample_mat_rejection_parallel(n, kappa, alpha):
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(sample_single_rejection, kappa, alpha) for _ in range(n)]
        results = [future.result() for future in concurrent.futures.as_completed(futures)]
    return tf.convert_to_tensor(results, dtype=tf.float64)

def ff_comp_rejection(m, kappa, alpha, loc):
    w = sample_mat_rejection_parallel(m, kappa, alpha)
    b = tf.random.uniform([m], 0, 2 * np.pi, dtype=tf.float64)
    ZX = tf.TensorArray(dtype=tf.float64, size=m)
    for i in range(m):
        row_value = tf.sqrt(tf.convert_to_tensor(2.0, tf.float64)) * tf.cos(w[i] * loc + b[i]) / tf.sqrt(tf.convert_to_tensor(m, tf.float64))
        ZX = ZX.write(i, row_value)
    ZX_matrix = ZX.stack()
    return tf.transpose(ZX_matrix)

def mat_spec_cdf(x, kappa, alpha):
    A = (gamma(alpha) * np.sqrt(4 * np.pi) * kappa**(2 * (alpha - 0.5))) / (2 * np.pi * gamma(alpha - 0.5))
    term1 = x * (1 + (x**2 / kappa**2))**alpha
    term2 = (x**2 + kappa**2)**(-alpha)
    hypergeo = np.array([hyp2f1(0.5, alpha, 1.5, -xi**2 / kappa**2) for xi in np.atleast_1d(x)])
    result = A * term1 * term2 * hypergeo + 0.5
    return result

def quant_spec(prob, kappa, alpha):
    def func_to_solve(x):
        return mat_spec_cdf(x, kappa=kappa, alpha=alpha) - prob
    sol = brentq(func_to_solve, 0,200000)
    return sol

def mat_spec_cdf_parallel(x_vals, kappa, alpha):
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(lambda x: mat_spec_cdf(x, kappa, alpha), x_vals))
    return np.array(results)

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
        nu = 0.4
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
        
        print("Tabulating")
        # print("Getting upper bound and lower bound")
        # First step getting a good resolution around the mode:
        ub_mode = quant_spec(prob = 0.75, kappa = kappa_val, alpha = alpha)
        lb_mode = -ub_mode
        factor = 1000
        n_terms_x_vals = min(math.ceil((ub_mode - lb_mode) * factor), 100000)
        x_vals = np.linspace(lb_mode, ub_mode, n_terms_x_vals, dtype="float64")
        
        # Taking this, due to numerical instability for small nu
        ub = 4*quant_spec(prob = 0.9, kappa = kappa_val, alpha = alpha)
        lb = -ub
        
        pub = mat_spec_cdf(ub, kappa = kappa_val, alpha = alpha)
        plb = mat_spec_cdf(lb, kappa = kappa_val, alpha = alpha)
        
        # print(f"Upper bound: {ub:.2f} Lower bound: {lb:.2f}")
        
        n_terms_away_mode = 200
        x_vals_away_modelb = np.linspace(lb, lb_mode, n_terms_away_mode, dtype="float64")
        x_vals_away_modeub = np.linspace(ub_mode, ub, n_terms_away_mode, dtype="float64")
        x_vals = np.concatenate([x_vals, x_vals_away_modelb, x_vals_away_modeub])
        x_vals = np.sort(np.unique(x_vals))
        
        # Compute CDF values on the grid
        # print("Tabulating the spectral density")
        cdf_vals = mat_spec_cdf_parallel(x_vals, kappa_val, alpha)
        print("Done")
        
        spline_interp = UnivariateSpline(x_vals, cdf_vals, s=0) # s=0 means interpolate
        
        
        # alpha > 0.7
        
        def find_root(prob):
            def func_to_solve(x):
                return spline_interp(x) - prob
            sol = brentq(func_to_solve, lb,ub)
            return sol

        def sample_single(lb, ub, plb, pub):
            U1 = np.random.uniform(0, 1)

            if U1 <= plb:
                return lb
            elif U1 >= pub:
                return ub
            else:
                return find_root(U1)  


        def sample_mat_parallel(n,lb,ub,plb,pub):
            with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = [executor.submit(sample_single,lb,ub,plb,pub) for _ in range(n)]
                results = [future.result() for future in concurrent.futures.as_completed(futures)]    
            return tf.convert_to_tensor(results, dtype=tf.float64)        
        
        def ff_comp(m, lb, ub, plb, pub, loc):
            # print("starting parallel sampling")
            w = sample_mat_parallel(m, lb, ub, plb, pub)
            # print("Finished parallel sampling")
            # print("Assembling result")
            b = tf.random.uniform([m], 0, 2 * np.pi, dtype=tf.float64)
            ZX = tf.TensorArray(dtype=tf.float64, size=m)
            for i in range(m):
                row_value = tf.sqrt(tf.convert_to_tensor(2.0, tf.float64)) * tf.cos(w[i] * loc + b[i]) / tf.sqrt(tf.convert_to_tensor(m, tf.float64))
                ZX = ZX.write(i, row_value)
            ZX_matrix = ZX.stack()
            # print("Finished assembling")
            return tf.transpose(ZX_matrix)

        def fourier_prediction(Y, obs_ind, mn, lb, ub, plb, pub, loc, sigma, sigma_e, kappa, alpha):
            sigma = tf.convert_to_tensor(sigma, dtype=tf.float64)
            sigma_e = tf.convert_to_tensor(sigma_e, dtype=tf.float64)
            if(alpha > 0.9):
                # print("starting ff_comp")
                K = ff_comp(m=mn, lb=lb, ub=ub, plb=plb, pub=pub, loc=loc) * sigma**2
                
            else:
                K = ff_comp_rejection(m=mn, kappa=kappa, alpha=alpha, loc=loc) * sigma**2
            D = tf.linalg.diag(tf.fill([tf.shape(K)[1]], tf.convert_to_tensor(1, dtype=tf.float64)))
            Bo = tf.gather(K, obs_ind, axis=0)
            # D is identity, so no need to invert
            # print("mult and solve part")
            Q_hat = D + tf.matmul(tf.transpose(Bo), Bo) / sigma_e**2
            mu_fourier = tf.matmul(K, tf.linalg.solve(Q_hat, tf.matmul(tf.transpose(Bo), Y) / sigma_e**2))
            # print("mult and solve done")                
            return mu_fourier        
        
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
                err_tmp = 0 
                loc_diff = loc[1] - loc[0]
                
                print("Starting Fourier")

                for jj in range(samples_fourier): 
                    print(f"Fourier Sample = {jj}")
                    mu_fourier = fourier_prediction(Y, obs_ind, mn, lb, ub, plb, pub, loc, sigma, sigma_e, kappa_val, alpha)
                    err_tmp += tf.sqrt(loc_diff * tf.reduce_sum(tf.square(mu - mu_fourier))) / (n_rep * samples_fourier)
                err_fourier[0, j] += err_tmp.numpy() 

            
        print(f"Time: {time.time() - time1}")
        print("Fourier Error :", err_fourier[0])
        all_err_fourier.append((nu, err_fourier[0])) 

        # Save partial results
        partial_save_path = os.path.join(folder_to_save, "pred_tables", f"{n}_{n_obs}", f"range_{range_val}", "fourier")
        os.makedirs(partial_save_path, exist_ok=True)
        
        np.savez(os.path.join(partial_save_path, f"partial_results_nu_{nu:.2f}.npz"), errors=err_fourier[0])

    # Save final results
    final_save_path = os.path.join(folder_to_save, "pred_tables")
    os.makedirs(final_save_path, exist_ok=True)
    
    final_results = {"nu": [], "errors": []}
    for nu, errors in all_err_fourier:  
        final_results["nu"].append(nu)
        final_results["errors"].append(errors)

    # Convert to DataFrame for saving
    final_df = pd.DataFrame(final_results)
    final_df.to_pickle(os.path.join(final_save_path, f"res_{n}_{n_obs}_range{range_val}_fourier.RDS"))

    return final_df

n = 10000
n_obs = 10000
range_val = 2
sigma = 1.0
sigma_e = 0.1
folder_to_save = os.getcwd()
n_rep = 100
samples_fourier = 100

get_fourier_errors(n, n_obs, range_val, n_rep, sigma, sigma_e, samples_fourier, folder_to_save)
