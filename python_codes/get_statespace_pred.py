import tensorflow as tf
import numpy as np
import math


@tf.function
def spec_coeff(kappa, alpha, n, nm=0):
    A = tf.exp(tf.math.lgamma(alpha)) * tf.sqrt(4 * math.pi) * kappa ** (2 * (alpha - 0.5)) / (2 * math.pi * tf.exp(tf.math.lgamma(alpha - 0.5)))

    ca = tf.cast(tf.math.ceil(alpha), tf.int32)
    k_values = tf.range(0, ca + n + 1, dtype=tf.float64)  
    
    k_factorial = tf.map_fn(lambda x: tf.exp(tf.math.lgamma(x + 1.0)), k_values)
    kappa_power = kappa ** (2 * (alpha + n) - 2 * k_values)

    prod_terms_range = tf.range(0, ca + n + 1, dtype=tf.float64)
    prod_terms_matrix = tf.range(0, ca + n + 1, dtype=tf.float64)[:, tf.newaxis] + tf.range(0, -tf.shape(prod_terms_range)[0], -1, dtype=tf.float64)
    
    mask = tf.cast(tf.sequence_mask(tf.range(0, ca + n + 1, dtype=tf.int32), maxlen=ca + n), dtype=tf.float64)
    
    masked_prod_matrix = (alpha + n + prod_terms_matrix) * mask
    prod_terms = tf.reduce_prod(masked_prod_matrix, axis=1)

    a = tf.where(k_values == 0, kappa_power / k_factorial, kappa_power / k_factorial * prod_terms)

    i_values = tf.range(0, n - nm + 1, dtype=tf.float64)
    n_choose_i = tf.map_fn(lambda i: tf.math.betainc(i, n - nm - i, 0.5), i_values)
    
    kappa_power_b = kappa ** (2 * (n - nm - i_values))
    b = n_choose_i * kappa_power_b
    
    return {'a': a, 'b': b * A}

def ab2spec(a, b, x, flim=2):
    a = tf.convert_to_tensor(a, dtype=tf.float64)
    b = tf.convert_to_tensor(b, dtype=tf.float64)
    x = tf.convert_to_tensor(x, dtype=tf.float64)
    
    nx = tf.size(x)
    x_max = x[-1]
    n = flim * nx - 1
    
    x_extended = np.linspace(0.0, x_max * flim, n, dtype='float64')
    x_extended = tf.convert_to_tensor(loc)
    
    ind = tf.range(0, n, dtype=tf.float64)
    dx = x_extended[1] - x_extended[0]
    
    ds = 2.0 * np.pi / (tf.convert_to_tensor(n, tf.float64) * dx)
    c = -(tf.convert_to_tensor(n, tf.float64) / 2.0) * ds
    s = c + ind * ds
    
    Snum = tf.zeros_like(s, dtype=tf.float64)
    Sden = tf.zeros_like(s, dtype=tf.float64)
    
    powers_num = tf.range(0, tf.size(b), dtype=tf.float64)
    s_powers_num = tf.pow(tf.expand_dims(s, axis=1), 2.0 * powers_num)  # Expand for broadcasting
    Snum = tf.reduce_sum(tf.expand_dims(b, axis=0) * s_powers_num, axis=1)
    
    powers_den = tf.range(0, tf.size(a), dtype=tf.float64)
    s_powers_den = tf.pow(tf.expand_dims(s, axis=1), 2.0 * powers_den)  # Expand for broadcasting
    Sden = tf.reduce_sum(tf.expand_dims(a, axis=0) * s_powers_den, axis=1)
    
    return Snum / Sden


def S2cov(S, x, flim=2):
    S = tf.convert_to_tensor(S, dtype=tf.complex128)  # S should be complex for FFT operations
    x = tf.convert_to_tensor(x, dtype=tf.float64)
    
    nx = tf.size(x)
    x_max = x[-1]
    n = flim * nx - 1
    
    x_extended = np.linspace(0.0, x_max * flim, n, dtype='float64')
    x_extended = tf.convert_to_tensor(loc)
    
    ind = tf.range(0, n, dtype=tf.float64)
    dx = x_extended[1] - x_extended[0]
    
    ds = 2.0 * np.pi / (tf.cast(n, tf.float64) * dx)
    c = -(tf.cast(n, tf.float64) / 2.0) * ds
    s = c + ind * ds
    
    s_min = tf.reduce_min(s)
    x1 = x[0]  
    fact_s = tf.exp(-1j * (s - s_min) * tf.cast(x1, tf.complex128))  
    
    phi_fft = tf.signal.fft(fact_s * S) 
    
    exp_term = tf.exp(-1j * tf.cast(c, tf.complex128) * tf.cast(x, tf.complex128))
    C = tf.math.real(ds * exp_term * phi_fft[:nx])  
    
    return C

def m_ss_fun(m, alpha):
    mn = tf.math.maximum(1.0, m - tf.floor(alpha))
    return mn

def statespace_prediction(kappa, nu, sigma, sigma_e, loc, Y, obs_ind, n, mn):
    ind = 1 + 100 * tf.range(0, n, dtype=tf.int32)
    h2 = np.linspace(0.0, tf.reduce_max(loc), 100 * (n - 1) + 1, dtype='float64')
    h2 = tf.convert_to_tensor(h2)
    
    a, b = spec_coeff(kappa=kappa, alpha=nu + 0.5, n=mn)
    
    S1 = ab2spec(a, b, h2, flim=2)
    r1 = S2cov(S1, h2, flim=2)
    
    acf = tf.gather(r1, ind) * sigma ** 2
    
    cov_mat = tf.signal.toeplitz(acf)

    acf_nugget = tf.identity(acf)
    acf_nugget = tf.tensor_scatter_nd_update(acf_nugget, [[0]], [acf_nugget[0] + sigma_e ** 2])
    
    toeplitz_op_nugget = tf.linalg.LinearOperatorToeplitz(acf_nugget)
    cov_mat_nugget = toeplitz_op_nugget.to_dense()
    cov_mat_nugget = tf.gather(cov_mat_nugget, obs_ind, axis=0)
    cov_mat_nugget = tf.gather(cov_mat_nugget, obs_ind, axis=1)
    
    d = tf.linalg.solve(cov_mat_nugget, Y)
    
    cov_mat_obs = tf.gather(cov_mat, obs_ind, axis=1)
    mu_ss = tf.matmul(cov_mat_obs, d)
    
    return mu_ss    

def get_statespace_errors(n, n_obs, range_val, n_rep, sigma, sigma_e, folder_to_save):
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
    all_err_ss = []

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
        err_ss = np.zeros((1, len(m_vec)))

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

            for j, m in enumerate(m_vec):
                mn = m_ss_fun(m, alpha, n, n_obs)
                mu_ss = statespace_prediction(kappa_val, nu, sigma, sigma_e, loc, Y, obs_ind, n, mn)
                
                loc_diff = loc[1] - loc[0]
                err_ss[0, j] += tf.sqrt(loc_diff * tf.reduce_sum(tf.square(mu - mu_ss))) / n_rep 

            
        print(f"Time: {time.time() - time1}")
        print("Statespace Error :", err_ss[0])
        all_err_ss.append((nu, err_ss[0])) 

        # Save partial results
        partial_save_path = os.path.join(folder_to_save, "pred_tables", f"{n}_{n_obs}", f"range_{range_val}", "statespace")
        os.makedirs(partial_save_path, exist_ok=True)
        
        np.savez(os.path.join(partial_save_path, f"partial_results_nu_{nu:.2f}.npz"), errors=err_ss[0])

    # Save final results
    final_save_path = os.path.join(folder_to_save, "pred_tables")
    os.makedirs(final_save_path, exist_ok=True)
    
    final_results = {"nu": [], "errors": []}
    for nu, errors in all_err_ss:  
        final_results["nu"].append(nu)
        final_results["errors"].append(errors)

    # Convert to DataFrame for saving
    final_df = pd.DataFrame(final_results)
    final_df.to_pickle(os.path.join(final_save_path, f"res_{n}_{n_obs}_range{range_val}_statespace.RDS"))

    return final_df

n = 10000
n_obs = 10000
range_val = 2
sigma = 1.0
sigma_e = 0.1
folder_to_save = os.getcwd()
n_rep = 100

get_statespace_errors(n, n_obs, range_val, n_rep, sigma, sigma_e, folder_to_save)
