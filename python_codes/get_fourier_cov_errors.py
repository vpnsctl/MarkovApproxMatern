import tensorflow as tf
import numpy as np

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

@tf.function
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

def compute_matern_covariance_toeplitz(n_points, kappa, sigma, nu, sigma_e, ret_operator=False):
    loc = np.linspace(0.0, n_points / 100.0, n_points, dtype='float64')
    loc = tf.convert_to_tensor(loc)
    Sigma_row = matern_covariance(loc, kappa=kappa, nu=nu, sigma=sigma)
    Sigma_row = tf.tensor_scatter_nd_update(Sigma_row, [[0]], [Sigma_row[0] + sigma_e**2])
    
    if ret_operator:
        return tf.linalg.LinearOperatorToeplitz(row=Sigma_row, col=Sigma_row)
    else:
        return Sigma_row

def sample_t_mat(kappa, alpha, m):
    return tf.convert_to_tensor(t.rvs(df=2 * alpha - 1, scale=kappa / tf.sqrt(2 * alpha - 1), size=m), dtype=tf.float64)

@tf.function
def ff_approx(m, loc, w):
    b = tf.random.uniform([m], 0, 2 * np.pi, dtype=tf.float64)
    loc = tf.convert_to_tensor(loc, dtype=tf.float64)  
    ZX = tf.sqrt(2.0) * tf.cos(tf.expand_dims(w, axis=1) * tf.expand_dims(loc, axis=0) + tf.expand_dims(b, axis=1)) / tf.sqrt(tf.cast(m, tf.float64))
    return tf.matmul(tf.transpose(ZX), ZX)

def compute_distances_fourier(N, n_obs, m_vec, nu_vec, range_val, sigma, samples):
    N = N[0]
    n_obs = n_obs[0]
    
    l2_err = np.zeros((len(nu_vec), len(m_vec)), dtype=np.float64)
    sup_err = np.zeros((len(nu_vec), len(m_vec)), dtype=np.float64)
    
    loc = np.linspace(0.0, N/100.0, N, dtype='float64')
    loc = tf.convert_to_tensor(loc)

    print("Fourier")
    print(f"N = {N}")
    print(f"n_obs = {n_obs}")
    print(f"Domain length = {loc.max()}")
    print(f"range = {range_val}")
    
    for i, nu in enumerate(nu_vec):
        print(f"{(i+1)/len(nu_vec):.2f}", end=" ")
        
        alpha = nu + 0.5
        kappa = np.sqrt(8 * nu) / range_val
        Sigma_row = compute_matern_covariance_toeplitz(n_points=n, kappa=kappa, sigma=sigma, nu=nu, sigma_e=0, ret_operator=False)
        Sigma_t = tf.linalg.LinearOperatorToeplitz(row=Sigma_row, col=Sigma_row).to_dense() # True covariance
        
        for j, m in enumerate(m_vec):
            mn = m_pca_fun(m, alpha)  
            Sigma_fou = np.zeros_like(Sigma_t)
            
            for k in range(samples):
                w = sample_t_mat(kappa = kappa, alpha = nu + 0.5, m = mn)
                Sigma_fou += ff_approx(m, loc, w).numpy() * sigma**2
            
            Sigma_fou /= samples
            
            if n_obs < N:
                l2_err[i, j] = np.sqrt(np.sum((Sigma_t[:n_obs, :] - Sigma_fou[:n_obs, :])**2)) * (loc[1] - loc[0])
                sup_err[i, j] = np.max(np.abs(Sigma_t[:n_obs, :] - Sigma_fou[:n_obs, :]))
            else:
                l2_err[i, j] = np.sqrt(np.sum((Sigma_t - Sigma_fou)**2)) * (loc[1] - loc[0])
                sup_err[i, j] = np.max(np.abs(Sigma_t - Sigma_fou))
    
    ret = {
        'L2': l2_err,
        'Linf': sup_err
    }
    
    return ret

import pickle  

# Define the parameters
nu_vec = tf.range(2.49, 0.01, -0.01, dtype = tf.float64)
sigma = 1
m = np.arange(1, 7)
samples_fourier = 10

n = 10000
n_obs = 10000
range_val = 0.5

# Assuming N, n_obs, range_, and m_fourier_fun are already defined in your script
dist_fourier = compute_distances_fourier(N=N, n_obs=n_obs, m_vec=m, nu_vec=nu_vec, range_=range_val, sigma=sigma, samples=samples_fourier)

# Create the filename dynamically using range, N, and n_obs
filename = f"distance_tables/raw_tables/dist_fourier_{n}_{n_obs}_range_{range_val}_calibrated.pkl"

# Save the result to a file using pickle
with open(filename, "wb") as f:
    pickle.dump(dist_fourier, f)

print(f"Results saved to {filename}")
