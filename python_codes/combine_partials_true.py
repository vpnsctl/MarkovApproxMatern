import h5py
import numpy as np
import os
import glob

def combine_partial_files(n, n_obs, range_value, sigma_e):
    # Create arrays to store all results
    n_nu = 249  # Number of nu values (from 0.01 to 2.49)
    n_rep = 100
    
    sim_data_combined = np.zeros((n_nu, n_rep, n_obs))
    true_mean_combined = np.zeros((n_nu, n_rep, n))
    true_sigma_combined = np.zeros((n_nu, n_rep, n))
    if n != n_obs:
        obs_ind_combined = np.zeros((n_nu, n_rep, n_obs))
    
    # Directory where partial files are stored
    base_dir = "python_codes/partial_results"
    
    # Generate nu values
    nu_values = np.arange(2.49, 0.00, -0.01)[:249]
    
    # Loop through nu values
    for idx, nu in enumerate(nu_values):
        # Construct filename
        filename = f"{base_dir}/simulation_results_n{n}_nobs{n_obs}_range{range_value}_nu{nu:.2f}_sigmae{sigma_e:.2f}.h5"
        
        if not os.path.exists(filename):
            print(f"Warning: File not found: {filename}")
            continue
            
        try:
            with h5py.File(filename, 'r') as f:
                sim_data_combined[idx, :, :] = f['sim_data_nu'][:]
                true_mean_combined[idx, :, :] = f['true_mean_nu'][:]
                true_sigma_combined[idx, :, :] = f['true_sigma_nu'][:]
                if n != n_obs:
                    obs_ind_combined[idx, :, :] = f['obs_ind_nu'][:]
                    
            print(f"Processed file for nu = {nu:.2f}")
            
        except Exception as e:
            print(f"Error processing file {filename}: {str(e)}")
            continue
    
    # Save combined results
    output_filename = f"combined_results_n{n}_nobs{n_obs}_range{range_value}_sigmae{sigma_e:.2f}.h5"
    
    with h5py.File(output_filename, 'w') as f:
        f.create_dataset('sim_data_result', data=sim_data_combined)
        f.create_dataset('true_mean_result', data=true_mean_combined)
        f.create_dataset('true_sigma_result', data=true_sigma_combined)
        f.create_dataset('nu_vec', data=nu_values)
        if n != n_obs:
            f.create_dataset('obs_ind_result', data=obs_ind_combined)
    
    print(f"\nCombined results saved to {output_filename}")
    print(f"Final array shapes:")
    print(f"sim_data_result: {sim_data_combined.shape}")
    print(f"true_mean_result: {true_mean_combined.shape}")
    print(f"true_sigma_result: {true_sigma_combined.shape}")

# Example usage
if __name__ == "__main__":
    # Set your parameters here
    n = 5000
    n_obs = 5000
    range_value = 2
    sigma_e = 0.1
    
    combine_partial_files(n, n_obs, range_value, sigma_e)