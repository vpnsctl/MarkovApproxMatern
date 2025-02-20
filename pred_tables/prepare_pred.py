import numpy as np
import pandas as pd

nu_values = np.arange(0.01, 2.50, 0.01)
all_data = []

for nu in nu_values:
    filename = f"pred_tables//5000_5000/range_2//sigmae_0.32//fourier//partial_results_nu_{nu:.2f}.npz"
    
    try:
        data = np.load(filename)
        
        mu_errors = data['errors_mu']
        sigma_errors = data['errors_sigma']
        
        for m in range(6):
            row = {
                'nu': nu,
                'method': 'fourier',
                'mu_error': mu_errors[m],  
                'sigma_error': sigma_errors[m],  
                'm': m+1
            }
            all_data.append(row)
            
    except Exception as e:
        print(f"Error processing file for nu = {nu}: {str(e)}")

df = pd.DataFrame(all_data)

df.to_csv('results.csv', index=False)


import numpy as np
import pandas as pd

nu_values = np.arange(0.01, 2.50, 0.01)
all_data = []

for nu in nu_values:
    filename = f"pred_tables//5000_5000/range_2//sigmae_0.32//statespace//partial_results_nu_{nu:.2f}_statespace_sigmae_0.32.npz"
    
    try:
        data = np.load(filename)
        
        mu_errors = data['errors_mu']
        sigma_errors = data['errors_sigma']
        
        for m in range(6):
            row = {
                'nu': nu,
                'method': 'fourier',
                'mu_error': mu_errors[m],  
                'sigma_error': sigma_errors[m],  
                'm': m+1
            }
            all_data.append(row)
            
    except Exception as e:
        print(f"Error processing file for nu = {nu}: {str(e)}")

df = pd.DataFrame(all_data)

df.to_csv('results.csv', index=False)
