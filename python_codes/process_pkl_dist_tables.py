import os
import h5py
import numpy as np
import pickle

input_dir = "distance_tables/raw_tables/"
output_dir = "distance_tables/raw_tables/"

os.makedirs(output_dir, exist_ok=True)
for filename in os.listdir(input_dir):
    if filename.endswith(".pkl"):
        with open(os.path.join(input_dir, filename), 'rb') as file:
            dict_table = pickle.load(file)
        if isinstance(dict_table, dict):
            h5_filename = os.path.join(output_dir, f"{filename.replace('.pkl', '.h5')}")
            
            with h5py.File(h5_filename, 'w') as h5file:
                for key, array in dict_table.items():
                    if isinstance(array, np.ndarray):
                        h5file.create_dataset(key, data=array)      
            print(f"Converted {filename} to {h5_filename}")
print("Conversion complete.")

