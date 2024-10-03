import os
import pandas as pd
import tensorflow as tf
import pickle

folder_path = 'pred_tables'
for filename in os.listdir(folder_path):
    if filename.endswith('.pkl'):
        file_path = os.path.join(folder_path, filename)
        new_filename = filename.replace('.pkl', '.csv')
        new_file_path = os.path.join(folder_path, new_filename)
        if not os.path.exists(new_file_path):
            with open(file_path, 'rb') as file:
                data = pickle.load(file)
            if 'errors' in data.columns:
                data['nu'] = data['nu'].apply(lambda x: x.numpy() if isinstance(x, tf.Tensor) else x)
                errors_df = pd.DataFrame(data['errors'].tolist(), columns=[f'error_{i+1}' for i in range(len(data['errors'][0]))])
                result_df = pd.concat([data[['nu']], errors_df], axis=1)
                result_df.to_csv(new_file_path, index=False)
                print(f"Processed and saved: {new_file_path}")
            else:
                print(f"'errors' column not found in {filename}. Skipping.")

print("Processing complete.")

