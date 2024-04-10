import sys
import os
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import pandas as pd
import time
import scipy as sp
import pickle
import musical

def process_and_save_model_data(working_dir, subdirectory, file_name, output_dir, dataset_prefix):
    base_name = file_name.replace("_model.pkl", "")
    
    file_path = os.path.join(working_dir, subdirectory, file_name)
    with open(file_path, 'rb') as f:
        model = pickle.load(f)
    
    nums = model.n_components

    H_df = model.H_df
    H_df.to_csv(os.path.join(output_dir, f'{dataset_prefix}_H_{base_name}_{nums}.csv'))
    
    W_df = model.W_df
    W_df.to_csv(os.path.join(output_dir, f'{dataset_prefix}_W_{base_name}_{nums}.csv'))
    
    model.plot_selection()
    plt.savefig(os.path.join(output_dir, f'{dataset_prefix}_{base_name}_{nums}_selection.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    fig = musical.sigplot_bar(model.W, sig_type="")
    plt.savefig(os.path.join(output_dir, f'{dataset_prefix}_{base_name}_{nums}_sig.png'), dpi=300, bbox_inches='tight')
    plt.close()

def process_files(working_dir, subdirectory, dataset_prefix):
    subdirectory_path = os.path.join(working_dir, subdirectory)
    output_dir = os.path.join(working_dir, "output", dataset_prefix)
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    pkl_files = [f for f in os.listdir(subdirectory_path) if f.endswith('_model.pkl')]
    
    for file_name in pkl_files:
        process_and_save_model_data(working_dir, subdirectory, file_name, output_dir, dataset_prefix)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python script.py <working_dir> <subdirectory> <dataset_prefix>")
        sys.exit(1)
    
    working_dir = sys.argv[1]
    subdirectory = sys.argv[2]
    dataset_prefix = sys.argv[3]  
    process_files(working_dir, subdirectory, dataset_prefix)
