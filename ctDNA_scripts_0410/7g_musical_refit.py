import numpy as np
import pandas as pd
import pickle
import musical
import os
import re
import sys
from musical.utils import get_sig_indices_associated
from musical.refit import refit

def refit_matrices(target_path, model_path, thresh_refit, save_path):
    sigs = [f for f in os.listdir(model_path) if f.endswith('.csv')]

    target_files = os.listdir(target_path)
    if 'length_matrix.csv' in target_files:
        matrices = ['length_matrix.csv', 'out5p_matrix.csv', 'in5p_matrix.csv']
    else:
        matrices = ['out5p_norm_matrix.csv', 'in5p_norm_matrix.csv']
   
    for matrix in matrices:
        df_x = pd.read_csv(os.path.join(target_path, matrix), sep=',', index_col=0)
        print(os.path.join(target_path, matrix))
        
        #col_sums = df_x.sum(axis=0)
        #selected_columns = col_sumsi[col_sums > 100].index
        #df_x = df_x[selected_columns]#
 
        if (df_x.sum(axis=0) >= 500).all():
            print(f"All columns in {matrix} have sums >= 500.")
        else:
            print(f"Not all columns in {matrix} meet the sum >= 500 criteria.")

        pattern = ""
        if "length" in matrix:
            pattern = "length"
        elif "out5p" in matrix:
            pattern = "out5p"
        elif "in5p" in matrix:
            pattern = "in5p"

        sig_file = next((sig for sig in sigs if re.search(pattern, sig)), None)
        if sig_file:
            df_sigs = pd.read_csv(os.path.join(model_path, sig_file), sep=',', index_col=0)
            print(os.path.join(model_path, sig_file))
            W_catalog = df_sigs
            df_x.index = W_catalog.index

            #thresh_refit = 0.01
            #thresh_refit = 0

            H_s, model_H = refit(X=df_x, W=W_catalog, thresh=thresh_refit, connected_sigs=True)

            output_file_name = f'refit_{pattern}.csv'
            H_s.to_csv(os.path.join(save_path, output_file_name))
            print(output_file_name)
        else:
            print(f"No signature file found for pattern: {pattern}")

if __name__ == "__main__":
    target_path = sys.argv[1] # directory where matrix files are 
    model_path = '/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/shared_tools/fixed_sig'
    #model_path = sys.argv[2] # fixed signature 
    thresh_refit = float(sys.argv[2]) # 0 or 0.01 try both
    save_path = sys.argv[3] # directory where things will be saved 
    refit_matrices(target_path, model_path, thresh_refit, save_path)
