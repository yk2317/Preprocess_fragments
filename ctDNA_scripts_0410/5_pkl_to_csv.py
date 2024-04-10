import os
import sys
import pickle
import pandas as pd

#input_dir = '/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/extraction/output_gerberg'  
#output_dir = '/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/extraction/output_gerberg/csv' 

input_dir = sys.argv[1]
output_dir = sys.argv[2]

os.makedirs(output_dir, exist_ok=True)

pkl_files = [file for file in os.listdir(input_dir) if file.endswith('.pkl')]

for pkl_file in pkl_files:
    pkl_path = os.path.join(input_dir, pkl_file)
    
    with open(pkl_path, 'rb') as file:
        data = pickle.load(file)
    
    for key in ['length', 'in5p', 'out5p']:
        if key in data:
            df = pd.DataFrame(list(data[key].items()), columns=['Key', 'Count'])
            df['Count'] = df['Count'].astype(int)
            
            base_name = os.path.splitext(pkl_file)[0]  
            csv_filename = f"{base_name}_{key}.csv"
            csv_path = os.path.join(output_dir, csv_filename)
            
            df.to_csv(csv_path, index=False, header=False)
            print(f"{csv_path} has been saved.")
