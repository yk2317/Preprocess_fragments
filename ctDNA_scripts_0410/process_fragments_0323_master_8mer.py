import sys
from itertools import product
import numpy as np
import pickle

def generate_all_8mers():
    bases = ['A', 'C', 'G', 'T']
    return [''.join(p) for p in product(bases, repeat=8)]

def save_count_from_master(master_file, output_dir):
    all_8mers = generate_all_8mers()
    motif_counts = {mer: 0 for mer in all_8mers}
    bins = np.array([50, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 700])
    length_counts = np.zeros(len(bins) - 1)

    base_name = master_file.split('/')[-1].replace('_pymaster.txt', '')

    with open(master_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            motif = parts[3]
            length = int(parts[2]) - int(parts[1]) + 1

            motif_counts[motif] = motif_counts.get(motif, 0) + 1

            bin_index = np.digitize(length, bins, right=False) - 1
            bin_index = min(bin_index, len(length_counts) - 1)
            length_counts[bin_index] += 1

    length_counts_dict = {f"{bins[i]}-{bins[i+1]}": length_counts[i] for i in range(len(bins)-1)}

    all_counts = {
        'length': length_counts_dict,
        'motif_counts': motif_counts
    }

    output_pkl_file = f"{output_dir}/{base_name}_counts.pkl"
    with open(output_pkl_file, 'wb') as pkl_file:
        pickle.dump(all_counts, pkl_file)

    print(f"All counts have been saved to {output_pkl_file}")

if __name__ == "__main__":
    master_file = sys.argv[1]
    output_dir = sys.argv[2]

    save_count_from_master(master_file, output_dir)
