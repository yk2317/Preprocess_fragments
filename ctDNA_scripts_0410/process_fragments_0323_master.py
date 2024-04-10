import sys
from itertools import product
import numpy as np
import pickle

def generate_all_4mers():
    bases = ['A', 'C', 'G', 'T']
    return [''.join(p) for p in product(bases, repeat=4)]

def flip(seq):
    return seq[::-1]

def base_flip(seq):
    base_flip_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join([base_flip_map[base] for base in seq])

def process_8mer(motif, is_positive_strand, counts):
    if is_positive_strand:
        out5p = flip(motif[:4])
        in5p = motif[4:]
    else:
        out5p = base_flip(motif[4:])
        in5p = flip(base_flip(motif[:4]))
    counts['out5p'][out5p] = counts['out5p'].get(out5p, 0) + 1
    counts['in5p'][in5p] = counts['in5p'].get(in5p, 0) + 1

def save_count_from_master(master_file, output_dir):
    all_4mers = generate_all_4mers()
    motif_counts = {'out5p': {mer: 0 for mer in all_4mers}, 'in5p': {mer: 0 for mer in all_4mers}}
    bins = np.array([50, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 700])
    length_counts = np.zeros(len(bins) - 1)

    base_name = master_file.split('/')[-1].replace('_pymaster.bed', '')
    is_positive_strand = "pos" in base_name

    with open(master_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            motif = parts[3]
            length = int(parts[2]) - int(parts[1]) + 1

            process_8mer(motif, is_positive_strand, motif_counts)

            bin_index = np.digitize(length, bins, right=False) - 1
            bin_index = min(bin_index, len(length_counts) - 1)
            length_counts[bin_index] += 1

    length_counts_dict = {f"{bins[i]}-{bins[i+1]}": length_counts[i] for i in range(len(bins)-1)}

    all_counts = {
        'length': length_counts_dict,
        'in5p': motif_counts['in5p'],
        'out5p': motif_counts['out5p']
    }

    output_pkl_file = f"{output_dir}/{base_name}_counts.pkl"
    with open(output_pkl_file, 'wb') as pkl_file:
        pickle.dump(all_counts, pkl_file)

    print(f"All counts have been saved to {output_pkl_file}")

if __name__ == "__main__":
    master_file = sys.argv[1]
    output_dir = sys.argv[2]

    save_count_from_master(master_file, output_dir)
