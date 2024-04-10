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

def process_8mer(motif, base_name, counts):
    if "pos" in base_name:
        out5p = flip(motif[:4])
        in5p = motif[4:]
    else:
        out5p = base_flip(motif[4:])
        in5p = base_flip(flip(motif[:4]))
    counts['out5p'][out5p] = counts['out5p'].get(out5p, 0) + 1
    counts['in5p'][in5p] = counts['in5p'].get(in5p, 0) + 1

def save_master(input_file, output_dir, genome):
    base_name = input_file.split('/')[-1].replace('_pass.bed', '')
    output_file = f"{output_dir}/{base_name}_pymaster.txt"
    print(f"Starting: {base_name}")

    processed_lines = []
    with open(input_file, 'r') as f:
        for line in f:
            chr, start, end, name, gc_content = line.strip().split('\t')
            start, end = int(start), int(end)
            #start_adjusted = start if "pos" in name else end

            if "pos" in name:
                start_adjusted = start -1
                end += 1
            else:
                start_adjusted = end -1
                start += 1
            motif = genome[chr][start_adjusted-4:start_adjusted+4] if chr in genome and start_adjusted-4 >= 0 else "NNNNNNNN"

            if "N" not in motif:
                length = abs(end - start) + 1
                processed_line = [chr, str(start), str(end), name, gc_content, motif]
                processed_lines.append('\t'.join(processed_line))

    with open(output_file, 'w') as f:
        for line in processed_lines:
            f.write(line + '\n')

    print("Master file saved")
    return processed_lines, base_name

def save_count(processed_lines, output_dir, base_name):
    all_4mers = generate_all_4mers()
    motif_counts = {'out5p': {mer: 0 for mer in all_4mers}, 'in5p': {mer: 0 for mer in all_4mers}}
    bins = np.array([0, 50, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 700])
    length_counts = np.zeros(len(bins) - 1)

    for line in processed_lines:
        parts = line.split('\t')
        motif = parts[5] 
        length = int(parts[2]) - int(parts[1])  

        process_8mer(motif, 'some_base_name', motif_counts)

        bin_index = np.digitize(length, bins, right=False) - 1
        bin_index = min(bin_index, len(length_counts) -1)
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
    input_file=sys.argv[1]
    output_dir=sys.argv[2]
    genome_pickle = "/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/shared_tools/genome.hg38.pickle"

    with open(genome_pickle, 'rb') as handle:
        genome_data = pickle.load(handle)

    processed_lines, base_name = save_master(input_file, output_dir, genome_data)
    save_count(processed_lines, output_dir, base_name)
