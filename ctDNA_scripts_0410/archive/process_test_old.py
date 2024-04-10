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

def process_8mer(motif, is_positive_strand):
    if is_positive_strand:
        out5p = flip(motif[:4])
        in5p = motif[4:]
    else:
        out5p = base_flip(motif[4:])
        in5p = flip(base_flip(motif[:4]))
    return out5p, in5p

def save_master(input_file, output_dir, genome):
    base_name = input_file.split('/')[-1].replace('_pass.txt', '')
    output_file = f"{output_dir}/{base_name}_pymaster.txt"
    print(f"Starting: {base_name}")

    is_positive_strand = "pos" in base_name

    processed_lines = []
    with open(input_file, 'r') as f:
        for line in f:
            chr, start, end, name, gc_content = line.strip().split('\t')
            start, end = int(start), int(end)
            start_adjusted = start if is_positive_strand else end

            start += 1  # Adjust for 0-based indexing if necessary

            # Fetch motif based on adjusted start position
            motif = genome[chr][start_adjusted-4:start_adjusted+4] if chr in genome and start_adjusted-4 >= 0 else "NNNNNNNN"

            if "N" not in motif:
                out5p, in5p = process_8mer(motif, is_positive_strand)
                length = abs(end - start) + 1
                processed_line = [chr, str(start), str(end), name, gc_content, motif, out5p, in5p, str(length)]
                processed_lines.append('\t'.join(processed_line))

    with open(output_file, 'w') as f:
        for line in processed_lines:
            f.write(line + '\n')

    print("Master file saved")
    return processed_lines, base_name

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    genome_pickle = "/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/shared_tools/genome.hg38.pickle"

    with open(genome_pickle, 'rb') as handle:
        genome_data = pickle.load(handle)

    processed_lines, base_name = save_master(input_file, output_dir, genome_data)
