from itertools import product

import csv
import glob
import sys

def count_8mer(directory_path, output_name):
    nucleotides = ['A', 'C', 'G', 'T']
    all_8mers = [''.join(mer) for mer in product(nucleotides, repeat=8)]

    summed_counts = {mer: 0 for mer in all_8mers}

    txt_files = glob.glob(f'{directory_path}/*.bed')

    for txt_file in txt_files:
        counts = {mer: 0 for mer in all_8mers}

        with open(txt_file, 'r') as file:
            for line in file:
                parts = line.strip().split('\t')
                if len(parts) >= 2:  
                    eight_mer = parts[3]  
                    if eight_mer in counts:  
                        counts[eight_mer] += 1

        for eight_mer in all_8mers:
            summed_counts[eight_mer] += counts[eight_mer]

    output_path = f'{directory_path}/{output_name}'

    sorted_8mers = sorted(summed_counts.items(), key=lambda item: item[1], reverse=True)

    with open(output_path, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        for eight_mer, count in sorted_8mers:
            csvwriter.writerow([eight_mer, count])

    print(f"Summed 8-mer counts have been saved to {output_path}, sorted in decreasing order of counts.")

if __name__ == "__main__":
    directory_path=sys.argv[1]
    output_name=sys.argv[2]
    count_8mer(directory_path, output_name)
