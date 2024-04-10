#!/bin/bash

# Capture command line arguments
file=$1  # The input file
num=$2   # The subsampling frequency

# Extract the base name of the file without the path
base_name=$(basename "$file")

# Generate the output file name
output_file="subsampled_${num}_${base_name}"

# Perform the subsampling
awk -v num="$num" 'BEGIN {srand()} !((NR-1) % num) && rand() <= 1/num' "$file" > "$output_file"

echo "Subsampling complete. Output file: $output_file"
