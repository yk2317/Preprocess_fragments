#!/bin/bash

TXT_DIR="$1"
OUTPUT_DIR="$2"
number="$3"

for file in "$TXT_DIR"/*.txt; do
base_name=$(basename "$file")
output_file="$OUTPUT_DIR/subsampled_$base_name"

sbatch -c 1 -p park -A park_contrib -t 30:00 --mem=300 --wrap="bash ./subsample_0213.sh "$file" "$output_file" "$number""
done
