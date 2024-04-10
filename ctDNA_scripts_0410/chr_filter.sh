#!/bin/bash

file="$1"
outdir="$2"

# Ensure the output directory exists
mkdir -p "$outdir"

# Extract the base name without the .bed extension
base_name=$(basename "$file" .bed)

# Sort and merge the filtered file, saving the output with a new name in the specified output directory
bedtools sort -i "$file" | bedtools merge > "${outdir}/merged_${base_name}.bed"
