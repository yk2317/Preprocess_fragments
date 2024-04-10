#!/bin/bash

input_file=$1
bed_file=$2
output_dir=$3

filename=$(basename "$input_file" .txt)
output_file="${output_dir}/filtered_${filename}.bed"

bedtools intersect -a "$input_file" -b "$bed_file" -wa -u | awk '{print $1"\t"$2"\t"$3"\t"$6}' > "$output_file"

echo "finished $input_file"
