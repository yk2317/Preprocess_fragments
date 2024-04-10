#!/bin/bash

file=$1
output_file=$2
num=$3

total_lines=$(wc -l < "$file")

lines_to_sample=$((total_lines / num))

shuf -n "$lines_to_sample" "$file" > "$output_file"
