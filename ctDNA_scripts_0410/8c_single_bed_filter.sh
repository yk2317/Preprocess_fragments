#!/bin/bash

SCRIPT_DIR=$1
OUTPUT_DIR=$2
MANIFEST_FILE="/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/shared_tools/pan_v2_Manifest.wo_iSNP.wo_decoy.merged.hg38.bed"

for input_file in $SCRIPT_DIR/*.txt; do
    filename=$(basename "$input_file" .txt)
    output_file="${OUTPUT_DIR}/filtered_${filename}.bed"

    if [ ! -f "$output_file" ]; then
        bedtools intersect -a "$input_file" -b "$MANIFEST_FILE" -wa -u > "$output_file"
    fi
done
