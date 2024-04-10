#!/bin/bash

compiled_vcf_file=$1
bam_list_file=$2
out_directory=$3

while IFS= read -r bam_file; do
    base_name=$(basename "$bam_file" ".bam")

    if [ ! -f "$bam_file" ]; then
        echo "BAM file not found: $bam_file"
        continue
    fi

    sbatch --job-name="V" -c 1 -p park -A park_contrib -t 30:00 --mem=300 --output="${out_directory}/${base_name}_variant.txt" --wrap="python /home/yok929/ctDNA_scripts/find_substitution_0204.py $compiled_vcf_file $bam_file"

done < "$bam_list_file"
