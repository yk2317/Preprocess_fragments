#!/bin/bash

vcf_file=$1
bam_file_path=$2
out_dir=$3

while read -r file_name; do
    echo "Processing BAM file: ${file_name}"
    base_name=$(basename -- "$file_name")
    sbatch --job-name="V_${base_name}" -c 1 -p priopark -A park_contrib -t 30:00 --mem=300 --mail-type=FAIL --mail-user=yk2317@gmail.com --wrap="python /home/yok929/ctDNA_scripts/find_0218.py ${vcf_file} ${bam_file_path} ${out_dir}"
done < "${bam_file_path}"
