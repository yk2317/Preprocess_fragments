#!/bin/bash

BAM_LIST=$1
OUT_DIR=$2 

while read -r file_name; do
    echo "Processing BAM file: ${file_name}"
    base_name=$(basename -- "$file_name")
    sbatch --job-name="E_${base_name}" -c 1 -p priopark -A park_contrib -t 12:00:00 --mem=300 --mail-type=FAIL --mail-user=yk2317@gmail.com --wrap="python /home/yok929/ctDNA_scripts/extraction_allcols_umi.py ${file_name} ${OUT_DIR}"
done < "${BAM_LIST}"
