#!/bin/bash

BAM_LIST=$1
OUT_DIR=$2 

# gerberg 12hr, MSK 2 days, Yonsei 3 days, Wyatt 5-7 days

while read -r file_name; do
    echo "Processing BAM file: ${file_name}"
    base_name=$(basename -- "$file_name")
    sbatch --job-name="E_${base_name}" -c 1 -p medium -t 3-00:00 --mem=300 --mail-type=FAIL --mail-user=yk2317@gmail.com --wrap="python /home/yok929/ctDNA_scripts/extraction_0209_noconsensus.py ${file_name} ${OUT_DIR}"
done < "${BAM_LIST}"
