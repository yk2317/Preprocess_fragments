#!/bin/bash

BAM_LIST=$1
OUT_DIR=$2

while read -r file_name; do
    echo "Processing BAM file: ${file_name}"
    base_name=$(basename -- "$file_name")
    output_file="${OUT_DIR}/${base_name%.bam}_pos_pass.txt"

    if [[ -f "$output_file" ]]; then
        echo "Output file already exists: ${output_file}"
        continue
    fi

    echo "not there"
    sbatch --job-name="W_${base_name%.bam}" -c 1 -p park -A park_contrib -t 7-00:00 --mem=500 --mail-type=FAIL --mail-user=yk2317@gmail.com --wrap="python /home/yok929/ctDNA_scripts/extraction_0209_wgs.py ${file_name} ${OUT_DIR}"
done < "${BAM_LIST}"
