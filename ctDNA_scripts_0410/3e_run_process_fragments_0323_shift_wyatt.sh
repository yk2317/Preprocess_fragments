#!/bin/bash

FILE_PATHS=$1
OUTPUT_DIR=$2
# 99GB 4hr 130GB / 50GB 2hr 70GB 
while IFS= read -r file; do
    base_name=$(basename $file .bed)

    sbatch --job-name="MS_${base_name}" -c 1 -p park -A park_contrib -t 4:00:00 --mem=16G --mail-type=FAIL --mail-user=yk2317@gmail.com --wrap="python /home/yok929/ctDNA_scripts/process_fragments_0323_shift_wyatt.py $file $OUTPUT_DIR"

done < "$FILE_PATHS"
