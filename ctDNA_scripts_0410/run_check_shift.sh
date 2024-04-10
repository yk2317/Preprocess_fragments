#!/bin/bash

FILE_PATHS=$1
OUTPUT_DIR=$2

while IFS= read -r file; do
    base_name=$(basename "$file" .txt)

    sbatch --job-name="PLT_${base_name}" -c 1 -p park -A park_contrib -t 3:00:00 --mem=30G --wrap="python /home/yok929/ctDNA_scripts/check_shift.R \"$file\" \"$OUTPUT_DIR\""

done < "$FILE_PATHS"
