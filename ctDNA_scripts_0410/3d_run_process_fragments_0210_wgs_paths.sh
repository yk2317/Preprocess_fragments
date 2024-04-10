#!/bin/bash

FILE_PATHS=$1  
OUTPUT_DIR=$2

while IFS= read -r file; do
    base_name=$(basename "$file" .txt)
    output_file="$OUTPUT_DIR/${base_name}_pymaster.txt"

    if [ ! -f "$output_file" ]; then
        sbatch --job-name="FB_${base_name}" -c 1 -p park -A park_contrib -t 18:00:00 --mem=230G --mail-type=FAIL --mail-user=yk2317@gmail.com --wrap="python /home/yok929/ctDNA_scripts/process_fragments_0311.py \"$file\" \"$OUTPUT_DIR\""
    else
        echo "File $output_file already exists"
    fi
done < "$FILE_PATHS"
