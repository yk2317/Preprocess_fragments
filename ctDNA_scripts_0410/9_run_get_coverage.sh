#!/bin/bash

FILE_PATHS=$1
OUTPUT_DIR=$2

while IFS= read -r file; do
    base_name=$(basename "$file" .bam)
    output_file="$OUTPUT_DIR/${base_name}_coverage.bed"

    if [ ! -f "$output_file" ]; then
        echo "running $file"
        sbatch --job-name="C_${base_name}" -c 1 -p park -A park_contrib -t 12:00:00 --mem=500 --mail-type=FAIL --mail-user=yk2317@gmail.com --wrap="bash /home/yok929/ctDNA_scripts/coverage.sh $file $OUTPUT_DIR"
    else
        echo "File $output_file already exists"
    fi
done < "$FILE_PATHS"
