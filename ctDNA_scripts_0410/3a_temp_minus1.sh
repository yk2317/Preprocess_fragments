#!/bin/bash

SCRIPT_DIR=$1
OUTPUT_DIR=$2

for file in $SCRIPT_DIR/*.bed
do
    base_name=$(basename $file .bed)
    #output_file="$OUTPUT_DIR/${base_name}_master.txt"

    #if [ ! -f "$output_file" ]; then
    sbatch --job-name="F_${base_name}" -c 1 -p priopark -A park_contrib -t 8:00:00 --mem=20G --mail-type=FAIL --mail-user=yk2317@gmail.com --wrap="python /home/yok929/ctDNA_scripts/process_fragments_0210_panel_minus1.py $file $OUTPUT_DIR"
    #fi
done
