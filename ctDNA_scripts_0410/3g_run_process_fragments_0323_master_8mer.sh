#!/bin/bash

SCRIPT_DIR=$1
OUTPUT_DIR=$2
# MSK 2hr, Wyatt 5-8 hrs, memory 1G is enough regardless of file size
for file in $SCRIPT_DIR/*.txt
do
    base_name=$(basename $file .txt)

    sbatch --job-name="8M_${base_name}" -c 1 -p park -A park_contrib -t 8:00:00 --mem=1G --mail-type=FAIL --mail-user=yk2317@gmail.com --wrap="python /home/yok929/ctDNA_scripts/process_fragments_0323_master_8mer.py $file $OUTPUT_DIR"
done
