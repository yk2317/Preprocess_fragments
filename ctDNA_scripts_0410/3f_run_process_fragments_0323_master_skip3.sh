#!/bin/bash

SCRIPT_DIR=$1
OUTPUT_DIR=$2
# let's try 50G for input size up to 45G 
for file in $SCRIPT_DIR/*.txt
do
    base_name=$(basename $file .txt)

    sbatch --job-name="skp3_${base_name}" -c 1 -p park -A park_contrib -t 5:00:00 --mem=500 --mail-type=FAIL --mail-user=yk2317@gmail.com --wrap="python /home/yok929/ctDNA_scripts/process_fragments_0323_master_skip3.py $file $OUTPUT_DIR"
done
