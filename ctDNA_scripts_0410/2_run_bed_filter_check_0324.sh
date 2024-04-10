#!/bin/bash

SCRIPT_DIR=$1
BED_DIR=$2
OUTPUT_DIR=$3

for file in $SCRIPT_DIR/*.txt
do
    base_name=$(basename $file .txt)
    output_file="$OUTPUT_DIR/filtered_${base_name}.bed"

    if [ ! -f "$output_file" ]; then
        sbatch --job-name="BED_${base_name}" -c 1 -p short -t 3:00:00 --mem=300 --mail-type=FAIL --mail-user=yk2317@gmail.com --wrap="bash /home/yok929/ctDNA_scripts/bed_filter_jobs_0209.sh $file $BED_DIR $OUTPUT_DIR"
    else
        echo "output file $output_file already exists"
    fi
done
