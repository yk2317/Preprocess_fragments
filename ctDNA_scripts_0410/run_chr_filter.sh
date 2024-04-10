#!/bin/bash

SCRIPT_DIR=$1
outdir=$2

for file in "$SCRIPT_DIR"/*.bed; do
    base_name=$(basename "$file" .bed)

    sbatch --job-name="chr_${base_name}" -c 1 -p short -t 12:00:00 --mem=500 --mail-type=FAIL --mail-user=yk2317@gmail.com --wrap="bash ./chr_filter.sh $file $outdir"
done
