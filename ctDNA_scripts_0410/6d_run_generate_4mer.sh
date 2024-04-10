#!/bin/bash

path=$1
outdir=$2
direction=$3

for file in $path/*.bed
do
    base_name=$(basename $file .bed)

    sbatch --job-name="10kb_${base_name}" -c 1 -p short -t 9:00:00 --mem=5G --mail-type=FAIL --mail-user=yk2317@gmail.com --wrap="Rscript /home/yok929/ctDNA_scripts/generate_4mer.R $file $outdir $direction"
done
