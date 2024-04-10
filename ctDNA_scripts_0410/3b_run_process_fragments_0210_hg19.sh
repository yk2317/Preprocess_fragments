#!/bin/bash

#SCRIPT_DIR=$1
#OUTPUT_DIR=$2
# has to have chr prefix 
SCRIPT_DIR="/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/extraction/yonsei_stitch/filtered/filter2"
OUTPUT_DIR="/n/data1/hms/dbmi/park/ctDNA_loci_project/fragmentomics/extraction/output_yonsei_stitch"

for file in $SCRIPT_DIR/*.bed
do
    base_name=$(basename $file .bed)
    #output_file="$OUTPUT_DIR/${base_name}_master.txt"

    #if [ ! -f "$output_file" ]; then
    sbatch --job-name="YF_${base_name}" -c 1 -p short -t 4:00:00 --mem=100G --mail-type=FAIL --mail-user=yk2317@gmail.com --wrap="python ./process_fragments_0212_hg19_panel.py $file $OUTPUT_DIR"
    #fi
done
