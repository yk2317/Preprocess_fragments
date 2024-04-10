#!/bin/bash

SCRIPT_DIR=$1
OUTPUT_DIR=$2

for file in $SCRIPT_DIR/*.txt
do
    input_file=$file
    base_name=$(basename "$input_file" raw_variant.txt)

    pos_file="${OUTPUT_DIR}/${base_name}_variant_pos.txt"
    neg_file="${OUTPUT_DIR}/${base_name}_variant_neg.txt"

    awk -v pos=$pos_file -v neg=$neg_file '{
        value = $6 >= 0 ? $6 : -$6;
        if (value >= 70 && value <= 700) {
            if ($5 == "+")
                print $0 >> pos;
            else if ($5 == "-")
                print $0 >> neg;
        }
    }' "$input_file"

    echo "Data separated into $pos_file and $neg_file"
done
