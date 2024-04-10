#!/bin/bash

base_path=$1 #run this in base folder: ex) MSK_hg38/duplex/
analysis=$2
method=$3
min_n=$4
max_n=$5

declare -a types=("length" "in5p" "out5p")

if [ ! -f "${base_path}/matrix/${analysis}/length_matrix.csv" ]; then
    types=("in5p_norm" "out5p_norm")
fi

if [[ "$method" == "NMF" ]]; then
    mem="3G"
    partition="park"
    account="park_contrib"
    time="12:00:00"
elif [[ "$method" == "MVNMF" ]]; then
    mem="5G"
    partition="park"
    account="park_contrib"
    time="5-00:00"
else
    echo "Invalid method specified. Exiting."
    exit 1
fi

for type in "${types[@]}"; do
    input="${base_path}/matrix/${analysis}/${type}_matrix.csv"
    output="${base_path}/model/${analysis}/${type}_${method}_${min_n}_${max_n}_model.pkl"

    sbatch --job-name="mt" -c 1 -p $partition -A $account -t $time --mem=$mem --wrap="python /home/yok929/ctDNA_scripts/train_0201.py $input $output $method $min_n $max_n"
    echo "Submitted: $base_path analysis for $type with $method, min_n=$min_n, max_n=$max_n"
done
