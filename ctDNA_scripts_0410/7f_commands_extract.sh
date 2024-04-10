#!/bin/bash
# run this in folder where output should be created 
path=$1

cd "$path"

python /home/yok929/ctDNA_scripts/7b_musical_extract.py . allreads MSK_allreads
python /home/yok929/ctDNA_scripts/7b_musical_extract.py . allreads_norm MSK_allreads_norm
python /home/yok929/ctDNA_scripts/7b_musical_extract.py . variants MSK_variants
python /home/yok929/ctDNA_scripts/7b_musical_extract.py . variants_norm MSK_variants_norm

python /home/yok929/ctDNA_scripts/7b_musical_extract.py . shift_allreads MSK_shift_allreads
python /home/yok929/ctDNA_scripts/7b_musical_extract.py . shift_allreads_norm MSK_shift_allreads_norm

echo "musical extract in $path"
