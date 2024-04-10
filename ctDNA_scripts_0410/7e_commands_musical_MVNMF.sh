#!/bin/bash

path=$1

cd "$path"

bash /home/yok929/ctDNA_scripts/7a_run_musical.sh . allreads MVNMF 10 15
bash /home/yok929/ctDNA_scripts/7a_run_musical.sh . allreads_norm MVNMF 10 15
bash /home/yok929/ctDNA_scripts/7a_run_musical.sh . variants MVNMF 10 15
bash /home/yok929/ctDNA_scripts/7a_run_musical.sh . variants_norm MVNMF 10 15

echo "musical run in $path"
