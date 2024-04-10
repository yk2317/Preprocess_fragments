#!/bin/bash
#run this in base folder: ex) MSK_hg38/duplex/
path=$1

cd "$path"

#bash /home/yok929/ctDNA_scripts/7a_run_musical.sh . allreads NMF 10 20
bash /home/yok929/ctDNA_scripts/7a_run_musical.sh . allreads_norm NMF 10 20
#bash /home/yok929/ctDNA_scripts/7a_run_musical.sh . variants NMF 10 20
bash /home/yok929/ctDNA_scripts/7a_run_musical.sh . variants_norm NMF 10 20

#bash /home/yok929/ctDNA_scripts/7a_run_musical.sh . allreads NMF 10 20
#bash /home/yok929/ctDNA_scripts/7a_run_musical.sh . allreads_norm NMF 10 20
#bash /home/yok929/ctDNA_scripts/7a_run_musical.sh . shift_allreads NMF 10 20
bash /home/yok929/ctDNA_scripts/7a_run_musical.sh . shift_allreads_norm NMF 10 20

#bash /home/yok929/ctDNA_scripts/7a_run_musical.sh . combined NMF 10 20
#bash /home/yok929/ctDNA_scripts/7a_run_musical.sh . combined_norm NMF 10 20

echo "musical run in $path"
