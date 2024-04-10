sbatch -c 1 -p park -A park_contrib -t 3:00:00 --mem=20G --wrap="python count_8mer.py ./pos count_8mer_MSK_pos.csv"
sbatch -c 1 -p park -A park_contrib -t 3:00:00 --mem=20G --wrap="python count_8mer.py ./neg count_8mer_MSK_neg.csv"
