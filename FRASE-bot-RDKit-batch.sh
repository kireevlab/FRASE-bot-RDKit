#!/bin/bash

#SBATCH -t 10:00
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=1GB 
#SBATCH --job-name FRASESCREEN
#SBATCH --array 0-10598

module load miniconda3
eval "$(conda shell.bash hook)"
conda activate rdkit-tools

start_time=$(date +%s)
folder_id=0
arr=(`ls ../FRASE_database_PDB${folder_id}/`)
inp=${arr[$SLURM_ARRAY_TASK_ID]}

python FRASE-bot-RDKit.py $inp $folder_id CIB1-Ana-prot.pdb

end_time=$(date +%s)
elapsed_sec=$((end_time - start_time))

echo "Elapsed time: $elapsed_sec s"
