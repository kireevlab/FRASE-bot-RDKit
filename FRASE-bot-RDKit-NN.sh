#!/bin/bash

#SBATCH --mem 20GB
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --job-name NeuralNetwork

module load miniconda3
conda activate Tensorflow

python FRASE-bot-RDKit-NN.py
