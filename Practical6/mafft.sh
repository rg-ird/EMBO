#!/bin/bash
#
#SBATCH --job-name=concensus
#SBATCH --ntasks=12
#SBATCH --mem-per-cpu=10000
#SBATCH --cpus-per-task=1
#SBATCH --partition=long

module load mafft/7.407
module load fasttree/2.1.10
#module load muscle/3.8.1551

FastTree ALL.fas.maaft > ALL.fas.maaft.tree
#muscle -in ALL.cln -out ALL.cln.muscle
#mafft --thread 12 ALL.fas > ALL.fas.maaft  
