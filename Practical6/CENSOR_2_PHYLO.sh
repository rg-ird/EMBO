#!/bin/bash
#
#SBATCH --job-name=censor
#SBATCH --ntasks=12
#SBATCH --mem-per-cpu=5000
#SBATCH --cpus-per-task=1
#SBATCH --partition=long
#
module load emboss/6.6.0
module load mafft/7.515
module load fasttree/2.1.10
#conda activate wise2

time python censor_2_phylolihgt6.py GCA_036785865.1_ASM3678586v1_genomic.fna RTcores-database-wickercode.LineageBianca_REXDB.fa GCA_036785865.1_ASM3678586v1_genomic.fna.map output_folder_HD2 final_HD2_base 0.50 50

