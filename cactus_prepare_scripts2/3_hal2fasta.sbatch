#!/bin/bash
#SBATCH -J 3_hal2fasta
#SBATCH -o cactus-logs/%x.log
#SBATCH -t 1-0
#SBATCH --mem=22G
#SBATCH -c 16

date

hal2fasta cactus-steps-output/Anc1.hal Anc1 --hdf5InMemory > cactus-steps-output/Anc1.fa
