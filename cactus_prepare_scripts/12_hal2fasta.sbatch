#!/bin/bash
#SBATCH -J 12_hal2fasta
#SBATCH -o cactus-logs/%x.log
#SBATCH -t 1-0
#SBATCH --mem=22G
#SBATCH -c 16

date

hal2fasta cactus-steps-output/s_lat_alignment.hal Anc0 --hdf5InMemory > cactus-steps-output/Anc0.fa
