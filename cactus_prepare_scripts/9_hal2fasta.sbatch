#!/bin/bash
#SBATCH -J 9_hal2fasta
#SBATCH -o %x.log
#SBATCH -t 1-0
#SBATCH --mem=22G
#SBATCH -c 16

date

hal2fasta cactus-steps-output/Anc2.hal Anc2 --hdf5InMemory > cactus-steps-output/Anc2.fa
