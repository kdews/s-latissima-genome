#!/bin/bash
#SBATCH -J 6_hal2fasta
#SBATCH -o %x.log
#SBATCH -t 1-0
#SBATCH --mem=22G
#SBATCH -c 16

date

hal2fasta cactus-steps-output/Anc3.hal Anc3 --hdf5InMemory > cactus-steps-output/Anc3.fa
