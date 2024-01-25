#!/bin/bash
#SBATCH -J 14_halAppendSubtree
#SBATCH -o %x.log
#SBATCH -t 1-0
#SBATCH --mem=22G
#SBATCH -c 16

date

halAppendSubtree cactus-steps-output/s_lat_alignment.hal cactus-steps-output/Anc1.hal Anc1 Anc1 --merge --hdf5InMemory
