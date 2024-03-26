#!/bin/bash
#SBATCH -J 13_halAppendSubtree
#SBATCH -o cactus-logs/%x.log
#SBATCH -t 1-0
#SBATCH --mem=22G
#SBATCH -c 16

date

halAppendSubtree cactus-steps-output/s_lat_alignment_new_filt.hal cactus-steps-output/Anc2.hal Anc2 Anc2 --merge --hdf5InMemory
