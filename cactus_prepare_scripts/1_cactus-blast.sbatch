#!/bin/bash
#SBATCH -J 1_cactus-blast
#SBATCH -o cactus-logs/%x.log
#SBATCH -t 1-0
#SBATCH --mem=72G
#SBATCH -c 16

date

cactus-blast cactus_prepare_steps_jobstore/1 cactus-steps-output/s_lat_alignment_new_filt.txt cactus-steps-output/Anc1.paf --root Anc1  --maxCores 14 --maxMemory 60.5Gi --logFile cactus-steps-output/logs/blast-Anc1.log
