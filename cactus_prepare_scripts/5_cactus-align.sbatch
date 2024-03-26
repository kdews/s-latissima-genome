#!/bin/bash
#SBATCH -J 5_cactus-align
#SBATCH -o cactus-logs/%x.log
#SBATCH -t 1-0
#SBATCH --mem=333G
#SBATCH -c 16

date

cactus-align cactus_prepare_steps_jobstore/4 cactus-steps-output/s_lat_alignment_new_filt.txt cactus-steps-output/Anc3.paf cactus-steps-output/Anc3.hal --root Anc3  --maxCores 14 --maxMemory 279.4Gi --logFile cactus-steps-output/logs/align-Anc3.log
