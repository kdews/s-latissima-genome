#!/bin/bash
#SBATCH -J 11_cactus-align
#SBATCH -o cactus-logs/%x.log
#SBATCH -t 1-0
#SBATCH --mem=333G
#SBATCH -c 16

date

cactus-align cactus_prepare_steps_jobstore/8 cactus-steps-output/s_lat_alignment.txt cactus-steps-output/Anc0.paf cactus-steps-output/s_lat_alignment.hal --root Anc0  --maxCores 14 --maxMemory 279.4Gi --logFile cactus-steps-output/logs/align-Anc0.log
