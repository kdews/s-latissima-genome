#!/bin/bash
#SBATCH -J 4_cactus-blast
#SBATCH -o cactus-logs/%x.log
#SBATCH -t 1-0
#SBATCH --mem=72G
#SBATCH -c 16

date

cactus-blast cactus_prepare_steps_jobstore/3 cactus-steps-output/s_lat_alignment.txt cactus-steps-output/Anc3.paf --root Anc3  --maxCores 14 --maxMemory 60.5Gi --logFile cactus-steps-output/logs/blast-Anc3.log
