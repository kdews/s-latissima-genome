#!/bin/bash
#SBATCH -J 0_cactus-preprocess
#SBATCH -o cactus-logs/%x.log
#SBATCH -t 1-0
#SBATCH --mem=22G
#SBATCH -c 16

date

cactus-preprocess cactus_prepare_steps_jobstore/0 s-latissima-genome/s_lat_alignment_new_filt.txt cactus-steps-output/s_lat_alignment_new_filt.txt --inputNames Saccharina_japonica Saccharina_latissima Macrocystis_pyrifera Undaria_pinnatifida Ectocarpus_siliculosus  --maxCores 14 --maxMemory 18.6Gi --logFile cactus-steps-output/logs/preprocess-Saccharina_japonica.log
