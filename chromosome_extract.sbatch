#!/bin/bash
#SBATCH -J chromosome_extract
#SBATCH --mem=100mb
#SBATCH --time=01:00:00
#SBATCH -o %x.log

# Print date and time
date

# Load conda
cond=~/.conda_for_sbatch.sh
[[ -a "$cond" ]] && source "$cond" || \
{ echo "Error on source of $cond"; exit 1; }

# Load samtools conda env
conda activate samtools
samtools --version

# Input
species="$1"
[[ -z "$species" ]] && species="species.txt"
fastas=($(awk '{print $2}' "$species"))

# Run samtools index tool on each FASTA
for fasta in "${fastas[@]}"
do
    cmd="samtools faidx $fasta"
    echo "$cmd"
    $cmd
done

