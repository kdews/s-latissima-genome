#!/bin/bash
#SBATCH -J lastz
#SBATCH -p htcondor
#SBATCH --mem=0
#SBATCH --time=2-0
#SBATCH -o %x.log

# Print date and time
date

# Load conda
cond=~/.conda_for_sbatch.sh
[[ -a "$cond" ]] && source "$cond" || \
{ echo "Error on source of $cond"; exit 1; }

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Aligns two genomes using LastZ.

Usage: sbatch ${SLURM_JOB_NAME}.sbatch <genome1.fasta> <genome2.fasta>
Requires:
 - LastZ (https://lastz.github.io/lastz)"
  exit 0
fi

# Load LastZ conda env
conda activate lastz
lastz --version
# Input data
genome1="$1"
genome1_no_ext="$(basename ${genome1%.*})"
genome2="$2"
genome2_no_ext="$(basename ${genome2%.*})"
out="${genome1_no_ext}_vs_${genome2_no_ext}.maf"

# Run LastZ alignment of two genomes
cmd="lastz ${genome1}[multiple] ${genome2}[multiple] --format=maf --step=20 --nogapped --output=${out}"
echo "$cmd"
$cmd
