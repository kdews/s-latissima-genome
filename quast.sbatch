#!/bin/bash
#SBATCH -J quast
#SBATCH -p oneweek
#SBATCH --mem=15gb
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=12
#SBATCH -o %x.log

# Print date and time
date

# Load conda and BUSCO env
cond=~/.conda_for_sbatch.sh
if [[ -a "$cond" ]]
then
  source "$cond"
else
  echo "Error on source of $cond"
  exit 1
fi
conda activate quast
# Load GNU Parallel
module load parallel

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Runs QUAST on a given assembly 
Usage: sbatch quast.sbatch <assembly.fasta>
Requires: QUAST (https://quast.sourceforge.net)"
  exit 0
fi


# Define variables
assembly="$1"
assembly_no_ext="$(basename "${assembly%.*}")"
# Shorten SLURM output variables
job="$SLURM_JOB_NAME"
nthrd="$SLURM_CPUS_PER_TASK"
# Output directory name
out="${job}_${assembly_no_ext}"

# quast="quast --large -k -b --report-all-metrics --threads $nthrd -o $out"
quast="quast --large -k -b --threads $nthrd -o $out"
cmd="$quast $assembly"
echo "$cmd"
$cmd
