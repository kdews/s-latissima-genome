#!/bin/bash
#SBATCH -J last-dotplot
#SBATCH -p shared
#SBATCH --mem=5gb
#SBATCH --time=01:00:00
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
Plots genome alignment using last-dotplot.

Usage: sbatch ${SLURM_JOB_NAME}.sbatch <alignment_file>
Requires:
 - Last (https://github.com/lpryszcz/last/)"
  exit 0
fi

# Load Last conda env
conda activate last
last-dotplot -h
# Input data
aln="$1"
aln_no_ext="$(basename ${aln%.*})"
out="${aln_no_ext}.png"

# Create last-dotplot image
cmd="last-dotplot $aln $out"
echo "$cmd"
$cmd
