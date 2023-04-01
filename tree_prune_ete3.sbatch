#!/bin/bash
#SBATCH -J tree_prune_ete3
#SBATCH -p main
#SBATCH --mem=5gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=12
#SBATCH -o %x.log

# Print date and time
date

# Source conda and activate ete3 env
cond=~/.conda_for_sbatch.sh
[[ -a "$cond" ]] && source "$cond" || \
{ echo "Error on source of $cond"; exit 1; }

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Runs ETE3-based Python script to prune a given Newick tree to contain only
desired species.
Usage: sbatch ${SLURM_JOB_NAME}.sbatch
Requires: 
 - Python v2.7.15 (https://www.python.org/download/releases/2.7/)
 - ETE Toolit v3.1.2 (http://etetoolkit.org/)
 - ETE dependencies..."
  exit 0
fi

# Define variables
script_dir="s-latissima-genome"
[[ -d $script_dir ]] && script_dir="${script_dir}/"
py_script=${SLURM_JOB_NAME}.py
# Cactus-formatted seqFile, i.e.,
# NEWICK tree
# name1 path1
# name2 path2
# ...
# nameN pathN
seq_file="s_lat_alignment.txt"


# Run tree pruning Python script
cmd="python ${script_dir}${py_script} ${script_dir}${seq_file}"
echo "$cmd"
$cmd