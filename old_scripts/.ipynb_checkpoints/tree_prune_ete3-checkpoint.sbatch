#!/bin/bash
#SBATCH -J tree_prune_ete3
#SBATCH -p oneweek
#SBATCH --constraint=xeon-2640v3
#SBATCH --mem=50gb
#SBATCH --time=7-0
#SBATCH --cpus-per-task=16
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
# Newick tree file
tree=assemblies/brown_algae_newick_tree.txt
# List of desired species in tree
spec_list=assembles/pretty_names.txt

# Run tree pruning Python script
cmd="python ${script_dir}${py_script} $tree $spec_list"
echo "$cmd"
$cmd