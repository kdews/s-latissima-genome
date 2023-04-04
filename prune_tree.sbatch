#!/bin/bash
#SBATCH -J prune_tree
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
Takes Newick tree and list of species, and uses R's ape to prune it to contain
only desired species. 
Usage: sbatch ${SLURM_JOB_NAME}.sbatch <newick_tree> <species_list.txt>
Requires: 
 - R v4.2.2
 - ape v5.7-1 (https://cran.r-project.org/package=ape)"
  exit 0
fi

# Define variables
script_dir="s-latissima-genome"
[[ -d $script_dir ]] && script_dir="${script_dir}/"
r_script=${SLURM_JOB_NAME}.R
newick_tree=$1
species_list=$2




# Run tree pruning Python script
cmd="Rscript ${script_dir}${r_script} ${newick_tree} ${species_list}"
echo "$cmd"
$cmd