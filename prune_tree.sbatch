#!/bin/bash
#SBATCH -J prune_tree
#SBATCH -p main
#SBATCH --mem=5gb
#SBATCH --time=01:00:00
#SBATCH -o %x.log

# Print date and time
date

# Load R module
module purge
module load gcc/11.3.0 openblas/0.3.20 r/4.2.2 

# Kelp Newick tree from https://doi.org/10.1016/j.ympev.2019.04.012
tree="https://ars.els-cdn.com/content/image/1-s2.0-S1055790319300892-mmc1.txt"

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]] || (( $# < 1 ))
then
  echo "\
Takes Newick tree and list of species, and uses R's ape to prune it to contain
only desired species. 
Usage: sbatch ${SLURM_JOB_NAME}.sbatch <species> [tree]

Requires: 
 - R v4.2.2
 - ape v5.7-1 (https://cran.r-project.org/package=ape)
 
Default Newick tree: Kelp phylogeny from Starko et al. 2019 ($tree)"
  exit 0
fi

# Define variables
script_dir="s-latissima-genome"
[[ -d $script_dir ]] && script_dir="${script_dir}/"
r_script=${SLURM_JOB_NAME}.R
species="$1"
(( $# > 1 )) && tree="$2"
if [[ ! -f $tree ]] && [[ $tree == http?(s)://*.txt ]]
then
    if [[ -f "$(basename $tree)" ]]
    then
      tree="$(basename $tree)"
    else
      wget "$tree"
      tree="$(basename $tree)"
    fi
fi
echo "Species list: $species."
echo "Input Newick tree: $tree"
[[ -f $species ]] || { echo "Error: species list not found."; exit 1; }
[[ -f $tree ]] || { echo "Error: tree not found."; exit 1; }

# Run tree pruning R script
cmd="Rscript ${script_dir}${r_script} ${species} ${tree}"
echo "$cmd"
$cmd