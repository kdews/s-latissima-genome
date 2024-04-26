#!/bin/bash
#SBATCH -J prune_tree
#SBATCH -p main
#SBATCH --mem=5gb
#SBATCH --time=01:00:00
#SBATCH -o %x.log

# Initialize script
init="/home1/kdeweese/scripts/s-latissima-genome/script_init.sh"
if [[ -a "$init" ]]; then source "$init"; else { echo "Init err"; exit 1; }; fi

# Load R module
module purge
module load gcc/11.3.0 openblas/0.3.20 r

# Kelp Newick tree from https://doi.org/10.1016/j.ympev.2019.04.012
tree="https://ars.els-cdn.com/content/image/1-s2.0-S1055790319300892-mmc1.txt"

# Help message
if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]] || (( $# < 1 ))
then
  echo "\
Takes Newick tree and assembly file, and uses R's ape to prune it to contain
only desired species.
Usage: sbatch ${job_name}.sbatch <assembly_file> [tree]

Requires: 
 - R
 - ape v5.7-1 (https://cran.r-project.org/package=ape)
 
Default Newick tree: Kelp phylogeny from Starko, S. et al. 2019 ($tree)"
  exit 0
fi

# Define variables
seqFile="${scripts_dir}s_lat_alignment.txt"
r_script="${scripts_dir}${job_name}.R"
assembly_file="$1"
(( $# > 1 )) && tree="$2"
if [[ ! -f "$tree" ]] && [[ "$tree" == http?(s)://*.txt ]]
then
    if [[ -f $(basename "$tree") ]]
    then
      tree=$(basename "$tree")
    else
      wget "$tree"
      tree=$(basename "$tree")
    fi
fi
echo "Assembly file: $assembly_file"
echo "Input Newick tree: $tree"
[[ -f "$assembly_file" ]] || { echo "Error: Assembly file not found."; exit 1; }
[[ -f "$tree" ]] || { echo "Error: Tree not found."; exit 1; }

# Run tree pruning R script
# Generates formatted seqFile, i.e.,
# NEWICK tree
# name1 path1
# name2 path2
# ...
# nameN pathN
cmd="Rscript $r_script $assembly_file $tree $seqFile $scripts_dir"
echo "$cmd"
$cmd
