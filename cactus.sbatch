#!/bin/bash
#SBATCH -J cactus
#SBATCH -p oneweek
#SBATCH --constraint=xeon-2640v3
#SBATCH --mem=50gb
#SBATCH --time=7-0
#SBATCH --cpus-per-task=16
#SBATCH -o %x.log

# Print date and time
date

# Load Python module and Cactus virtualenv
module load gcc/11.3.0 python/3.9.12 git/2.36.1
env_path=/home1/kdeweese/bin/cactus-bin-v2.4.4/cactus_env/bin/activate
[[ -a "$env_path" ]] && source "$env_path" || \
{ echo "Error on source of $env_path"; exit 1; }

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Runs Cactus alignment of multiple genomes from different species to generate a
HAL file. Can create a Cactus-formatted seqFile from a Newick tree and species
list, if provided.
Usage: sbatch ${SLURM_JOB_NAME}.sbatch [newick_tree] [species_list]
Requires: 
 - virtualenv (https://virtualenv.pypa.io/en/latest)
 - Toil (https://toil.readthedocs.io/en/latest)
 - Cactus (https://github.com/ComparativeGenomicsToolkit/cactus"
  exit 0
fi


# If needed, generates formatted seqFile, i.e.,
# NEWICK tree
# name1 path1
# name2 path2
# ...
# nameN pathN
seqfile="s_lat_alignment.txt"
(( $# > 0 )) && { newick_tree="$1"; species_list="$2"; }
if [[ -f $seqfile ]]
then
  echo "Found seqfile: $seqfile"
elif [[ -v newick_tree ]] && [[ -v species_list ]]
  echo "Newick tree: ($newick_tree)"
  echo "Species list: ($species_list)."
  [[ -f $newick_tree ]] || { echo "Error: tree not found."; exit 1; }
  [[ -f $species_list ]] || { echo "Error: species list not found."; exit 1; }
  echo "Generating formatted seqFile ($seqfile)..."
  cp $newick_tree $seqfile
  cat $species_list >> $seqfile
else
  echo "\
Error: no seqfile ($seqfile), Newick tree ($newick_tree) or species list
($species_list) found."
  exit 1
fi


# Define variables
seqfile_no_ext=$(basename "${seqfile%.*}")
tmp="${SLURM_JOB_NAME}_tmp"
js="${SLURM_JOB_NAME}_jobstore"
# Output HAL file name
hal="${seqfile_no_ext}.hal"
# Create temporary work directory
mkdir -p $tmp

# Set SLURM arguments for Toil
export TOIL_SLURM_ARGS="-t 1-0 -q normal"
# Set name of partition to use for Toil
export TOIL_SLURM_PE='main'

# Run Cactus multialignment
opts="--batchSystem slurm --binariesMode local --workDir $tmp --clean onError"
slurm_opts="--consCores ${SLURM_CPUS_PER_TASK}"
cmd="cactus $slurm_opts $opts $js $seqfile $hal"
echo "$cmd"
# $cmd