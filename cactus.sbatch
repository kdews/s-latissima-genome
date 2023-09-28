#!/bin/bash
#SBATCH -J cactus
#SBATCH -p shared
#SBATCH --mem=0
#SBATCH --time=7-0
#SBATCH -o %x.log

## NOTE: if oneweek Discovery partition, use: --constraint=xeon-2640v3

# Print date and time
date

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Runs Cactus alignment of multiple genomes from different species to generate a
HAL file. Can create a Cactus-formatted seqFile from a Newick tree and species
list, if provided.
Usage: 
  Base: sbatch cactus.sbatch [species] [tree]
  Slurm-enabled: sbatch -J cactus_slurm --mem=1gb cactus.sbatch [species] [tree]

Requires:
 - virtualenv (https://virtualenv.pypa.io/en/latest)
 - Toil (https://toil.readthedocs.io/en/latest)
 - Cactus (https://github.com/ComparativeGenomicsToolkit/cactus)"
  exit 0
fi

# Load Python module and Cactus virtualenv
module purge
module load gcc/12.3.0 python/3.11.4 git/2.42.0 py-pip/23.0
# module load gcc/11.3.0 python/3.9.12 git/2.36.1
env_path="/home1/kdeweese/bin/cactus-bin-v2.6.7/cactus_env/bin/activate"
[[ -a "$env_path" ]] && source "$env_path" || \
{ echo "Error on source of $env_path"; exit 1; }

# Input file(s)
seqfile="s_lat_alignment.txt"
# If needed, generates formatted seqFile, i.e.,
# NEWICK tree
#
# name1 path1
# name2 path2
# ...
# nameN pathN
(( $# > 0 )) && {  species="$1"; tree="$2"; }
if [[ -f $seqfile ]]
then
  echo "Found seqfile: $seqfile"
elif [[ -v tree ]] && [[ -v species ]]
then
  echo "Species list: $species."
  echo "Newick tree: $tree"
  [[ -f $species ]] || { echo "Error: species list not found."; exit 1; }
  [[ -f $tree ]] || { echo "Error: tree not found."; exit 1; }
  echo "Generating formatted seqFile ($seqfile)..."
  cp $tree $seqfile
  echo "" >> $seqfile
  cat $species >> $seqfile
else
  echo "Error: no seqFile ($seqfile), tree ($tree) or species list ($species)."
  exit 1
fi

# Define variables
js="${SLURM_JOB_NAME}_jobstore"
seqfile_no_ext=$(basename "${seqfile%.*}")
hal="${seqfile_no_ext}.hal"

# Run Cactus multialignment
if [[ "$SLURM_JOB_NAME" == "cactus" ]]
then
  # Base Cactus
  cmd="cactus $js $seqfile $hal --realTimeLogging"
  echo "$cmd"
  $cmd
elif [[ "$SLURM_JOB_NAME" == "cactus_slurm" ]]
then
  # SLURM-enabled Cactus
  # Set SLURM arguments for Toil
  export TOIL_SLURM_ARGS="-t 1-0 -q normal -p shared"
  # Set name of partition to use for parallel Toil jobs
  export TOIL_SLURM_PE="${SLURM_JOB_PARTITION}"
  # max memory per node on shared partition = 60 gb
  slurm_opts1="--batchSystem slurm --binariesMode local"
  slurm_opts2="--consMemory 60gb --consCores 16 --configFile config-slurm.xml"
  morelogs="--logDebug"
  nohardlink="--noMoveExports --noLinkImports"
  toilopts="--disableCaching --disableAutoDeployment"
  res="--restart"
  # opts="$slurm_opts1 $slurm_opts2 $morelogs $nohardlink $toilopts $res"
  opts="--realTimeLogging $slurm_opts1 $slurm_opts2"
  cmd="cactus $js $seqfile $hal $opts"
  echo "$cmd"
  $cmd
else
  echo "Error: job name (${SLURM_JOB_NAME}) not recognized."
fi


# Extra options
# Create local temporary work directory
# tmp="${SLURM_JOB_NAME}_tmp"
# mkdir -p $tmp
# tmp_opt="--workDir $tmp"
# might lead to max allocated mem errors: --maxLocalJobs 100
# msgs="${SLURM_JOB_NAME}_msgs"
# logs2="--writeMessages $msgs --disableCaching --noStdOutErr --writeLogs $wklg"
# cleans="--clean never --stats --cleanWorkDir never"