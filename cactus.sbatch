#!/bin/bash
#SBATCH -J cactus
#SBATCH -p oneweek
#SBATCH --mem=150gb
#SBATCH --time=7-0
#SBATCH --cpus-per-task=16
#SBATCH -o %x.log

# Print date and time
date

# Load conda and Cactus env
cond=~/.conda_for_sbatch.sh
[[ -a "$cond" ]] && source "$cond" || \
{ echo "Error on source of $cond"; exit 1; }
conda activate cactus
which cactus

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Runs Cactus alignment of multiple genomes from different species to generate a HAL file.
Usage: sbatch ${SLURM_JOB_NAME}.sbatch
Requires: Cactus (https://github.com/ComparativeGenomicsToolkit/cactus"
  exit 0
fi

# Define variables
script_dir="s-latissima-genome"
[[ -d $script_dir ]] && script_dir="${script_dir}/"
seqfile="s_lat_alignment.txt"
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
export TOIL_SLURM_PE='oneweek'

# Run Cactus multialignment
opts="--batchSystem Slurm --binariesMode local --workDir $tmp --clean onError"
# slurm_opts="--maxNodes ${SLURM_CPUS_PER_TASK} --maxMemory ${SLURM_MEM_PER_NODE}"
cmd="cactus $opts $js ${script_dir}${seqfile} $hal"
echo "$cmd"
$cmd