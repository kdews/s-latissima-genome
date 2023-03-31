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
Runs Cactus alignment of multiple genomes from different species to generate a HAL file.
Usage: sbatch ${SLURM_JOB_NAME}.sbatch
Requires: 
 - virtualenv (https://virtualenv.pypa.io/en/latest)
 - Toil (https://toil.readthedocs.io/en/latest)
 - Cactus (https://github.com/ComparativeGenomicsToolkit/cactus"
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
export TOIL_SLURM_PE='main'

# Run Cactus multialignment
opts="--batchSystem slurm --binariesMode local --workDir $tmp --clean onError"
slurm_opts="--consCores ${SLURM_CPUS_PER_TASK}"
cmd="cactus $slurm_opts $opts $js ${script_dir}${seqfile} $hal"
echo "$cmd"
$cmd