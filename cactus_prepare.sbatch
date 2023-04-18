#!/bin/bash
#SBATCH -J cactus_prepare
#SBATCH -p shared
#SBATCH --mem=0
#SBATCH --time=7-0
#SBATCH --cpus-per-task=16
#SBATCH -o %x.log

## NOTE: if oneweek Discovery partition, use: --constraint=xeon-2640v3

# Print date and time
date

# Load Python module and Cactus virtualenv
module purge
module load gcc/11.3.0 python/3.9.12 git/2.36.1
env_path="/home1/kdeweese/cactus-bin-v2.5.0/cactus_env/bin/activate"
[[ -a "$env_path" ]] && source "$env_path" || \
{ echo "Error on source of $env_path"; exit 1; }

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
ncor="${SLURM_CPUS_PER_TASK}"
seqfile_no_ext=$(basename "${seqfile%.*}")
tmp="${SLURM_JOB_NAME}_tmp"
js="${SLURM_JOB_NAME}_jobstore"
msgs="${SLURM_JOB_NAME}_msgs"
# Output HAL file name
hal="${seqfile_no_ext}.hal"
# Create local temporary work directory
mkdir -p $tmp

# Set SLURM arguments for Toil
export TOIL_SLURM_ARGS="-t 1-0 -q normal -p shared"
# Set name of partition to use for parallel Toil jobs
export TOIL_SLURM_PE="${SLURM_JOB_PARTITION}"

# Run Cactus multialignment
slurm_opts1="--batchSystem slurm --consCores $ncor --disableAutoDeployment"
slurm_opts2="--defaultMemory 4G"
logs="--workDir $tmp --logDebug --writeMessages $msgs --realTimeLogging"
opts="--binariesMode local --outHal $hal"
cmd="cactus-prepare-toil $opts $slurm_opts1 $logs file:./$js $seqfile"
echo "$cmd"
$cmd
