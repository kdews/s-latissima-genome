#!/bin/bash
#SBATCH -J halSynteny
#SBATCH -p htcondor
#SBATCH --mem=20gb
#SBATCH --time=1-0
#SBATCH -o %x.log

# Print date and time
date

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Runs halSynteny on given HAL file and outputs a PSL formatted file of synteny
blocks.
Usage: sbatch ${SLURM_JOB_NAME}.sbatch <inHal> [outPsl]

Requires:
 - HAL (https://github.com/ComparativeGenomicsToolkit/hal/)"
  exit 0
fi

# Load Python module and Cactus virtualenv (contains HAL)
module purge
module load gcc/12.3.0 python/3.11.4 git/2.42.0 py-pip/23.0
env_path="/home1/kdeweese/bin/cactus-bin-v2.6.7/cactus_env/bin/activate"
if [[ -a "$env_path" ]]
then
  source "$env_path"
else
  echo "Error on source of $env_path"
  exit 1
fi

# Input variable(s)
inHal="cactus-steps-output/s_lat_alignment.hal"
[[ -z "$1" ]] || inHal="$1"
inHal_no_ext=$(basename "${inHal%.*}")

# Extract non-ancestral genomes from HAL
genomes=($(halStats --genomes "$inHal" | sed 's/ /\n/g' | grep -v "Anc"))

# Run halSynteny for each combination of genomes
for g1 in "${genomes[@]}"
do
  for g2 in "${genomes[@]}"
  do
    if [[ "$g1" != "$g2" ]]
    then
      outPsl="${inHal_no_ext}_${g1}_vs_${g2}.psl"
      cmd="halSynteny --queryGenome $g1 --targetGenome $g2 $inHal $outPsl"
      echo "$cmd"
      $cmd
    fi
  done
done