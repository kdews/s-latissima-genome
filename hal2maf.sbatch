#!/bin/bash
#SBATCH -J hal2maf
#SBATCH -p shared
#SBATCH --mem=20gb
#SBATCH --time=05:00:00
#SBATCH -o %x.log

# Print date and time
date

# Initialize script
init="s-latissima-genome/script_init.sh"
if [[ -a "$init" ]]; then source "$init"; else { echo "Init err"; exit 1; }; fi

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Extracts MAF from given HAL file using hal2maf from HAL tools.
Usage: sbatch ${SLURM_JOB_NAME}.sbatch [inHal] [outMaf]

Requires:
 - HAL (https://github.com/ComparativeGenomicsToolkit/hal/)"
  exit 0
fi

# Load Python module and Cactus virtualenv (contains HAL)
module purge
module load gcc/12.3.0 python/3.11.4 git/2.42.0 py-pip/23.0
env_path="/home1/kdeweese/bin/cactus-bin-v2.6.7/cactus_env/bin/activate"
[[ -a "$env_path" ]] && source "$env_path" || \
{ echo "Error on source of $env_path"; exit 1; }

# Input variable(s)
inHal="cactus-steps-output/s_lat_alignment.hal"
[[ -z "$1" ]] || inHal="$1"
inHal_no_ext=$(basename "${inHal%.*}")
outMaf="${inHal_no_ext}.maf"

cmd="hal2maf --hdf5InMemory $inHal $outMaf"
echo "$cmd"
$cmd
