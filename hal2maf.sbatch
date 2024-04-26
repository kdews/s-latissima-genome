#!/bin/bash
#SBATCH -J hal2maf
#SBATCH -p gpu
#SBATCH --mem=10gb
#SBATCH --time=05:00:00
#SBATCH -o %x.log

# Initialize script
init="/home1/kdeweese/scripts/s-latissima-genome/script_init.sh"
if [[ -a "$init" ]]; then source "$init"; else { echo "Init err"; exit 1; }; fi

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Extracts MAF from given HAL file using hal2maf from HAL tools.
Usage: sbatch ${job_name}.sbatch <inHal> [outMaf]

Requires:
 - HAL (https://github.com/ComparativeGenomicsToolkit/hal/)"
  exit 0
fi

# Input
inHal="cactus-steps-output/s_lat_alignment.hal"
[[ -n "$1" ]] && inHal="$1"
[[ -f "$inHal" ]] || { echo "Error: HAL ($inHal) not found."; exit 1; }
inHal_no_ext=$(basename "${inHal%.*}")
outMaf="${inHal_no_ext}.maf"
[[ -n "$2" ]] && outMaf="$2"

# Convert input HAL to MAF file
# Load Python module and Cactus virtualenv (contains HAL)
module purge
$load_cactus_mods
source_file "$cactus_env"
which hal2maf
echo
# Run hal2maf
cmd="hal2maf --hdf5InMemory $inHal $outMaf"
echo "$cmd"
$cmd
# Clear Cactus env (for pipeline)
deactivate
