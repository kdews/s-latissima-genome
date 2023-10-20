#!/bin/bash
#SBATCH -J hal2paf
#SBATCH -p shared
#SBATCH --mem=20gb
#SBATCH --time=05:00:00
#SBATCH -o %x.log

# Print date and time
date

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Extracts PAF from given HAL file using hal2paf from HAL tools.
Usage: sbatch ${SLURM_JOB_NAME}.sbatch <inHal> [outPaf]

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
outPaf="${inHal_no_ext}.paf"

cmd="hal2paf --hdf5InMemory $inHal"
echo "$cmd > $outPaf"
$cmd > "$outPaf"

# # Extract non-ancestral genomes from HAL
# genomes=($(halStats --genomes "$inHal" | sed 's/ /\n/g' | grep -v "Anc"))

# # Run halSynteny for each combination of genomes
# for g1 in "${genomes[@]}"
# do
#   for g2 in "${genomes[@]}"
#   do
#     if [[ "$g1" != "$g2" ]]
#     then
#       outPsl="${inHal_no_ext}_${g1}_vs_${g2}.psl"
#       cmd="halSynteny --queryGenome $g1 --targetGenome $g2 $inHal $outPsl"
#       echo "$cmd"
#       $cmd
#     fi
#   done
# done