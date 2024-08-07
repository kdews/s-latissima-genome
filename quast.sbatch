#!/bin/bash
#SBATCH -J quast
#SBATCH -p oneweek
#SBATCH --mem=15gb
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=12
#SBATCH -o %x.log

# Initialize script
init="/home1/kdeweese/scripts/s-latissima-genome/script_init.sh"
if [[ -a "$init" ]]; then source "$init"; else { echo "Init err"; exit 1; }; fi

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Runs QUAST on a given assembly 
Usage: sbatch ${job_name}.sbatch <assembly.fasta>
Requires: QUAST (https://quast.sourceforge.net)"
  exit 0
fi

# Load QUAST env
conda activate quast
# Load GNU Parallel
module load parallel

# Define variables
assembly="$1"
gff="$2"
spc="$3"
assembly_no_ext="$(basename "${assembly%.*}")"
echo "Species: $spc"
# Shorten SLURM output variable(s)
nthrd="$SLURM_CPUS_PER_TASK"
# Output directory name
out="${job_name}_${assembly_no_ext}"

# Run QUAST
# quast="quast --large -k -b --report-all-metrics --threads $nthrd -o $out"
cmd=(
  "quast"
  "--large"
  "-k"
  "-s"
  "-f"
  "--threads"
  "$nthrd"
  "-o"
  "$out"
  "-g"
  "$gff"
  "-l"
  "$spc"
  "$assembly"
)
echo "${cmd[*]}"
"${cmd[@]}"
