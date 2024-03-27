#!/bin/bash
#SBATCH -J halSynteny
#SBATCH -p oneweek
#SBATCH --mem=5gb
#SBATCH --time=05:00:00
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=1
#SBATCH -o %x.log

# Initialize script
init="s-latissima-genome/script_init.sh"
if [[ -a "$init" ]]; then source "$init"; else { echo "Init err"; exit 1; }; fi

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Runs halSynteny on given HAL file and outputs PSL formatted files of synteny
blocks for each combination of genomes.
Usage: sbatch ${job_name}.sbatch <inHal>

Requires:
 - HAL tools (https://github.com/ComparativeGenomicsToolkit/hal/)"
  exit 0
fi

# Load Python module and Cactus virtualenv (contains HAL tools)
module purge
# module load gcc/12.3.0 python git py-pip
module load gcc/11.3.0 python git py-pip parallel
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
inHal_no_ext="$(basename "${inHal%.*}")"

# Extract non-ancestral genomes from HAL
mapfile -t gen < <(halStats --genomes "$inHal" | sed 's/ /\n/g' | grep -v "Anc")

# Run halSynteny for each combination of genomes
# Count number of jobs that will be produced
parallel --plus "echo {uniq}" ::: "${gen[@]}" ::: "${gen[@]}"
njobs="$(parallel --plus "echo {uniq}" ::: "${gen[@]}" ::: "${gen[@]}" | wc -l)"
# Use GNU parallel to parallelize jobs
pll="parallel --plus -j $njobs --joblog ${job_name}.joblog"
echo "$pll"
# Toggle for parallel's "dry run" option to preview generated commands
# pll="$pll --dry-run"
# Use SLURM srun to schedule each task
sr="srun --exclusive -N 1 -n 1"
# Output PSL filename
outPsl="${inHal_no_ext}_{1uniq}_vs_{2uniq}.psl"
# halSynteny command
cmd="$sr halSynteny --queryGenome {1uniq} --targetGenome {2uniq} $inHal $outPsl"
# Generate all combinations of species with parallel's "uniq" replacement string
$pll "echo $cmd && $cmd" ::: "${gen[@]}" ::: "${gen[@]}"

