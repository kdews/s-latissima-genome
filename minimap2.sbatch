#!/bin/bash
#SBATCH -p gpu
#SBATCH --mem=70g
#SBATCH --cpus-per-task=4
#SBATCH --time=07:00:00
#SBATCH -J minimap2
#SBATCH -o %x_%a.log
#SBATCH -a 1-5

# Initialize script
init="s-latissima-genome/script_init.sh"
if [[ -a "$init" ]]; then source "$init"; else { echo "Init err"; exit 1; }; fi

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Aligns two genomes with minimap2 into PAF file.

Usage: sbatch ${job_name}.sbatch <ref_list> <query>

Requires:
 - minimap2 (https://github.com/lh3/minimap2/)"
  exit 0
fi

# Input variables
ref_list="$1"
query="$2"
# Default input
[[ ! -f "$ref_list" ]] && ref_list="s-latissima-genome/filt_species_table.txt"
[[ ! -f "$query" ]] && query="ragout-out/Saccharina_latissima_scaffolds.fasta"
ref="$(sed -n "${SLURM_ARRAY_TASK_ID}"p "$ref_list" | awk -F '\t' '{print $2}')"
r_base="$(basename "${ref%%.*}")"
q_base="$(basename "${query%%.*}")"
base_name="${q_base}_vs_${r_base}"
outPaf="${base_name}.paf"
thr="$((SLURM_CPUS_PER_TASK - 1))"

# Load minimap2 conda environment
conda activate minimap2
# Verify program installation
echo "minimap2 version $(minimap2 --version)"
echo

# Run minimap2 genome-genome alignment
cmd="minimap2 -t $thr -x asm20 -o $outPaf $ref $query"
echo "$cmd"
$cmd
