#!/bin/bash
#SBATCH --mem=70g
#SBATCH --cpus-per-task=4
#SBATCH --time=07:00:00
#SBATCH -J minimap2
#SBATCH -o logs/%x_%j.log

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
ref_list="s-latissima-genome/filt_species_table.txt"
# ref_list="$1"
ref="$(sed -n "${SLURM_ARRAY_TASK_ID}"p $ref_list | awk '{print $2}')"
query="ragout-out/Saccharina_latissima_scaffolds.fasta"
r_base="$(basename "$ref" | sed 's/\..*//g')"
q_base="$(basename "$query" | sed 's/\..*//g')"
base_name="${q_base}_vs_${r_base}"
outPaf="${base_name}.paf"
thr="$((SLURM_CPUS_PER_TASK - 1))"

# Load minimap2 conda environment
conda activate minimap2
# Verify program installation
minimap2 --version


cmd="minimap2 -t $thr -x asm20 $ref $query"
echo "$cmd > $outPaf"
$cmd > "$outPaf"
