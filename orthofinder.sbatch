#!/bin/bash
#SBATCH -J orthofinder
#SBATCH -p gpu
#SBATCH --mem=25gb
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=32
#SBATCH -o %x.log

# Initialize script
init="s-latissima-genome/script_init.sh"
if [[ -a "$init" ]]; then source "$init"; else { echo "Init err"; exit 1; }; fi

# Script variables
indir="ortho_prots"

# Help message
if [[ "$1" = "-h" ]] || [[ "$1" = "--help" ]]
then
  echo "\
Takes as input protein FASTAs intended for OrthoFinder.
Usage: sbatch ${job_name}.sbatch <assembly_file>

Requires:
 - OrthoFinder (https://github.com/davidemms/OrthoFinder)"
  exit 0
fi

# Load OrthoFinder conda env
conda activate orthofinder
orthofinder --help | grep "version"
echo

# Parallelization
# Low-mem threads
thr_t="$SLURM_CPUS_PER_TASK"
# High-mem threads
thr_a="$(( SLURM_CPUS_PER_TASK / 8 ))"

# Input
# List of assemblies (s-latissima-genome/species_table.txt)
assembly_file="$1"
if [[ -f "$assembly_file" ]]
then
  echo "Using assembly file: $assembly_file"
  mapfile -t prot_list < <(grep -v "#" "$assembly_file" | awk -F '\t' '{print $4}')
  echo "Found ${#prot_list[@]} protein FASTAs."
else
  echo "Error: assembly file ($assembly_file) not found."
fi

# Create input directory
mkdir -p -v "$indir"
# Copy protein FASTAs to $indir
for i in "${!prot_list[@]}"
do
  prot="${prot_list[i]}"
  if [[ -f "$prot" ]]
  then
    # Handle protein FASTA extension not recognized by OrthoFinder (.tfa)
    prot_ext="${prot##*.}"
    prot_base="$(basename "${prot%%."$prot_ext"}")"
    if [[ "$prot_ext" == "tfa" ]]
    then
      new_prot="${indir}/${prot_base}.faa"
    else
      new_prot="${indir}/${prot_base}.${prot_ext}"
    fi
    cp -v -u "$prot" "$new_prot"
  fi
done

# Run OrthoFinder
cmd="orthofinder -t $thr_t -a $thr_a -f $indir"
echo "$cmd"
$cmd
