#!/bin/bash
#SBATCH -J split_chr0
#SBATCH --mem=25gb
#SBATCH --time=01:00:00
#SBATCH -o %x.log

# Initialize script
init="/home1/kdeweese/scripts/s-latissima-genome/script_init.sh"
if [[ -a "$init" ]]; then source "$init"; else { echo "Init err"; exit 1; }; fi

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Extracts and splits artificial chromosomes from ORCAE FASTAs.

Usage: sbatch ${job_name}.sbatch <assembly.fa>

Requires:
 - Samtools (https://www.htslib.org/)
 - miscperlscripts (https://github.com/MadsAlbertsen/miscperlscripts/)"
  exit 0
fi

# Input and output variables
assembly="$1"
assembly_base="$(basename "${assembly%%.*}")"
assembly_ext="$(basename "${assembly##*.}")"
out1="${assembly_base}_chr0.${assembly_ext}"
out2="${assembly_base}_chr0_split.${assembly_ext}"

# Extract chr0 from assembly
conda activate samtools
samtools --version | head n1
cmd="samtools faidx -o $out1 $assembly chr0"
echo "$cmd"
$cmd
conda deactivate

# Split chr0 scaffold into contigs
cmd="split.scaffolds.to.contigs.pl -i $out1 -o $out2"
echo "$cmd"
$cmd
