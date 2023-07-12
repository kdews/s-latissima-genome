#!/bin/bash
#SBATCH -J seqtk_subseq
#SBATCH --mem=10gb
#SBATCH --time=01:00:00
#SBATCH -o %x.log

# Print date and time
date

# Load conda
cond=~/.conda_for_sbatch.sh
[[ -a "$cond" ]] && source "$cond" || \
{ echo "Error on source of $cond"; exit 1; }

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Uses a list of FASTA file IDs to subset from a multi-fasta file. Output FASTA is
named after list and FASTA file prefixes.
Usage: sbatch seqtk_subseq.sbatch <list> <multifasta>
Requires: seqtk (https://github.com/lh3/seqtk)"
  exit 0
fi

# Set variables
seqlist="$1"
seqlist_no_ext="$(basename ${seqlist%.*})"
assembly="$2"
assembly_no_ext="$(basename ${assembly%.*})"
prefix="${seqlist_no_ext}_${assembly_no_ext}"
subset="${prefix}.fasta"

# Load conda environment for seqtk
conda activate seqtk

# Subset list of FASTA IDs into $outfa
cmd="seqtk subseq $assembly $seqlist"
echo "$cmd > $subset"
$cmd > $subset

# Separate contig multifasta into unique files
int=$(grep -c ">" "$subset")
cmd="seqtk split -n $int $prefix $subset"
echo "$cmd"
$cmd