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
seqlist="busco_contigs.txt"
seqlist="$1"
seqlist_no_ext="$(basename ${seqlist%.*})"
fasta="assemblies/SugarkelpHifiasm.raconPolished.fasta"
fasta="$2"
fasta_no_ext="$(basename ${fasta%.*})"
prefix="${seqlist_no_ext}_${fasta_no_ext}"
outfa="${prefix}.fasta"

# Load conda environment for seqtk
conda activate seqtk

# Subset list of FASTA IDs into $outfa
cmd="seqtk subseq $fasta $seqlist"
echo "$cmd > $outfa"
$cmd > $outfa

# Separate contig multifasta into unique files
int=$(grep -c ">" "$outfa")
cmd="seqtk split -n $int $prefix $outfa"
echo "$cmd"
$cmd