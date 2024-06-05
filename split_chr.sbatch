#!/bin/bash
#SBATCH -J split_chr
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

Usage: sbatch ${job_name}.sbatch <assembly.fa> <chr_ID> <minsize>

Requires:
 - Samtools (https://www.htslib.org/)
 - miscperlscripts (https://github.com/MadsAlbertsen/miscperlscripts/)"
  exit 0
fi

# Input and output variables
assembly="$1"
# ID of chromosome to split
chr_ID="$2"
# Minimum size of produced contigs (S. japonica=500bp, E. siliculosus=2002bp)
minsize="$3"
assembly_idx="${assembly}.fai"
assembly_no_ext="${assembly%%.*}"
assembly_ext="${assembly##*.}"
assembly_no_chr_list="${assembly_no_ext}_no_${chr_ID}.list"
assembly_no_chr="${assembly_no_ext}_no_${chr_ID}.${assembly_ext}"
out1="${assembly_no_ext}_${chr_ID}.${assembly_ext}"
out2="${assembly_no_ext}_${chr_ID}_split.${assembly_ext}"
out3="${assembly_no_ext}_${chr_ID}_split_in_genome.${assembly_ext}"

# Activate Samtools
conda activate samtools
samtools --version | head -n1
# Index assembly
cmd="samtools faidx $assembly"
echo "$cmd"
$cmd
# Extract chr_ID from assembly
cmd="samtools faidx -o $out1 $assembly $chr_ID"
echo "$cmd"
$cmd
conda deactivate
# Subset all other IDs not chr_ID from assembly
awk '{print $1}' "$assembly_idx" | grep -v "$chr_ID" > "$assembly_no_chr_list"
conda activate seqtk
cmd="seqtk subseq $assembly $assembly_no_chr_list"
echo "$cmd > $assembly_no_chr"
$cmd > "$assembly_no_chr"
conda deactivate

# Split chr_ID into contigs
cmd="split.scaffolds.to.contigs.pl -i $out1 -o $out2 -m $minsize"
echo "$cmd"
$cmd

# Create full assembly FASTA version with split chr
cmd="cat $assembly_no_chr $out2"
echo "$cmd > $out3"
$cmd > "$out3"

# Index newly created assembly FASTA
conda activate samtools
cmd="samtools faidx $out3"
echo "$cmd"
$cmd

# Clean up intermediate files
cmd="rm $out1 $out2 $assembly_no_chr $assembly_no_chr_list"
echo "$cmd"
$cmd

# Results
echo "Results in: $out3"
