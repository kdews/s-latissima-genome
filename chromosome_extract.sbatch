#!/bin/bash
#SBATCH -J chromosome_extract
#SBATCH --mem=100mb
#SBATCH --time=01:00:00
#SBATCH -o %x.log

# Print date and time
date

# Load conda
cond=~/.conda_for_sbatch.sh
[[ -a "$cond" ]] && source "$cond" || \
{ echo "Error on source of $cond"; exit 1; }

# Index assemblies
# Input
r_script="s-latissima-genome/faidx_plots.R"
for_seqtk="for_seqtk.txt"
species_file="$1"
[[ -z "$species_file" ]] && species_file="species.txt"
fastas=($(awk '{print $2}' "$species_file"))
# Load samtools conda env
conda activate samtools
samtools --version
# Run samtools index tool on each FASTA
for fasta in "${fastas[@]}"
do
    cmd="samtools faidx $fasta"
    echo "$cmd"
    $cmd
done

# Plot and filter assembly indices
# Load R module
module purge
module load gcc/11.3.0 openblas/0.3.20 r/4.3.1
Rscript --version
# Run Rscript to plot and filter contigs
cmd="Rscript --vanilla $r_script $species_file $for_seqtk"
echo "$cmd"
$cmd

# Extract "chromosomes" from assembly FASTAs
# Load conda environment for seqtk
conda activate seqtk
# Create subsetted assembly FASTAs
assemblies=($(awk -F '\t' '{print $2}' "$for_seqtk"))
contig_lists=($(awk -F '\t' '{print $3}' "$for_seqtk"))
for i in "$(seq 0 $((${#assemblies[@]}-1)))"
do
    cmd="seqtk subseq ${assemblies[i]} ${contig_lists[i]}"
    echo "$cmd"
    $cmd
done