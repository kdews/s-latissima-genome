#!/bin/bash
#SBATCH -J chromosome_extract
#SBATCH -p main
#SBATCH --mem=2gb
#SBATCH --time=01:00:00
#SBATCH -o %x.log

# Print date and time
date

# Load conda
cond=~/.conda_for_sbatch.sh
if [[ -a "$cond" ]]
then
  source "$cond"
else
  echo "Error on source of $cond"
  exit 1
fi

# Index assemblies
# Input
r_script="s-latissima-genome/faidx_plots.R"
for_seqtk="for_seqtk.txt"
# Species table file
species_file="$1"
# Species of interest
spc_int="$2"
if [[ -z "$species_file" ]] 
then
    species_file="assemblies/species_table.txt"
fi
if [[ -z "$spc_int" ]] 
then
    spc_int="Saccharina_latissima"
fi

# fastas=( $(awk '{print $2}' "$species_file") )
mapfile -t fastas < <(awk -F '\t' '{print $2}' "$species_file")
# Load samtools conda env
conda activate samtools
samtools --version
# Run samtools index tool on each FASTA
for fasta in "${fastas[@]}"
do
    if [[ -f "${fasta}.fai" ]]
    then
        echo "Found file ${fasta}.fai, skipping samtools faidx step."
    else
        cmd="samtools faidx $fasta"
        echo "$cmd"
        $cmd
    fi
done

# Plot and filter assembly indices
# Load R module
module purge
module load gcc/11.3.0 openblas/0.3.20 r/4.3.1
Rscript --version
# Run Rscript to plot and filter contigs
if [[ -f "$for_seqtk" ]]
then
    echo "Found file $for_seqtk, skipping Rscript step."
else
    cmd="Rscript --vanilla $r_script $species_file $for_seqtk $spc_int"
    echo "$cmd"
    $cmd
fi

# Extract "chromosomes" from assembly FASTAs
# Load conda environment for seqtk
conda activate seqtk
# Create subsetted assembly FASTAs
assemblies=( $(awk -F '\t' '{print $2}' "$for_seqtk") )
contig_lists=( $(awk -F '\t' '{print $3}' "$for_seqtk") )
for i in "${!assemblies[@]}"
do
    out="${contig_lists[i]%.*}.fasta"
    echo "$out"
    cmd="seqtk subseq ${assemblies[i]} ${contig_lists[i]}"
    echo "$cmd > $out"
    $cmd > "$out"
done