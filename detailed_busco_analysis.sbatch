#!/bin/bash
#SBATCH -J detailed_busco_analysis
#SBATCH -p gpu
#SBATCH --mem=10gb
#SBATCH -c 4
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
Analysis of BUSCO results that plots a comparison of BUSCO scores across genomes
and highlights specific conserved BUSCO genes in heatmaps. Also outputs a list
('busco_contigs.txt') of unfiltered contig IDs that contain BUSCO genes.
Usage: sbatch detailed_busco_analysis.sbatch
Requires: BUSCO (https://busco.ezlab.org)"
  exit 0
fi

# Load modules needed for Rscript
scripts_dir="s-latissima-genome"
mods="$(cat ${scripts_dir}/modules_for_Rscript.txt)"
module load $mods
cmd="Rscript ${scripts_dir}/${SLURM_JOB_NAME}.R"
echo "$cmd"
$cmd

# Load conda environment for seqtk
conda activate seqtk
fasta="assemblies/SugarkelpHifiasm.raconPolished.fasta"
fasta_no_path="$(basename $fasta)"
seqlist="busco_contigs.txt"
seqlist_no_ext="$(basename ${seqlist%.*})"
outfa="${seqlist_no_ext}_${fasta_no_path}"
cmd="seqtk subseq $fasta $seqlist"
echo "$cmd > $outfa"
$cmd > $outfa