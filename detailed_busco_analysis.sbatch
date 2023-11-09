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
if [[ -a "$cond" ]]
then
  source "$cond"
else
  echo "Error on source of $cond"
  exit 1
fi

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Analysis of BUSCO results that plots a comparison of BUSCO scores across genomes
and highlights specific conserved BUSCO genes in heatmaps. Also outputs a list
('busco_contigs.txt') of unfiltered contig IDs that contain BUSCO genes. If 
'busco_contigs.txt' exists, delete it to rerun the analysis.
Usage: sbatch ${SLURM_JOB_NAME}.sbatch
Requires:
 - R (https://www.r-project.org)
 - seqtk (https://github.com/lh3/seqtk)
 - MUMmer (https://github.com/mummer4/mummer)"
  exit 0
fi

# Analysis of BUSCO results
# Set variables
scripts_dir="s-latissima-genome"
r_script="${SLURM_JOB_NAME}.R"
contig_ids="busco_contigs.txt"
assembly="assemblies/SugarkelpHifiasm.raconPolished.fasta"
# Only execute if result not present
if [[ -f "$contig_ids" ]]
then
  echo "Skipping $r_script run, $contig_ids found in working directory."
else
  # Load modules needed for Rscript
  mods="$(cat ${scripts_dir}/modules_for_Rscript.txt)"
  module load $mods
  # Run detailed BUSCO analysis
  cmd="Rscript ${scripts_dir}/$r_script $contig_ids"
  echo "$cmd"
  $cmd
fi
# Count resulting contigs
int="$(wc -l $contig_ids | awk '{print $1}')"

# Subset contigs
# Set variables
contig_ids_no_ext="$(basename ${contig_ids%.*})"
assembly_no_ext="$(basename ${assembly%.*})"
prefix="${contig_ids_no_ext}_${assembly_no_ext}"
subset="${prefix}.fasta"
# Only execute if result not present
if [[ -f $subset ]]
then
  echo "Skipping subsetting step, $subset exists in working directory."
else
  # Load conda environment for seqtk
  conda activate seqtk
  # Subset list of FASTA IDs into multi-FASTA
  cmd="seqtk subseq $assembly $contig_ids"
  echo "$cmd > $subset"
  $cmd > $subset
fi

# Generate mums with MUMmer for contigs versus reference genome
# Set variables
ref="assemblies/Macpyr2_AssemblyScaffolds_Repeatmasked.fasta"
# Run all contigs versus entire genome
cmd="sbatch -a 0 ${scripts_dir}/mummer.sbatch $ref $subset"
echo "$cmd"
$cmd
# Set array size based on number of contigs
cmd="sbatch -a 1-$int ${scripts_dir}/mummer.sbatch $ref $assembly $contig_ids"
echo "$cmd"
$cmd

