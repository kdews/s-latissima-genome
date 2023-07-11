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
('busco_contigs.txt') of unfiltered contig IDs that contain BUSCO genes. If 
'busco_contigs.txt' exists, delete it to rerun the analysis.
Usage: sbatch detailed_busco_analysis.sbatch
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

# Count number of contigs in result
int="$(wc -l $contig_ids | awk '{print $1}')"

# Subset contigs and split into single FASTAs
cmd="sbatch --parsable ${scripts_dir}/seqtk_subseq.sbatch $contig_ids"
echo "$cmd"
job1="$($cmd)"

# Generate mums with MUMmer for contigs versus reference genome
# Set variables
ref="assemblies/Macpyr2_AssemblyScaffolds_Repeatmasked.fasta"
query="busco_contigs_SugarkelpHifiasm.raconPolished.fasta"
cmd="sbatch -d afterok:${job1} -n $int ${scripts_dir}/mummer.sbatch $ref $query"
echo "$cmd"
$cmd

