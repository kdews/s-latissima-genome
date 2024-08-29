#!/bin/bash
#SBATCH -J scaffold_eval
#SBATCH -p main
#SBATCH --mem=2gb
#SBATCH --time=01:00:00
#SBATCH -o %x.log

# Initialize script
init="/home1/kdeweese/scripts/s-latissima-genome/script_init.sh"
if [[ -a "$init" ]]; then source "$init"; else { echo "Init err"; exit 1; }; fi

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Given a species of interest and an outgroup species, generates plots comparing 
assemblies on basis of scaffolding.
Usage: sbatch ${job_name}.sbatch <assembly_file> <spc_int> <out_spc>

Required input:
  assembly_file    tabulated species names and file PATHs
    Format:
    Species_name	Genome_PATH	Annotation_PATH	Proteins_PATH	Gene_Info_PATH
  spc_int          species of interest (string, formatted 'Genus_specificname')
  out_spc          outgroup species (string, formatted 'Genus_specificname')

Requires:
 - samtools (https://www.htslib.org)
 - R (https://www.r-project.org)"
  exit 0
fi

# Input
# Redirect to $scripts_dir, if exists
r_script="${scripts_dir}${job_name}.R"
# Species table file
assembly_file="$1"
# Species of interest
spc_int="$2"
# Outgroup species
out_spc="$3"
if [[ -z "$assembly_file" ]] 
then
    assembly_file="${scripts_dir}species_table.txt"
fi
if [[ -z "$spc_int" ]] 
then
    spc_int="Saccharina_latissima"
fi
if [[ -z "$out_spc" ]] 
then
    out_spc="Ectocarpus_sp."
fi

# Load samtools conda env
conda activate samtools
samtools --version | head -n4
# Run samtools index tool on each FASTA
mapfile -t fastas < <(awk -F '\t' '{print $2}' "$assembly_file")
# Index assemblies
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
echo
# Load R module
module purge
$load_R_mods
Rscript --version
# Run Rscript to plot and filter contigs
cmd="Rscript --vanilla $r_script $assembly_file $spc_int $out_spc $scripts_dir"
echo "$cmd"
$cmd
