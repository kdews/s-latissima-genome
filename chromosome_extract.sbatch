#!/bin/bash
#SBATCH -J chromosome_extract
#SBATCH -p main
#SBATCH --mem=2gb
#SBATCH --time=01:00:00
#SBATCH -o %x.log

# Initialize script
init="s-latissima-genome/script_init.sh"
if [[ -a "$init" ]]; then source "$init"; else { echo "Init err"; exit 1; }; fi

# Script environment variables
for_seqtk="for_seqtk.txt"
filt_assemb="filt_species_table.txt"

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Extracts longest representative scaffolds of each genome, given a species of
interest and an outgroup species.
Usage: sbatch ${job_name}.sbatch <assembly_file> <spc_int> <out_spc>

Required input:
  assembly_file    tabulated species names and file PATHs
    Format:
    Species_name	Genome_PATH	Annotation_PATH	Proteins_PATH	Gene_Info_PATH
  spc_int          species of interest (string, formatted 'Genus_specificname')
  out_spc          outgroup species (string, formatted 'Genus_specificname')

Output:
  $for_seqtk (table grouping FASTAs and list of contigs for filtering)
  $filt_assemb (filtered assembly_file table)

Requires:
 - samtools (https://www.htslib.org)
 - R (https://www.r-project.org)
 - seqtk (https://github.com/lh3/seqtk)"
  exit 0
fi

# Input
# Redirect to $scripts_dir, if exists
r_script="${scripts_dir}faidx_plots.R"
for_seqtk="${scripts_dir}${for_seqtk}"
filt_assemb="${scripts_dir}${filt_assemb}"
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
    out_spc="Ectocarpus_siliculosus"
fi

# Load samtools conda env
conda activate samtools
samtools --version
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
module load gcc/11.3.0 openblas/0.3.20 r/4.3.1
Rscript --version
# Run Rscript to plot and filter contigs
cmd="Rscript --vanilla $r_script $assembly_file $for_seqtk $spc_int $out_spc $scripts_dir"
echo "$cmd"
$cmd
# if [[ -f "$for_seqtk" ]]
# then
#     echo "Found file $for_seqtk, skipping Rscript step."
# else
#     cmd="Rscript --vanilla $r_script $assembly_file $for_seqtk $spc_int $out_spc $scripts_dir"
#     echo "$cmd"
#     $cmd
# fi

# # Extract "chromosomes" from assembly FASTAs
# echo
# # Load conda environment for seqtk
# conda activate seqtk
# if [[ -f "$for_seqtk" ]]
# then
#     mapfile -t species < <(awk -F '\t' '{print $1}' "$for_seqtk")
#     mapfile -t assemblies < <(awk -F '\t' '{print $2}' "$for_seqtk")
#     mapfile -t contig_lists < <(awk -F '\t' '{print $3}' "$for_seqtk")
# else
#     echo "Error: no input file ($for_seqtk) found for seqtk step."
#     exit 1
# fi
# # Create subsetted assembly FASTAs
# for i in "${!assemblies[@]}"
# do
#     out="${contig_lists[i]%.*}.fasta"
#     echo "$out"
#     cmd="seqtk subseq ${assemblies[i]} ${contig_lists[i]}"
#     echo "$cmd > $out"
#     $cmd > "$out"
#     # Tabulate species and filtered FASTAs
#     if [[ "$i" -eq 0 ]]
#     then
#         paste <(echo "${species[i]}") <(echo "$out") > "$filt_assemb"
#     else
#         paste <(echo "${species[i]}") <(echo "$out") >> "$filt_assemb"
#     fi
# done
# echo "Filtered assembly file: $filt_assemb"