#!/bin/bash
#SBATCH -J chromosome_extract
#SBATCH -p main
#SBATCH --mem=2gb
#SBATCH --time=01:00:00
#SBATCH -o %x.log

# Print date and time
date
echo

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
scripts_dir_name="s-latissima-genome"
[[ -d $scripts_dir_name ]] && scripts_dir="${scripts_dir_name}/"
r_script="${scripts_dir}faidx_plots.R"
# Species table file
assembly_file="$1"
# Species of interest
spc_int="$2"
if [[ -z "$assembly_file" ]] 
then
    assembly_file="${scripts_dir}species_table.txt"
fi
if [[ -z "$spc_int" ]] 
then
    spc_int="Saccharina_latissima"
fi
# Output
for_seqtk="${scripts_dir}for_seqtk.txt"
filt_assemb="${scripts_dir}filt_species_table.txt"

# Load samtools conda env
conda activate samtools
samtools --version
# Run samtools index tool on each FASTA
mapfile -t fastas < <(awk -F '\t' '{print $2}' "$assembly_file")
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
if [[ -f "$for_seqtk" ]]
then
    echo "Found file $for_seqtk, skipping Rscript step."
else
    cmd="Rscript --vanilla $r_script $assembly_file $for_seqtk $spc_int $scripts_dir"
    echo "$cmd"
    $cmd
fi

# Extract "chromosomes" from assembly FASTAs
echo
# Load conda environment for seqtk
conda activate seqtk
if [[ -f "$for_seqtk" ]]
then
    mapfile -t species < <(awk -F '\t' '{print $1}' "$for_seqtk")
    mapfile -t assemblies < <(awk -F '\t' '{print $2}' "$for_seqtk")
    mapfile -t contig_lists < <(awk -F '\t' '{print $3}' "$for_seqtk")
else
    echo "Error: no input file ($for_seqtk) found for seqtk step."
    exit 1
fi
# Create subsetted assembly FASTAs
for i in "${!assemblies[@]}"
do
    out="${contig_lists[i]%.*}.fasta"
    echo "$out"
    cmd="seqtk subseq ${assemblies[i]} ${contig_lists[i]}"
    echo "$cmd > $out"
    $cmd > "$out"
    # Tabulate species and filtered FASTAs
    if [[ "$i" -eq 0 ]]
    then
        paste <(echo "${species[i]}") <(echo "$out") > "$filt_assemb"
    else
        paste <(echo "${species[i]}") <(echo "$out") >> "$filt_assemb"
    fi
done
echo "Filtered assembly file: $filt_assemb"