#!/bin/bash
#SBATCH -J split_artificial
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
Extracts and splits artificial chromosomes (chr 0) from ORCAE FASTAs.

Usage: sbatch ${job_name}.sbatch <species_table>

Requires:
 - samtools (https://www.htslib.org/)
 - miscperlscripts (https://github.com/MadsAlbertsen/miscperlscripts/)"
  exit 0
fi

# Species table file
assembly_file="$1"
# Patterns to search for species with artificial chromosomes (chromosome 0)
spc_pats=("Saccharina japonica" "Ectocarpus")
# Minimum size of produced contigs (S. japonica=500bp, E. siliculosus=2002bp)
minsizes=(500 2002)
for i in "${!spc_pats[@]}"
do
  # Iterate through species and associated min sizes
  spc_pat="${spc_pats[i]}"
  minsize="${minsizes[i]}"
  # Match current species to assembly in $assembly_file
  spc="$(grep "$spc_pat" "$assembly_file" | awk -F '\t' '{print $1}')"
  assembly="$(grep "$spc" "$assembly_file" | awk -F '\t' '{print $2}')"
  if echo "$assembly" | grep -q "split_artificial"
  then
    echo "Skipping split for $spc"
  else
    echo "Species: $spc"
    echo "Assembly: $assembly"
    echo "Minimum contig size: $minsize"

    # Index assembly (if needed) to find artificial chromosome ID
    assembly_idx="${assembly}.fai"
    if [[ -f "$assembly_idx" ]]
    then
      echo "Using found index: $assembly_idx"
    else
      # Activate samtools
      conda activate samtools
      samtools --version | head -n1
      # Index assembly
      cmd="samtools faidx $assembly"
      echo "$cmd"
      $cmd
      conda deactivate
    fi

    # Find chromosome IDs where only number present is 0
    chr_ID="$(awk -F '\t' '{print $1}' "$assembly_idx"| grep -v "[1-9]" | grep "0\+$")"
    echo "Artificial chromosome ID: $chr_ID"

    # Temp filenames
    assembly_no_ext="${assembly%%.*}"
    ext="${assembly##*.}"
    temp_assembly_file="${assembly_file}.temp"
    # Assembly FASTA minus chr_ID
    assembly_no_chr_list="${assembly_no_ext}_no_${chr_ID}.list"
    assembly_no_chr="${assembly_no_ext}_no_${chr_ID}.${ext}"
    # FASTA of just chr_ID to be split into contigs
    chr_subset="${assembly_no_ext}_subset_${chr_ID}.${ext}"
    chr_split="${assembly_no_ext}_subset_${chr_ID}_split.${ext}"
    # Final output FASTA with artificial chromosome split
    assembly_chr_split="${assembly_no_ext}_split_artificial.${ext}"

    # Extract chr_ID from assembly
    echo "Extracting $chr_ID"
    # Activate samtools
    conda activate samtools
    samtools --version | head -n1
    cmd="samtools faidx -o $chr_subset $assembly $chr_ID"
    echo "$cmd"
    $cmd
    conda deactivate

    # Subset all other IDs not chr_ID from assembly
    conda activate seqtk
    echo "seqtk $(seqtk |& grep -i "version")"
    awk '{print $1}' "$assembly_idx" | grep -v "^$chr_ID$" > "$assembly_no_chr_list"
    cmd="seqtk subseq $assembly $assembly_no_chr_list"
    echo "$cmd > $assembly_no_chr"
    $cmd > "$assembly_no_chr"
    conda deactivate

    # Split chr_ID into contigs
    echo "Splitting $chr_ID into contigs."
    cmd="split.scaffolds.to.contigs.pl -i $chr_subset -o $chr_split -m $minsize"
    echo "$cmd"
    $cmd

    # Create full assembly FASTA version with split chr and index
    cmd="cat $assembly_no_chr $chr_split"
    echo "$cmd > $assembly_chr_split"
    $cmd > "$assembly_chr_split"
    conda activate samtools
    cmd="samtools faidx $assembly_chr_split"
    echo "$cmd"
    $cmd
    conda deactivate

    # Clean up temp files
    echo "Removing temporary files."
    cmd="rm $chr_subset $chr_split $assembly_no_chr $assembly_no_chr_list"
    echo "$cmd"
    $cmd

    echo "Full assembly with $chr_ID split into contigs: $assembly_chr_split"

    # Create temporary species table file 
    if [[ -f "$temp_assembly_file" ]]
    then
      # Edit existing temp file if found
      sed -i "s|$assembly|$assembly_chr_split|g" "$temp_assembly_file"
    else
      sed "s|$assembly|$assembly_chr_split|g" "$assembly_file" > "$temp_assembly_file"
    fi
  fi
  echo
done
echo "Finished."

# Update species table with new assembly FASTA
if [[ -f "$temp_assembly_file" ]]
then
  echo
  echo "Previous species table version backed up in ${assembly_file}.og"
  cp "$assembly_file" "${assembly_file}.og"
  echo "Updating species table: $assembly_file"
  mv "$temp_assembly_file" "$assembly_file"
fi
