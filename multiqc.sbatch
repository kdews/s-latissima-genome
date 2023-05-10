#!/bin/bash
#SBATCH -J multiqc
#SBATCH -p main
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=05:00:00
#SBATCH -o %x.log

# Print date and time
date

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Generates MultiQC report from files in directories named ${prefixes[@]}...,
where by default prefixes=("busco" "quast").

Usage: sbatch busco_compare.sbatch [prefixes]
Requires: MultiQC (https://multiqc.info)"
  exit 0
fi

# Load conda
cond=~/.conda_for_sbatch.sh
[[ -a "$cond" ]] && source "$cond" || \
{ echo "Error on source of $cond"; exit 1; }

# Load multiqc
conda activate multiqc
multiqc --version

# Create text file to rename samples
samples_file="assemblies/pretty_names.txt"
mqc_samples="multiqc_$(basename $samples_file)"
if [[ -f $mqc_samples ]]
then
  echo "Found file for sample naming: $mqc_samples"
else
  echo "Creating sample name file: $mqc_samples"
  echo "MultiQC Names	Proper Names" > $mqc_samples
  awk -F '\t' -v OFS='\t' '{print $2, $1}' $samples_file | sed "s/#//g" | \
  sed "s/\.fasta//g" >> $mqc_samples
fi


# Arguments to MultiQC
opts="-f --sample-names $mqc_samples"
# Optional arguments
# Specify directory containing scripts
scripts_dir="s-latissima-genome"
[[ -d "$scripts_dir" ]] && scripts_dir="${scripts_dir}/"
# Configuration file
config="-c ${scripts_dir}multiqc_config.yaml"
# Prefixes of directories to search
prefixes=("aug_busco")
# Can optionally be specified by n positional arguments
(( $# > 0 )) && prefixes=${@}
prefixes=( "${prefixes[@]/%/*/short_*.txt}" )

# Run MultiQC
cmd="multiqc $opts $config ${prefixes[@]}"
echo "$cmd"
$cmd
