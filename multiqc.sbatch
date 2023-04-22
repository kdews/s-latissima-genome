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

# Optional arguments
# Specify directory containing scripts
scripts_dir="s-latissima-genome"
[[ -d "$scripts_dir" ]] && scripts_dir="${scripts_dir}/"
# Configuration file
config="-c ${scripts_dir}multiqc_config.yaml"
# Prefixes of directories to search
prefixes=("aug_busco" "quast")
# Can optionally be specified by n positional arguments
(( $# > 0 )) && prefixes=${@}
prefixes=( "${prefixes[@]/%/*/}" )

# Run MultiQC
cmd="multiqc $config $prefixes"
echo "$cmd"
$cmd
