#!/bin/bash
#SBATCH -J multiqc
#SBATCH --mem=2gb
#SBATCH --time=05:00:00
#SBATCH -o %x.log

# Initialize script
init="s-latissima-genome/script_init.sh"
if [[ -a "$init" ]]; then source "$init"; else { echo "Init err"; exit 1; }; fi

# Prefixes of directories to search
prefixes=("aug_busco*/short_*" "quast*/report*")
# prefixes=("aug_busco*/short_*.txt" "quast*/report*")

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Generates MultiQC report from files in directories named <prefix>….
Usage: sbatch ${job_name}.sbatch [prefixes]
Options:
 prefixes     space-separated list of prefixes to use when searching for
              directories to be used in analysis, named with the convention
              <prefix>… (default: ${prefixes[*]})
Requires:
 - MultiQC (https://multiqc.info)"
  exit 0
fi

# Load multiqc
conda activate multiqc
multiqc --version

# Create text file to rename samples
samples_file="${scripts_dir}species_table.txt"
mqc_samples="${scripts_dir}multiqc_$(basename "$samples_file")"
if [[ -f $mqc_samples ]]
then
  echo "Found file for sample naming: $mqc_samples"
else
  echo "Creating sample name file: $mqc_samples"
  echo "MultiQC Names	Proper Names" >> "$mqc_samples"
  # awk -F '\t' -v OFS='\t' '{print ".*"$2, $1}' $samples_file | sed "s/#//g" | \
  awk -F '\t' -v OFS='\t' '{print $2, $1}' "$samples_file" | sed "s/#//g" | \
  sed "s/\.fasta//g" >> "$mqc_samples"
fi


# Arguments to MultiQC
opts="-f"
# opts2="--replace-names $mqc_samples"
opts2="--sample-names $mqc_samples"
# Optional arguments
# Configuration file
config="-c ${scripts_dir}multiqc_config.yaml"
# Can optionally specify prefixes with positional arguments
(( $# > 0 )) && prefixes=("${@}")
# prefixes=( "${prefixes[@]/%/*/short_*.txt}" "${prefixes[@]/%/*/report*}")
# prefixes=( "${prefixes[@]/%/*/}" )

# Run MultiQC
# cmd="multiqc $opts $opts2 $config ${prefixes[@]}"
cmd="multiqc $opts $opts2 ${prefixes[*]}"
echo "$cmd"
$cmd
