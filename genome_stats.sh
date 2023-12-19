#!/bin/bash

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]] || (( $# < 2 ))
then
  echo "\
Runs BUSCO and QUAST on list of genome assemblies.
Usage:
  bash genome_stats.sh <assembly_file> [path/to/aug_busco.sbatch] [path/to/quast.sbatch]
  bash genome_stats.sh <assembly_file> [path/to/busco_compare.sbatch]
Requires: 
  - BUSCO (https://busco.ezlab.org)
  - QUAST (https://quast.sourceforge.net)"
  exit 0
fi

# Input
scripts_dir_name="s-latissima-genome"
[[ -d $scripts_dir_name ]] && scripts_dir="${scripts_dir_name}/"
# List of assemblies (s-latissima-genome/species_table.txt)
assembly_file="$1"
if [[ -f "$assembly_file" ]]
then
  echo "Using assembly file: $assembly_file"
  mapfile -t assembly_list < <(grep -v "#" "$assembly_file" | awk -F '\t' '{print $2}')
  mapfile -t annot_list < <(grep -v "#" "$assembly_file" | awk -F '\t' '{print $3}')
  echo "Found ${#assembly_list[@]} assemblies and ${#annot_list[@]} annotations."
else
  echo "Error: assembly file ($assembly_file) not found."
fi
# Scripts to run (all positional arguments after first one)
script_list=("${@:2}")
echo "Running ${#script_list[@]} scripts: ${script_list[*]}"

# Iterate through given genome assemblies
for i in "${!assembly_list[@]}"
do
  assembly="${assembly_list[i]}"
  annot="${annot_list[i]}"
  if [[ -f "$assembly" ]]
  then
    assembly_no_ext=$(basename "${assembly%.*}")
  else
    echo "Error: assembly ($assembly) not found."
    exit 1
  fi
  if ! [[ -f "$annot" ]]
  then
    echo "Error: annotation ($annot) not found."
    exit 1
  fi
  for script in "${script_list[@]}"
  do
    if [[ -f "$script" ]]
    then
      script_no_ext=$(basename "${script%.*}")
      if [[ "$script_no_ext" == "busco_compare" ]]
      then
        submission="sbatch $script $assembly_file"
        echo "$submission"
        $submission
        exit 0
      fi
    else
      echo "Error: script $script not found."
      exit 1
    fi
    log="-o ${script_no_ext}_${assembly_no_ext}.log"
    submission="sbatch $log $script $assembly $annot"
    echo "$submission"
    $submission
  done
done
