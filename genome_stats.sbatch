#!/bin/bash

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Runs BUSCO and QUAST on list of genome assemblies.
Usage: bash genome_stats.sbatch <assembly_file>
Requires: 
  - BUSCO (https://busco.ezlab.org)
  - QUAST (https://quast.sourceforge.net)"
  exit 0
fi

# List of assemblies
if [[ -z "$1" ]]
then
  assembly_file="assemblies/species_table.txt"
else
  assembly_file="$1"
fi
if [[ -f "$assembly_file" ]]
then
  assembly_list=( $(grep -v "#" "$assembly_file" | awk -F '\t' '{print $2}') )
  annot_list=( $(grep -v "#" "$assembly_file" | awk -F '\t' '{print $3}') )
else
  echo "Error: assembly file ($assembly_file) not found."
fi
script_dir="s-latissima-genome"
script_list=( "aug_busco.sbatch" "quast.sbatch" )
# script_list=( "busco_compare.sbatch" )
# Prepend script directory to all script names
script_list=( "${script_list[@]/#/$script_dir/}" )

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
        submission="sbatch $script busco_summaries"
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
