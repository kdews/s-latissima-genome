#!/bin/bash

# Initialize script
init="/home1/kdeweese/scripts/s-latissima-genome/script_init.sh"
if [[ -a "$init" ]]; then source "$init"; else { echo "Init err"; exit 1; }; fi

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]] || (( $# < 2 ))
then
  echo "\
Runs BUSCO and QUAST on list of genome assemblies.
Usage:
  bash ${job_name}.sh <assembly_file> [path/to/aug_busco.sbatch] [path/to/quast.sbatch]
  bash ${job_name}.sh <assembly_file> [path/to/busco_compare.sbatch]
Requires:
  - BUSCO (https://busco.ezlab.org)
  - QUAST (https://quast.sourceforge.net)
  - R (https://www.r-project.org)"
  exit 0
fi

# Input
# List of assemblies (s-latissima-genome/species_table.txt)
assembly_file="$1"
if [[ -f "$assembly_file" ]]
then
  echo "Using assembly file: $assembly_file"
  mapfile -t spc_list < <(grep -v "#" "$assembly_file" | awk -F '\t' '{print $1}')
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
  spc="${spc_list[i]}"
  # spc="${spc// /_}"
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
    log="${script_no_ext}_${assembly_no_ext}.log"
    submission=("sbatch" "-o" "$log" "$script" "$assembly" "$annot" "$spc")
    echo "${submission[*]}"
    "${submission[@]}"
  done
done
