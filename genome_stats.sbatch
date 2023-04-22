#!/bin/bash

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Runs BUSCO and QUAST on list of genome assemblies.
Usage: bash genome_stats.sbatch
Requires: 
  - BUSCO (https://busco.ezlab.org)
  - QUAST (https://quast.sourceforge.net)"
  exit 0
fi


# Iterate through all genome assemblies in directory
for assembly in $(ls assemblies | grep "\.fasta")
do
  if [[ -f assemblies/"$assembly" ]]
  then
    assembly_no_ext=$(basename "${assembly%.*}")
  else
    echo "Error: assembly $assembly not found."
    exit 1
  fi
  for script in "$@"
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
    submission="sbatch $log $script assemblies/$assembly"
    echo "$submission"
    $submission
  done
done
