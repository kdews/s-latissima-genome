#!/bin/bash
#SBATCH -J busco_compare
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=05:00:00
#SBATCH -o %x.log

# Initialize script
init="/home1/kdeweese/scripts/s-latissima-genome/script_init.sh"
if [[ -a "$init" ]]; then source "$init"; else { echo "Init err"; exit 1; }; fi

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Generates plot(s) with conda-installed BUSCO generate_plot.py.
Usage: sbatch ${job_name}.sbatch <assembly_file>
Requires: BUSCO (https://busco.ezlab.org)"
  exit 0
fi

# Define variables
r_script="${scripts_dir}${job_name}.R"
R_mod_file="${scripts_dir}modules_for_Rscript.txt"
sum_dir="busco_summaries"
lineages=("eukaryota_odb10" "stramenopiles_odb10")
if [[ -z "$1" ]]
then
  assembly_file="${scripts_dir}/species_table.txt"
else
  assembly_file="$1"
fi
if [[ -f "$assembly_file" ]]
then
  mapfile -t assemblies < <(grep -v "#" "$assembly_file" | awk -F '\t' '{print $2}')
else
  echo "Error: assembly file ($assembly_file) not found."
fi

# Compile BUSCO summary text files (if needed)
mkdir -p "$sum_dir"
for lineage in "${lineages[@]}"
do
  wd="${sum_dir}/${lineage}"
  mkdir -p "$wd"
  # Clear any old text files from working directory (if they exist)
  [[ $(ls "${wd}"/*.txt 2>/dev/null) ]] && rm "${wd}"/*.txt
  # Copy BUSCO results for each lineage and assembly into summary directory
  for assembly in "${assemblies[@]}"
  do
    assembly_base=$(basename "${assembly%.*}")
    busco_base="aug_busco_${assembly_base}_${lineage}"
    sum_file="$busco_base/short_summary.specific.${lineage}.${busco_base}.txt"
    # Clean up BUSCO summary filename
    new="${sum_file//*aug_busco_/}"
    new_file="short_summary.specific.${lineage}.${new}"
    if [[ -f "$sum_file" ]]
    then
      cmd="cp -u $sum_file ${wd}/${new_file}"
      echo "$cmd"
      $cmd
    else
      pwd
      echo "Error: $lineage BUSCO summary file ($sum_file) not found for $assembly_base"
      exit 1
    fi
  done
  # Run BUSCO's generate_plot.py on each lineage
  # Activate BUSCO conda env
  conda activate busco
  cmd="generate_plot.py -wd $wd --no_r"
  echo "$cmd"
  $cmd
  # Deactivate env
  conda activate
  # Clean up BUSCO graphs
  # Load modules for Rscript
  IFS=" " read -r -a R_mods <<< "$(cat "$R_mod_file")"
  module load "${R_mods[@]}"
  cmd="Rscript $r_script $wd $assembly_file $scripts_dir"
  echo "$cmd"
  $cmd
  # Unload modules
  module purge
done
