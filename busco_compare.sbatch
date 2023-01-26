#!/bin/bash
#SBATCH -J busco_compare
#SBATCH -p gpu
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=05:00:00
#SBATCH -o %x.log

# Print date and time
date

# Load conda
cond=~/.conda_for_sbatch.sh
[[ -a "$cond" ]] && source "$cond" || \
{ echo "Error on source of $cond"; exit 1; }

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Generates plot(s) with conda-installed BUSCO generate_plot.py.
Usage: sbatch busco_compare.sbatch
Requires: BUSCO (https://busco.ezlab.org)"
  exit 0
fi


# Define variables
sum_dir="$1"
lineages=("eukaryota_odb10" "stramenopiles_odb10" "quast")


# Compile BUSCO summary text files (if needed)
mkdir -p "$sum_dir"
for lineage in "${lineages[@]}"
do
  wd="${sum_dir}/${lineage}/"
  mkdir -p "$wd"
  if [[ -n $(ls *${lineage}*/short_summary*.txt) ]]
  then
    cmd="cp -u *${lineage}*/short_summary*.txt $wd"
  elif [[ -n $(ls *${lineage}*/busco_stats/short_summary*.txt) ]]
  then
    cmd="cp -u *${lineage}*/busco_stats/short_summary*.txt $wd"
  else
    echo "Error, $lineage BUSCO summary file not found."
    exit 1
  fi
  echo "$cmd"
  $cmd
  cd "$wd"
  for file in $(ls *.txt)
  do
    if [[ "$lineage" == "quast" ]]
    then
      cmd="mv $file short_summary_${new}.txt"
    else
      cmd="mv $file short_summary.specific.${lineage}.${new}.txt"
    fi
    echo "$cmd"
    #$cmd
  done
  cd ../..
done


# For directory in $sum_dir, generate BUSCO graph
for dir in $(ls "$sum_dir")
do
  wd="${sum_dir}/${dir}"
  cmd="generate_plot.py -wd $wd"
  if [[ "$dir" == "quast" ]]
  then
    conda activate busco_for_quast
  else
    conda activate busco
  fi
  echo "$cmd"
  #$cmd
done
