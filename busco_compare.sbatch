#!/bin/bash
#SBATCH -J busco_compare
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
[[ -z $1 ]] && sum_dir="busco_summaries" || sum_dir="$1"
lineages=("eukaryota_odb10" "stramenopiles_odb10" "quast")


# Compile BUSCO summary text files (if needed)
mkdir -p "$sum_dir"
for lineage in "${lineages[@]}"
do
  wd="${sum_dir}/${lineage}/"
  mkdir -p "$wd"
  # Clear any old text files from working directory
  rm ${wd}/*.txt
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
    new="${file//*short_summary_/}"
    new="${new//*aug_busco_/}"
    # new="${new//Macpyr2/Mac_pyr2}"
    # new="$(echo $new | awk -F_ '{print $1,$2}' | sed 's/ /_/g')"
    if [[ "$lineage" == "quast" ]]
    then
      cmd="mv $file short_summary_${new}"
    else
      cmd="mv $file short_summary.specific.${lineage}.${new}"
    fi
    echo "$cmd"
    $cmd
  done
  cd ../..
done


# For directory in $sum_dir, generate BUSCO graph
for dir in $(ls "$sum_dir")
do
  wd="${sum_dir}/${dir}"
  cmd="generate_plot.py -wd $wd --no_r"
  if [[ "$dir" == "quast" ]]
  then
    conda activate busco_for_quast
  else
    conda activate busco
  fi
  echo "$cmd"
  $cmd
  # Clean up BUSCO graph
  module load gcc/11.3.0 openblas/0.3.20 r/4.2.2
  # module load gcc/11.3.0 openblas/0.3.20 r/4.3.0
  cmd="Rscript s-latissima-genome/pretty_busco.R $wd"
  echo $cmd
  $cmd
done
