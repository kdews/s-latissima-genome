#!/bin/bash
#SBATCH -J fetch_assemblies
#SBATCH -p oneweek
#SBATCH --mem=5gb
#SBATCH --time=05:00:00
#SBATCH -o %x.log

# Print date and time
date

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Downloads assembly and associated feature file from given GCA accession(s).
Usage: sbatch fetch_assemblies.sbatch <portal_list>
Requires: NCBI Datasets CLI (https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/)"
  exit 0
fi

# Load conda
cond=~/.conda_for_sbatch.sh
if [[ -a "$cond" ]]
then
  source "$cond"
else
  echo "Error on source of $cond"
  exit 1
fi
# # Load unzip module
# module load gcc/11.3.0 gcc/12.3.0 gcc/8.3.0 intel/19.0.4 unzip

# Input accession(s)
if [[ -z "$1" ]]
then
  portal_file="assemblies/portal_list.txt"
else
  portal_file="$1"
fi
if [[ -f "$portal_file" ]]
then
  portal_list=( $(cat "$portal_file") )
else
  echo "Error: $portal_file not found."
  exit 1
fi
# Output directory
outdir="assemblies"
mkdir -p "$outdir"

# Download both genome and GFF3 feature file for each accession
for acc in "${portal_list[@]}"
do
  curl 'https://signon.jgi.doe.gov/signon/create' --data-urlencode 'login=kdeweese@usc.edu' --data-urlencode 'password=That1214!' -c cookies > /dev/null
  curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory?organism=Ectsil1' -b cookies > Ectsil1.xml
  fname="${outdir}/${acc}.zip"
  cmd="datasets download genome accession $acc $opts --filename $fname"
  echo "$cmd"
  $cmd
  newdir="${outdir}/${acc}"
  cmd="mkdir -p $newdir"
  echo "$cmd"
  $cmd
  cmd="unzip $fname -d $newdir"
  echo "$cmd"
  $cmd
  # cmd="mv ${newdir}/ncbi_dataset/data/${acc}/* $outdir"
  # echo "$cmd"
  # $cmd
  # cmd="rm -r $newdir"
  # echo "$cmd"
  # $cmd
done
