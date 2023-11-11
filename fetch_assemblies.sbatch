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

# Load curl module
module purge
module load gcc/12.3.0 curl/8.1.2

# Input portal(s)
portal_file="assemblies/portal_names.txt"
if [[ -f "$portal_file" ]]
then
  portal_list=( $(cat "$portal_file") )
else
  echo "Error: $portal_file not found."
  exit 1
fi
# Input username and password
if [[ "$#" -eq 0 ]]
then
  echo "Error: username and password requred."
  exit 1
else
  user="$1"
  pass="$2"
fi
# JGI login
login="jgi_login"
user="$1"
pass="$2"
opts="-L --data-urlencode login=$user --data-urlencode password=$pass -c $login"
if [[ -f "$login" ]]
then
  echo "Using login file: $login"
else
  cmd="curl $opts https://signon.jgi.doe.gov/signon/create"
  echo "$cmd > /dev/null"
  $cmd > /dev/null
fi
# Base JGI url
u1="https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory?organism="
# Output directory
outdir="assemblies"
mkdir -p "$outdir"

# Download both genome and GFF3 feature file for each accession
for portal in "${portal_list[@]}"
do
  # XML file
  xml="${outdir}/${portal}.xml"
  cmd="curl -L -b $login -o $xml ${u1}${portal}"
  echo "$cmd"
  $cmd
  # Masked genome file
  ptn="Repeatmasked\.fasta.*fileType=\"Assembly\""
  u2=$(grep "$ptn" "$xml" | sed 's/.* url="//g' | sed 's/".*//g' | sed 's/amp;//g')
  u2="https://genome.jgi.doe.gov${u2}"
  out="${outdir}"/$(basename "$u2")
  cmd="curl -L -b $login -o $out $u2"
  echo "$cmd"
  $cmd
  # Feature file
  ptn="GeneCatalog.*gff3.*fileType=\"Annotation\""
  u3=$(grep "$ptn" "$xml" | sed 's/.* url="//g' | sed 's/".*//g' | sed 's/amp;//g')
  u3="https://genome.jgi.doe.gov${u3}"
  out="${outdir}"/$(basename "$u3")
  cmd="curl -L -b $login -o $out $u3"
  echo "$cmd"
  $cmd
done

