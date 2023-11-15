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

# Input
# Portal names
portal_file="assemblies/portal_names.txt"
if [[ -f "$portal_file" ]]
then
  portal_list=( $(cat "$portal_file") )
else
  echo "Error: $portal_file not found."
  exit 1
fi
# Username and password or curl login file
# Default curl login filename
login="jgi_login"
if [[ "$#" -eq 0 ]]
then
  echo "Error: provide either username and password or curl login file."
  exit 1
elif [[ "$#" -eq 1 ]]
then
  echo "Curl login file set to: $1"
  login="$1"
  if ! [[ -f "$login" ]]
  then
    echo "Error: curl login file ($login) not found."
    exit 1
  fi
elif [[ "$#" -eq 2 ]]
then
  echo "Creating curl login file ($login) using given username and password."
  user="$1"
  pass="$2"
  # Options for curl login
  log="--data-urlencode login=$user --data-urlencode password=$pass -c $login"
  cmd="curl -L $log https://signon.jgi.doe.gov/signon/create"
  $cmd > /dev/null
else
  echo "Error: too many positional arguments given."
  exit 1
fi
# Base JGI url
u0="https://genome.jgi.doe.gov"
# XML base URL
uxml="/portal/ext-api/downloads/get-directory?organism="
# Output directory
outdir="assemblies"
mkdir -p "$outdir"
# Output table of species paired with filenames
outtab="${outdir}/species_table.txt"

# Functions
# Function to parse URL of a file based on pattern in JGI XML file
getUrl () {
  local ptn="$1"
  local xml="$2"
  local url
  url=$(grep "$ptn" "$xml" | sed 's/.* url="//g' | sed 's/".*//g' | sed 's/amp;//g')
  if [[ -z "$url" ]]
  then
    echo "Error: no URL found for pattern ($ptn) in XML file ($xml)."
    exit 1  
  fi
  url="${u0}${url}"
  echo "$url"
}
# Function to download file with secure curl login given URL and output file
getFile () {
  local url="$1"
  local out="$2"
  local cmd="curl -L -b $login -o $out $url"
  echo "$cmd"
  $cmd
}
# Function to check if file is zipped, and unzip if it is, and update filename
checkZip () {
  local fname="$1"
  local ext
  ext="$(basename "${fname##*.}")"
  if [[ "$ext" == "gz" ]]
  then
    if [[ -f "$fname" ]]
    then
      gunzip -f "$fname"
      fname="${fname%%\.gz}"
    else
      echo "Error: file ($fname) not found"
      exit 1
    fi
  fi
  echo "$fname"
}

# Download both genome and GFF3 feature file for each accession
for i in "${!portal_list[@]}"
do
  portal="${portal_list[i]}"
  # URL and file path of XML
  url="${u0}${uxml}${portal}"
  xml="${outdir}/${portal}.xml"
  getFile "$url" "$xml"
  # Masked genome file
  ptn='Repeatmasked\.fasta.*fileType="Assembly"'
  url=$(getUrl "$ptn" "$xml")
  gen="${outdir}"/$(basename "$url")
  getFile "$url" "$gen"
  gen="$(checkZip "$gen")"
  # Feature file
  ptn='GeneCatalog.*gff3.*fileType="Annotation"'
  url=$(getUrl "$ptn" "$xml")
  annot="${outdir}"/$(basename "$url")
  getFile "$url" "$annot"
  annot="$(checkZip "$annot")"
  # Species name
  ptn='<file label="'
  spc=$(grep "$ptn" "$xml" | sed "s/$ptn//g" | sed 's/".*//g' | sort -u)
  # Save specices name to tab-delimited table with genome and annotation paths
  if [[ "$i" -eq 0 ]]
  then
    printf "%s\t%s\t%s\n" "$spc" "$gen" "$annot" > "$outtab"
  else
    printf "%s\t%s\t%s\n" "$spc" "$gen" "$annot" >> "$outtab"
  fi
done

echo "Find table of species and files in: $outtab"
