#!/bin/bash
#SBATCH -J fetch_assemblies
#SBATCH -p oneweek
#SBATCH --mem=5gb
#SBATCH --time=05:00:00
#SBATCH -o %x.log

# Print date and time
date

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]] || (( $# < 2 ))
then
  echo "\
Downloads assembly and associated feature file from given GCA accession(s).
Usage:
  sbatch ${SLURM_JOB_NAME}.sbatch <portal_list> [username] [password]
  sbatch ${SLURM_JOB_NAME}.sbatch <portal_list> [curl_login_file]
Requires: curl"
  exit 0
fi

# Input
scripts_dir_name="s-latissima-genome"
[[ -d $scripts_dir_name ]] && scripts_dir="${scripts_dir_name}/"
# Username and password or curl login file
# Default curl login filename
login="jgi_login"
if [[ "$#" -eq 2 ]]
then
  portal_file="$1"
  login="$2"
  echo "Portal file set to: $portal_file"
  echo "Curl login file set to: $login"
elif [[ "$#" -eq 3 ]]
then
  portal_file="$1"
  user="$2"
  pass="$3"
  echo "Creating curl login file ($login) using given username and password."
  # Load curl module
  module purge
  module load gcc/12.3.0 curl
  # Options for curl login
  log="--data-urlencode login=$user --data-urlencode password=$pass -c $login"
  cmd="curl -L $log https://signon.jgi.doe.gov/signon/create"
  $cmd > /dev/null
else
  echo "Error: too many positional arguments given."
  exit 1
fi
[[ -f "$login" ]] || { echo "Error: curl login file ($login) not found."; exit 1; }
# Portal names (assemblies/portal_names.txt)
if [[ -f "$portal_file" ]]
then
  mapfile -t portal_list < <(cat "$portal_file")
else
  echo "Error: $portal_file not found."
  exit 1
fi

# Script variables
# grep search patterns
ptns=('Repeatmasked\.fasta' 'GeneCatalog\w+\.gff3' 'GeneCatalog.*\.aa\.fasta' 'KOG\.tab')
# Base JGI url
u0="https://genome.jgi.doe.gov"
# XML base URL
uxml="/portal/ext-api/downloads/get-directory?organism="
# Output directory
outdir="assemblies"
mkdir -p "$outdir"
# Output list of species names and list of paths
spc_list="${scripts_dir}species_list.txt"
path_list="${scripts_dir}path_list.txt"
# Output table of species paired with filenames
assembly_file="${scripts_dir}species_table.txt"

# Functions
# Function to parse URL of a file based on pattern in JGI XML file
getUrl () {
  local ptn="$1"
  local xml="$2"
  local url
  # url=$(grep "$ptn" "$xml" | sed 's/.* url="//g' | sed 's/".*//g' | sed 's/amp;//g')
  url="$(grep -P "$ptn" "$xml" | sed 's/.* url="//g' | sed 's/".*//g' | grep -v "&")"
  if [[ -z "$url" ]]
  then
    echo "Error: no URL found for pattern ($ptn) in XML file ($xml)."
    exit 1  
  fi
  url="${u0}${url}"
  echo "$url"
}
# Function to parse MD5 of a file based on pattern in JGI XML file
getMd5 () {
  local ptn="$1"
  local xml="$2"
  local md5
  md5="$(grep -P "$ptn" "$xml" | grep "md5" | sed 's/.* md5="//g' | sed 's/".*//g')"
  if [[ -z "$md5" ]]
  then
    echo "Error: no MD5 found for pattern ($ptn) in XML file ($xml)."
    exit 1  
  fi
  echo "$md5"
}
# Function to download file with secure curl login given URL and output file
getFile () {
  local url="$1"
  local out="$2"
  # Load curl module
  module purge
  module load gcc/12.3.0 curl
  local cmd="curl -L -b $login -o $out $url"
  echo "$cmd"
  $cmd
}
# Function to check if file is compressed, uncompress, and update filename
checkZip () {
  local fname="$1"
  local ext
  ext="$(basename "${fname##*.}")"
  if [[ "$ext" == "gz" ]]
  then
    if [[ -f "$fname" ]]
    then
      # Handles GFF3 tarballs and concatenates into genomic GFF3
      if echo "$fname" | grep -q "\.tar\.gz$"
      then
        mapfile -t gff_list < <(tar -xvf "$fname" | sort -V)
        fname="${fname/\.tar\.gz/}.gff3"
        grep -h "##\w" "${gff_list[@]}" | sort -uV > "${fname}" 
        grep -hv "#" "${gff_list[@]}" >> "${fname}"
        rm "${gff_list[@]}"
      else
        # Load gzip module
        module purge
        module load gcc/11.3.0 gzip
        gunzip -k -f "$fname"
        fname="${fname%%.gz}"
      fi
    else
      echo "Error: file ($fname) not found"
      exit 1
    fi
  fi
  echo "$fname"
}

# Iterate through given portal names
last_i="$((${#portal_list[@]}-1))"
for i in "${!portal_list[@]}"
do
  portal="${portal_list[i]}"
  # Handle explicit links
  if [[ "$portal" == http?(s)://* ]]
  then
    fname="${outdir}/$(basename "$portal")"
    curl -o "$fname" "$portal"
    fname="$(checkZip "$fname")"
    # Get species name from link
    last_sname="$sname"
    sname="$(echo "$portal" | sed 's/.*gdb\/\(.*\)\/.*/\1/g')"
    if [[ -z "$last_sname" ]]
    then
      temp_fnames=()
      temp_fnames+=("$fname")
      continue
    elif [[ "$last_sname" != "$sname" ]]
    then
      spc="$last_sname"
      fnames=("${temp_fnames[@]}")
      # Clear temporary array for next species
      temp_fnames=()
      temp_fnames+=("$fname")
    # Catch final file in array
    elif [[ "$i" -eq "$last_i" ]]
    then
      spc="$sname"
      temp_fnames+=("$fname")
      fnames=("${temp_fnames[@]}")
    elif [[ "$last_sname" == "$sname" ]]
    then
      temp_fnames+=("$fname")
      continue
    fi
  else
    # URL and file path of XML
    url="${u0}${uxml}${portal}"
    xml="${outdir}/${portal}.xml"
    getFile "$url" "$xml"
    # Download files for each accession based on patterns (see "ptns" array)
    fnames=()
    for n in "${!ptns[@]}"
    do
      # Get URL and MD5 from XML file with search pattern
      ptn="${ptns[n]}"
      url="$(getUrl "$ptn" "$xml")"
      md5="$(getMd5 "$ptn" "$xml")"
      # Download file with curl to output file $fname
      fname="${outdir}/$(basename "$url")"
      getFile "$url" "$fname"
      # Check downloaded file against MD5 in XML
      md5_test="$(md5sum "$fname" | awk '{print $1}')"
      echo "XML MD5: $md5"
      echo "$fname MD5: $md5"
      [[ "$md5" != "$md5_test" ]] && { echo "Error: MD5s differ"; exit 1; }
      # Species name
      ptn='<file label="'
      spc="$(grep -P "$ptn" "$xml" | sed "s/$ptn//g" | sed 's/".*//g' | sort -u)"
      # Uncompress downloaded file if needed
      fname="$(checkZip "$fname")"
      fnames+=("$fname")
    done
  fi
  if echo "${fnames[*]}" | grep -iq "ectsiv2"
  then
    spc="Ectocarpus siliculosus Ec 32 V2"
  elif echo "${fnames[*]}" | grep -iq "sacja"
  then
    spc="Saccharina japonica str. Ja"
  fi
  if [[ "$i" -eq 0 ]]
  then
    # printf "%s\t%s\t%s\n" "$spc" "${fnames[0]}" "${fnames[1]}" > "$assembly_file"
    # Save specices names
    echo "$spc"
    echo "$spc" > "$spc_list"
    # Save paths
    echo "${fnames[*]}" | sed 's/ /\t/g'
    echo "${fnames[*]}" | sed 's/ /\t/g' > "$path_list"
  else
    # printf "%s\t%s\t%s\n" "$spc" "${fnames[0]}" "${fnames[1]}" >> "$assembly_file"
    echo "$spc"
    echo "$spc" >> "$spc_list"
    echo "${fnames[*]}" | sed 's/ /\t/g'
    echo "${fnames[*]}" | sed 's/ /\t/g' >> "$path_list"
  fi
done


# Combine species and paths into tab-delimited table
paste "$spc_list" "$path_list" > "$assembly_file"
echo "Find table of species and files in: $assembly_file"

