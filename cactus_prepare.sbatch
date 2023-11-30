#!/bin/bash
#SBATCH -J cactus_prepare
#SBATCH -p shared
#SBATCH --mem=1gb
#SBATCH --time=00:30:00
#SBATCH -o %x.log

# Print date and time
date

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Prepares commands for Cactus alignment of multiple genomes. Requires a Cactus-
formatted seqFile (Newick tree and species list).
Usage: sbatch ${SLURM_JOB_NAME}.sbatch <seqFile>
Requires:
 - virtualenv (https://virtualenv.pypa.io/en/latest)
 - Toil (https://toil.readthedocs.io/en/latest)
 - Cactus (https://github.com/ComparativeGenomicsToolkit/cactus)"
  exit 0
fi

# Load Python module and Cactus virtualenv
module purge
module load gcc/12.3.0 python/3.11.4 git/2.42.0 py-pip/23.0
env_path="/home1/kdeweese/bin/cactus-bin-v2.6.7/cactus_env/bin/activate"
if [[ -a "$env_path" ]]
then
  source "$env_path"
else
  echo "Error on source of $env_path"
  exit 1
fi

# Input
# Requires Cactus-formatted seqFile, i.e.,
# NEWICK tree
# name1 path1
# name2 path2
# ...
# nameN pathN
# seqFile: (s_lat_alignment.txt)
if [[ "$#" -eq 1 ]]
then
    seqFile="$1"
    if [[ -f "$seqFile" ]]
    then
      echo "Using seqFile: $seqFile"
    else
      echo "Error: no seqFile found ($seqFile)."
      exit 1
    fi
else
  echo "Error: one positional argument expected."
  exit 1
fi
# Other required inputs
js="--jobStore ${SLURM_JOB_NAME}_steps_jobstore"
seqFile_no_ext=$(basename "${seqFile%.*}")
hal="${seqFile_no_ext}.hal"
outdir="--outDir cactus-steps-output"
outseq="--outSeqFile cactus-steps-output/${seqFile}"
outhal="--outHal cactus-steps-output/$hal"
# Memory (in GB) and CPU options
# Preprocessing memory: 20 GB
premem=20
# BLAST memory
# 3 genomes: 35 GB
# 5 genomes: 65 GB
blastmem=65
# Alignment memory
# 3 genomes: 150 GB
# 5 genomes: 
alnmem=300
cpus=14
opts1="--preprocessCores $cpus --blastCores $cpus --alignCores $cpus"
opts2="--preprocessMemory ${premem}G --blastMemory ${blastmem}G --alignMemory ${alnmem}G"
opts="$opts1 $opts2"
cmd_file="cactus_prepare_commands.txt"
sbatch_list_file="cactus_sbatch_list.txt"

# Run cactus-prepare
cmd="cactus-prepare $seqFile $js $outdir $outseq $outhal $opts"
echo "$cmd"
$cmd > $cmd_file

# Split each command into separate SBATCH file in $prep_dir
prep_dir="cactus_prepare_scripts"
mkdir -p "$prep_dir"
mapfile -t commands < <(grep -v '#' $cmd_file | sed '/^[[:space:]]*$/d')
declare -a sbatch_list
for i in "${!commands[@]}"
do
    # Iterate through each command
    cmd_text="${commands[i]}"
    cmd_name="${i}_${cmd_text%% *}"
    cmd_path="${prep_dir}/${cmd_name}.sbatch"
    # Default memory allocation
    mem=20
    # Choose correct memory allocation
    if echo "$cmd_name" | grep -q "preprocess"
    then
      mem="$premem"
    elif echo "$cmd_name" | grep -q "blast"
    then
      mem="$blastmem"
    elif echo "$cmd_name" | grep -q "align"
    then
      mem="$alnmem"
    fi
    # Add 10% more memory for overhead
    reqmem="$((mem * 10 / 9))G"
    # Add 2 cores of overhead
    reqcpus="$((cpus + 2))"
    # Add text to SBATCH file
    echo "Writing $cmd_path"
    {
      echo "#!/bin/bash"
      echo "#SBATCH -J $cmd_name"
      echo "#SBATCH -o %x.log"
      echo "#SBATCH -t 1-0"
      echo "#SBATCH --mem=$reqmem"
      echo "#SBATCH -c $reqcpus"
      echo ""
      echo "date"
      echo ""
      echo "$cmd_text"
    } > "$cmd_path"
    sbatch_list[i]="$cmd_path"
done

# Write ordered list of SBATCH files to text file
printf "%s\n" "${sbatch_list[@]}" > "$sbatch_list_file"
echo "List of SBATCH files for cactus_run_prepared in: $sbatch_list_file"
