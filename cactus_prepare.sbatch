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
Prepares commands for Cactus alignment of multiple genomes. Can create a Cactus-
formatted seqFile from a Newick tree and species list, if provided.
Usage:
    no seqFile mode: sbatch ${SLURM_JOB_NAME}.sbatch <species> <tree>
    seqFile mode: sbatch ${SLURM_JOB_NAME}.sbatch <seqFile>

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
[[ -a "$env_path" ]] && source "$env_path" || \
{ echo "Error on source of $env_path"; exit 1; }

# Input file(s)
# Default seqFile
seqFile="s_lat_alignment.txt"
# If needed, generates formatted seqFile, i.e.,
# NEWICK tree
#
# name1 path1
# name2 path2
# ...
# nameN pathN
if (( $# == 1 ))
then
    seqFile="$1"
    [[ -f $seqFile ]] && echo "Found seqFile: $seqFile"
elif (( $# == 2 ))
then
    species="$1"
    tree="$2"
    echo "Species list: $species."
    echo "Newick tree: $tree"
    [[ -f $species ]] || { echo "Error: Species list not found."; exit 1; }
    [[ -f $tree ]] || { echo "Error: Tree not found."; exit 1; }
    echo "Generating formatted seqFile ($seqFile)..."
    cp "$tree" "$seqFile"
    echo "" >> "$seqFile"
    cat "$species" >> "$seqFile"
elif [[ -f $seqFile ]]
then
    echo "Using default seqFile: $seqFile"
else
  echo "Error: no seqFile ($seqFile), tree ($tree) or species list ($species)."
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
reqmem=35
reqcpus=16
# Leave 20% of memory for overhead
mem="$((reqmem * 4 / 5))G"
# Leave 2 cores of overhead
cpus="$((reqcpus - 2))"
opts1="--defaultCores $cpus --preprocessCores $cpus --blastCores $cpus --alignCores $cpus"
opts2="--defaultMemory $mem --preprocessMemory $mem --blastMemory $mem --alignMemory $mem"
opts="$opts1 $opts2"
cmd_file="cactus_prepare_commands.txt"

# Run cactus-prepare
cmd="cactus-prepare $seqFile $js $outdir $outseq $outhal $opts"
echo "$cmd"
$cmd > $cmd_file

# Split each command into separate SBATCH file in $prep_dir
prep_dir="cactus_prepare_scripts"
mkdir -p "$prep_dir"
commands="$(grep -v '#' $cmd_file | sed '/^[[:space:]]*$/d')"
num=$(echo "$commands" | wc -l)
for i in $(seq 1 "$num")
do
    # Iterate through each command
    cmd_text=$(echo "$commands" | sed -n ${i}p)
    cmd_name="${i}_${cmd_text%% *}"
    cmd_path="${prep_dir}/${cmd_name}.sbatch"
    echo "Writing $cmd_path"
    # Add text to SBATCH file
    echo "#!/bin/bash" > "$cmd_path"
    echo "#SBATCH -J $cmd_name" >> "$cmd_path"
    echo "#SBATCH -o %x.log" >> "$cmd_path"
    echo "#SBATCH -t 1-0" >> "$cmd_path"
    echo "#SBATCH --mem=${reqmem}G" >> "$cmd_path"
    echo "#SBATCH -c $reqcpus" >> "$cmd_path"
    echo " " >> "$cmd_path"
    echo "date" >> "$cmd_path"
    echo "$cmd_text" >> "$cmd_path"
done
