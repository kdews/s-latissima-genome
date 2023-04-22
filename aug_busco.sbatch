#!/bin/bash
#SBATCH -J aug_busco
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=10:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH -o %x.log

# Print date and time
date

# Load conda and BUSCO env
cond=~/.conda_for_sbatch.sh
[[ -a "$cond" ]] && source "$cond" || \
{ echo "Error on source of $cond"; exit 1; }
conda activate busco
# Load GNU Parallel
module load parallel

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Runs BUSCO on given assembly using eukaryota_odb10 and stramenopiles_odb10
lineages.
Usage: sbatch busco_augustus.sbatch <assembly.fasta>
Requires: BUSCO (https://busco.ezlab.org)"
  exit 0
fi


# Define variables
mode=genome
L1=eukaryota_odb10
L2=stramenopiles_odb10
assembly=$1
assembly_no_ext=$(basename "${assembly%.*}")
# Shorten SLURM output variables
job="$SLURM_JOB_NAME"
ntsk="$SLURM_NTASKS"
nthrd="$SLURM_CPUS_PER_TASK"
# Output name
out="${job}_${assembly_no_ext}"


# Invoke parallel with jobs=$SLURM_NTASKS and a joblog file
# srun gives 1 core to the set of threads defined by $SLURM_CPUS_PER_TASK
parallel="parallel -j $ntsk --joblog ${out}.joblog --link"
srun="srun --exclusive -N 1 -n 1 --cpus-per-task=$nthrd"
busco="busco -f -c $nthrd -m $mode --augustus"
$parallel "$srun $busco -l {} -i $assembly -o ${out}_{}" ::: $L1 $L2
