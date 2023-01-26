#!/bin/bash
#SBATCH -J busco
#SBATCH -p shared
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=10:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH -o %x.log

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
Usage: sbatch busco.sbatch <assembly.fasta>
Requires: BUSCO (https://busco.ezlab.org)"
  exit 0
fi

# Define variables
mode=genome
L1=eukaryota_odb10
L2=stramenopiles_odb10
assembly=$1
assembly_no_ext=$(basename "${assembly%.*}")
# Output directory name
out1=${assembly_no_ext}_busco_${L1}
out2=${assembly_no_ext}_busco_${L2}
# Shorten SLURM output variables
job="${SLURM_JOB_NAME}"
ntsk="$SLURM_NTASKS"
nthrd="$SLURM_CPUS_PER_TASK"

# Invoke parallel with jobs=$SLURM_NTASKS and a joblog file
# srun gives 1 core to the set of threads defined by $SLURM_CPUS_PER_TASK
parallel="parallel -j $ntsk --joblog ${assembly_no_ext}_${job}.joblog --link"
srun="srun --exclusive -N 1 -n 1 --cpus-per-task=$nthrd"
busco="busco -f -c $nthrd -m $mode"
$parallel "$srun -o {2}.log $busco -l {1} -i $assembly -o {2}" ::: $L1 $L2 ::: $out1 $out2
