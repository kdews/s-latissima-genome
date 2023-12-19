#!/bin/bash
#SBATCH -J sdr_analysis
#SBATCH -p oneweek
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=05:00:00
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=12
#SBATCH -o %x.log

# Print date and time
date
echo

# Help message
if [ $1 = "-h" ] || [ $1 = "--help" ]
then
  echo "\
Takes as input protein FASTAs intended for OrthoFinder.
Usage: sbatch ${SLURM_JOB_NAME}.sbatch </path/to/ortho_dirs>

Requires:
 - OrthoFinder (https://github.com/davidemms/OrthoFinder)"
  exit 0
fi

# Load GNU Parallel
module purge
module load gcc/11.3.0 parallel 

# Change directory to where directories containing protein
# FASTAs are located named "ortho_etc/"
path_to_ortho_dirs=$1
cd $path_to_ortho_dirs

# Use GNU parallel to parallelize jobs
pll="parallel -j $SLURM_NTASKS --joblog ${SLURM_JOB_NAME}.joblog"
# Allocate single core to the set of threads defined by $SLURM_CPUS_PER_TASK
sr="srun --exclusive -N 1 -n $SLURM_CPUS_PER_TASK"
cmd="$sr orthofinder -t $SLURM_CPUS_PER_TASK -f {}"
$pll "echo $cmd && $cmd" ::: $(ls -d ortho*/)

# Run OrthoFinder in directories starting with "ortho" in path_to_ortho_dirs
# parallel -j "$SLURM_NTASKS" --joblog orthofinder.joblog \
# srun --exclusive -N 1 -n $SLURM_CPUS_PER_TASK \
# orthofinder -t $SLURM_CPUS_PER_TASK -f {} ::: \
# `ls -d ortho*/`
