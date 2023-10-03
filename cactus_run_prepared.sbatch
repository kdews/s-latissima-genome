#!/bin/bash
#SBATCH -J cactus_run_prepared
#SBATCH -p shared
#SBATCH --mem=50mb
#SBATCH --time=00:30:00
#SBATCH -o %x.log

# Print date and time
date

# Load Python module and Cactus virtualenv
module purge
module load gcc/12.3.0 python/3.11.4 git/2.42.0 py-pip/23.0
env_path="/home1/kdeweese/bin/cactus-bin-v2.6.7/cactus_env/bin/activate"
[[ -a "$env_path" ]] && source "$env_path" || \
{ echo "Error on source of $env_path"; exit 1; }

# Run cactus step-by-step with dependent SLURM jobs
prep_dir="cactus_prepare_scripts"
num=$(ls "$prep_dir" | wc -l)
sbatch="sbatch --parsable -p ${SLURM_JOB_PARTITION}"
#mem_opts="-t 1-0 --mem=20gb -c 10"
for i in $(seq 1 "$num")
do
    sfile=$(ls "$prep_dir" | sed -n ${i}p)
    sfile="${prep_dir}/${sfile}"
    if (( $i == 1 ))
    then
        cmd="$sbatch $sfile"
        echo "$cmd"
        jobid=$($cmd)
        echo "Submitted $sfile with Slurm job ID: $jobid"
    else
        cmd="$sbatch --dependency=afterok:$jobid $sfile"
        echo "$cmd"
        jobid=$($cmd)
        echo "Submitted $sfile with Slurm job ID: $jobid"
    fi
done