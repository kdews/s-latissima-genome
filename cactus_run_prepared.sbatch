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

# Define functions
# Function to check given logfile for prior successful run
check_prior () {
    local ok_mess="Successfully deleted the job store"
    # If already run, take old job ID
    if [[ -f "$1" ]] && [[ $(grep "$ok_mess" "$1") ]]
    then
        local id=$(grep "SLURM_JOB_ID" "$logfile" | sed 's/SLURM_JOB_ID = //g')
        echo "$id"
    fi
}
# Function to write sbatch command
submit () {
    part="$SLURM_JOB_PARTITION"
    [[ $(echo $* | grep "align") ]] && part="meb_largemem"
    cmd="sbatch --parsable -p $part $*"
    echo "$cmd"
}

# Run cactus step-by-step with dependent SLURM jobs
prep_dir="cactus_prepare_scripts"
num=$(ls "$prep_dir" | wc -l)
for i in $(seq 1 "$num")
do
    sfile=$(ls "$prep_dir" | sed -n ${i}p)
    echo "$sfile"
    cmd_name="${sfile%%.sbatch}"
    logfile="${cmd_name}.log"
    sfile="${prep_dir}/${sfile}"
    # Check if job has already run successfully
    oldid="$(check_prior $logfile)"
    [[ -z "$oldid" ]] || echo "check_prior found successful run $oldid for $cmd_name"
    dep="--dependency=afterok:$jobid"
    # No dependency for first job
    (( i == 1 )) && dep=""
    cmd="$(submit $dep $sfile)"
    [[ -z "$oldid" ]] && { echo "$cmd"; jobid="$($cmd)"; } || jobid="$oldid"
    [[ $? -ne 0 ]] && exit $?
    echo "$cmd_name job ID: $jobid"
    echo ""
done