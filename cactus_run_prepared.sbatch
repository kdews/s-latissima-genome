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
if [[ -a "$env_path" ]]
then
  source "$env_path"
else
  echo "Error on source of $env_path"
  exit 1
fi

# Define functions
# Function to check given logfile for prior successful run
check_prior () {
    local ok_mess="Successfully deleted the job store"
    # If already run, take old job ID
    if [[ -f "$1" ]] && grep -q "$ok_mess" "$1"
    then
      local id
      id="$(grep "SLURM_JOB_ID" "$logfile" | sed 's/SLURM_JOB_ID = //g')"
      echo "$id"
    else
      echo "no_id"
    fi
}
# Function to write sbatch command
submit () {
    local part
    local cmd
    part="$SLURM_JOB_PARTITION"
    # Choose larger partition for BLAST and alignment jobs
    if echo "$*" | grep -q "blast\|align"
    then 
      part="meb_largemem"
    fi
    cmd="sbatch --parsable -p $part $*"
    echo "$cmd"
}

# Run cactus step-by-step with dependent SLURM jobs
sbatch_list_file="cactus_sbatch_list.txt"
mapfile -t sbatch_list < <(cat "$sbatch_list_file")
for i in "${!sbatch_list[@]}"
do
    sfile="${sbatch_list[i]}"
    echo "Step: $sfile"
    cmd_name="$(basename "${sfile%%.sbatch}")"
    logfile="${cmd_name}.log"
    # Check if job has already run successfully
    oldid="$(check_prior "$logfile")"
    # [[ -z "$oldid" ]] || echo "check_prior found successful run $oldid for $cmd_name"
    # [[ "$oldid" == "no_id" ]] || echo "Successful run of ${cmd_name}: $oldid"
    dep="--dependency=afterok:$jobid"
    # No dependency for first job
    (( i == 0 )) && dep=""
    cmd="$(submit "$dep" "$sfile")"
    # [[ -z "$oldid" ]] && { echo "$cmd"; jobid="$($cmd)"; } || jobid="$oldid"
    if [[ "$oldid" == "no_id" ]]
    then
      echo "$cmd"
      # Submit new job
      if ! jobid="$($cmd)"
      then
        echo "Error submitting command: $cmd"
        exit 1
      fi
    else
      # Use old successful run
      echo "Using previous successful run of ${cmd_name}: $oldid"
      jobid="$oldid"
    fi
    # [[ "$?" -ne 0 ]] && exit "$?"
    echo "$cmd_name job ID: $jobid"
    echo ""
done