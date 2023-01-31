#!/bin/bash
#SBATCH --mem=15gb
#SBATCH --time=01:00:00

prefix=$1
qc_dir=$2
# Optional, to specify directory containing scripts
scripts_dir=$3

# Load conda
source $(which activate) || { echo "Error sourcing $cond"; exit 1; }
# Load multiqc
conda activate multiqc
multiqc --version
if [[ $? -ne 0 ]]
then
	echo "Error - check MultiQC installation."
fi

# If input directory exists, run MultiQC
if [[ -d $qc_dir ]]
then
	multiqc -c ${scripts_dir}multiqc_config.yaml \
-o $qc_dir $qc_dir
else
	echo "Error - input directory to MultiQC, \
${qc_dir}, not found."
fi

# Checkpoints
if [[ $? -eq 0 ]]
then
	if [[ -f ${qc_dir}/multiqc_report.html  ]]
	then
		if [[ ! -d checkpoints ]]
		then
			mkdir checkpoints
		fi
		touch checkpoints/${prefix}.checkpoint
		exit 0
	else
		echo "Error - MultiQC report not created."
		exit 1
	fi
else
	echo "Error - non-zero exit code returned by MultiQC."
	exit 1
fi

