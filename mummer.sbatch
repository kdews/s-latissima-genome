#!/bin/bash
#SBATCH -J mummer
#SBATCH --mem=10gb
#SBATCH -c 10
#SBATCH --time=01:00:00
#SBATCH -o %x_%a.log
#SBATCH -p htcondor

# Print date and time
date

# Load conda
cond=~/.conda_for_sbatch.sh
if [[ -a "$cond" ]]
then
  source "$cond"
else
  echo "Error on source of $cond"
  exit 1
fi

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Aligns with MUMmer's NUCmer in two modes:
  1. Multi-FASTA query:  Uses entire given query multi-FASTA to align to ref...
  2. Single-FASTA query: Uses specific query contig ID from given query assembly
                         to align to ref, specified by array ID number and list
                         of contigs...
...and generates a dotplot of the alignment.
Usage:
  Multi-FASTA mode
    sbatch -a 0 ${SLURM_JOB_NAME}.sbatch <ref.fasta> <query.fasta>
  Single-FASTA mode
    sbatch -a 1-<int> ${SLURM_JOB_NAME}.sbatch <ref.fasta> <query.fasta> <contig_list.txt>
Requires:
 - seqtk (https://github.com/lh3/seqtk)
 - MUMmer (https://github.com/mummer4/mummer)
 - gnuplot (http://www.gnuplot.info)"
  exit 0
fi

# Set variables
num="$SLURM_ARRAY_TASK_ID"
ref="$1"
assembly="$2"
ref_no_ext="$(basename ${ref%.*})"
assembly_no_ext="$(basename ${assembly%.*})"
if [[ "$num" -eq 0 ]]
then
  query="$assembly"
  query_prefix="${num}_${assembly_no_ext}"
else
  contig_ids="$3"
  contig_id="$(sed -n ${num}p $contig_ids)"
  query_prefix="${num}_${assembly_no_ext}_${contig_id}"
  tempquery="${query_prefix}.fasta"
  query="$tempquery"
  tempfile="mummer_temp_${num}.txt"
fi
mum_prefix="mummer_${query_prefix}_vs_${ref_no_ext}"

# Subset contig
if [[ -z "$contig_id" ]]
then
  echo "No subsetting, running MUMmer with multi-FASTA query: $query"
else
  # Load conda environment for seqtk
  conda activate seqtk
  # Subset list of FASTA IDs into multi-FASTA
  # Create temporary list file for seqtk
  echo "$contig_id" > "$tempfile"
  cmd="seqtk subseq $assembly $tempfile"
  echo "$cmd > $query"
  $cmd > $query
fi

# Alignment of contigs to related reference genome
# Load gnuplot and GNU parallel
conda activate
module purge
module load gcc/11.3.0 gnuplot/5.4.3 libgd/2.2.4
# Load MUMmer conda env
conda activate mummer
# Don't rerun completed analysis
if [[ -f "${mum_prefix}.delta" ]]
then
  echo "Skipping MUMmer alignment, found ${mum_prefix}.delta"
else
  # Run alignment
  cmd="nucmer -p $mum_prefix -t $SLURM_CPUS_PER_TASK -g 1000 $ref $query"
  # cmd="mummer -mum -b -n -l 100 $ref $query"
  echo "$cmd"
  $cmd
fi

# Output coordinates from delta file
cmd="show-coords -Trcl ${mum_prefix}.delta"
echo "$cmd > ${mum_prefix}.coords"
$cmd > "${mum_prefix}.coords"

# Generate MUMmerplot 
opts="-f -l --png --large -R $ref -Q $query"
cmd="mummerplot $opts -p $mum_prefix -t $mum_prefix ${mum_prefix}.delta"
echo "$cmd"
$cmd

# Remove temporary files unless error occured
if [[ $? -eq 0 ]]
then
  [[ -f "$tempfile" ]] && rm "$tempfile"
  [[ -f "$tempquery" ]] && rm "$tempquery"
else
  echo "Error code: $?"
fi
