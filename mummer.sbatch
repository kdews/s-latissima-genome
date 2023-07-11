#!/bin/bash
#SBATCH -J mummer
#SBATCH -n 13
#SBATCH --mem=25gb
#SBATCH --time=02:00:00
#SBATCH -o %x.log
#SBATCH -p htcondor

# Print date and time
date

# Load conda
cond=~/.conda_for_sbatch.sh
[[ -a "$cond" ]] && source "$cond" || \
{ echo "Error on source of $cond"; exit 1; }

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Aligns FASTA file of contigs to genome reference and outputs MUMmer results
(contigs_vs_genome.out), as well as gnuplot of MUMmer alignments 
(contigs_vs_genome.gp, contigs_vs_genome.png, ...).
Usage: sbatch mummer.sbatch <contigs.fasta> <genome.fasta>
Requires:
 - MUMmer (https://github.com/mummer4/mummer)
 - gnuplot (http://www.gnuplot.info)"
  exit 0
fi

# Set variables
ref="$1"
ref="assemblies/Macpyr2_AssemblyScaffolds_Repeatmasked.fasta"
ref_no_ext="$(basename ${ref%.*})"
query="$2"
query="busco_contigs_SugarkelpHifiasm.raconPolished.fasta"
query_no_ext="$(basename ${query%.*})"
prefix="${query_no_ext}_vs_${ref_no_ext}"
q_ids=($(ls ${query_no_ext}*.fa))

# Load gnuplot and GNU parallel
module load gcc/11.3.0 gnuplot/5.4.3 libgd/2.2.4 parallel
# Load MUMmer conda env
conda activate mummer

# Invoke parallel with jobs=$SLURM_NTASKS and a joblog file
thr="$SLURM_CPUS_PER_TASK"
par="parallel -j $SLURM_NTASKS --link"
# srun arguments allocate a single core to the set
# of threads defined by $SLURM_CPUS_PER_TASK
srun="srun --exclusive"
mum_opts="--jl mummer.joblog --res {/.}.out"
cmd="mummer -mum -b -n $ref {}" 
$par $mum_opts $srun $cmd ::: ${q_ids[@]}


# # Don't rerun completed analysis
# if [[ -f "${prefix}.out" ]]
# then
#   echo "Skipping MUMmer alignment, found ${prefix}.out"
# else
#   # Run MUMmer on each split query FASTA
#   q_ids=($(ls ${query_no_ext}*.fa))
#   num=$((${#q_ids[@]}-1))
#   for i in $(seq 0 $num)
#   do
#     contig=${q_ids[i]}
#     cmd="mummer -mum -b -n $ref $contig"
#     echo "$cmd > ${prefix}.out"
#     $cmd > "${prefix}.out"
#   done
# fi

# # Generate MUMmerplot 
# opts="--postscript --large -R $ref"
# for i in $(seq 0 $num)
# do
#   contig=${q_ids[i]}
#   title="$(grep '>' $contig | sed 's/>//g' | sed 's/ .*//g')"
#   echo "$title"
#   cmd="mummerplot $opts -p $title -q $contig -t $title ${prefix}.out"
#   echo "$cmd"
#   $cmd
# done