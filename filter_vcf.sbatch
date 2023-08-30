#!/bin/bash
#SBATCH -J filter_vcf
#SBATCH -p oneweek
#SBATCH --mem=150gb
#SBATCH -c 4
#SBATCH --time=02:00:00
#SBATCH -o %x.log

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
Filters input VCF using the following specifications:
- GATK4
   - QD < 2.0
   - MQ < 40.0
   - FS > 60.0
   - HaplotypeScore > 13.0
   - MQRankSum < -12.5
   - ReadPosRankSum < -8.0
- VCFtools
   - Remove indels
   - Biallelic SNPs only
   - Quality threshold = 30
   - Site called in ≥90% of individuals
   - Mean depth per site ≥3 reads

Usage: sbatch ${SLURM_JOB_NAME}.sbatch <in.vcf> [reference.fasta]
Requires:
 - GATK (https://gatk.broadinstitute.org/)
 - VCFtools (https://vcftools.github.io/)"
  exit 0
fi

# GATK filtering
# Load GATK conda env
conda activate gatk
gatk --version
# Input data
vcf="$1"
vcf_no_ext="$(basename ${vcf%.*})"
# Check for existence of VCF index file
if [[ -f "${vcf}.tbi" ]]
then
   echo "Using index file: ${vcf}.tbi"
else
   echo "Creating index for VCF file: $vcf"
   cmd="gatk IndexFeatureFile -I $vcf"
   echo "$cmd"
   $cmd
fi
temp_vcf="temp_${vcf_no_ext}.filtered.vcf.gz"
data="-V $vcf -O $temp_vcf"
# Give reference flag if supplied
[[ $# -eq 2 ]] && data="$data -R $2"
gatk VariantFiltration $data \
--filter-name "QD_2" --filter-expression "QD < 2.0" \
--filter-name "MQ_40" --filter-expression "MQ < 40.0" \
--filter-name "FS_60" --filter-expression "FS > 60.0" \
--filter-name "Haplo_13" --filter-expression "HaplotypeScore > 13.0" \
--filter-name "MQRankSum_-12.5" --filter-expression "MQRankSum < -12.5" \
--filter-name "ReadPosRankSum_-8" --filter-expression "ReadPosRankSum < -8.0"

# VCFtools filtering
# Load VCFtools conda env
conda activate vcftools
vcftools --version
out="${vcf_no_ext}.filtered.vcf.gz"
data2="--gzvcf $temp_vcf --out $out --recode --recode-INFO-all"
# Filtering options
opts="--remove-indels --min-alleles 2 --max-alleles 2"
opts2="--minQ 30 --max-missing 0.1 --min-meanDP 3"
cmd="vcftools $data2 $opts $opts2"
echo "$cmd"
$cmd
