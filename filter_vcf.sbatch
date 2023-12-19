#!/bin/bash
#SBATCH -J filter_vcf
#SBATCH -p htcondor
#SBATCH --mem=5gb
#SBATCH --time=02:00:00
#SBATCH -o %x.log

# Initialize script
init="s-latissima-genome/script_init.sh"
if [[ -a "$init" ]]; then source "$init"; else { echo "Init err"; exit 1; }; fi

# Help message
if [[ $1 == "-h" ]] || [[ $1 == "--help" ]]
then
  echo "\
Filters input VCF using the following specifications:
- GATK4
   - QD < 2.0
   - MQ < 40.0
   - FS > 60.0
   - MQRankSum < -12.5
   - ReadPosRankSum < -8.0
- VCFtools
   - Remove indels
   - Biallelic SNPs only
   - Quality threshold = 30
   - Site called in ≥90% of individuals
   - Mean depth per site ≥3 reads

Usage: sbatch ${job_name}.sbatch <in.vcf> [reference.fasta]
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
vcf_no_ext="$(basename "${vcf%.*}")"
# Check for existence of VCF index file, and index if necessary
if [[ -f "${vcf}.tbi" ]]
then
   echo "Using index file: ${vcf}.tbi"
else
   echo "Creating index for VCF file: $vcf"
   gatk IndexFeatureFile -I "$vcf"
fi
temp1="temp1_${vcf_no_ext}.filtered.vcf.gz"
temp2="temp2_${vcf_no_ext}.filtered.vcf.gz"
data="-V $vcf -O $temp1"
# Give reference flag if supplied
[[ $# -eq 2 ]] && data="$data -R $2"
# Mark variants that don't pass given filters
echo "Marking variants that don't pass filters."
gatk VariantFiltration --verbosity ERROR "$data" \
--filter-name "QD_2" --filter-expression "QD < 2.0" \
--filter-name "MQ_40" --filter-expression "MQ < 40.0" \
--filter-name "FS_60" --filter-expression "FS > 60.0" \
--filter-name "MQRankSum_-12.5" --filter-expression "MQRankSum < -12.5" \
--filter-name "ReadPosRankSum_-8" --filter-expression "ReadPosRankSum < -8.0"
# Remove variants that don't pass fitlers (FILTER ≠ PASS)
echo "Applying filters."
gatk FilterVcf -I "$temp1" -O "$temp2"

# VCFtools filtering
# Load VCFtools conda env
conda activate vcftools
vcftools --version
out="${vcf_no_ext}.filtered.vcf.gz"
data2="--gzvcf $temp2 --out $out --recode --recode-INFO-all"
# Filtering options
opts="--remove-indels --min-alleles 2 --max-alleles 2"
opts2="--minQ 30 --max-missing 0.1 --min-meanDP 3"
cmd="vcftools $data2 $opts $opts2"
echo "$cmd"
$cmd

if [[ $? -eq 0 ]]
then
   echo "Successfully filtered $vcf"
   echo "Removing temporary files (${temp1}, ${temp2})"
   echo "Removing temporary file indices (${temp1}.tbi, ${temp2}.tbi)"
   rm "$temp1" "$temp2"
   echo "Find results in $out"
fi