# Clear environment
rm(list = ls())
# Required packages
library(tidyverse)
library(ape)
library(pegas)
library(seqinr)
library(adegenet)
library(poppr)
if (require(showtext, quietly = TRUE)) {
  showtext_auto()
  if (interactive()) showtext_opts(dpi = 100) else showtext_opts(dpi = 300)
}
# Set working directory
setwd("/scratch1/kdeweese/latissima/genome_stats/")

# Input
vcffile <- "master_SlaSLCT1FG3_1_AssemblyScaffolds_Repeatmasked.vcf.filtered.vcf.gz.recode.vcf"
vcf <- read.vcf(vcffile, to = 1729853)
class(vcf)
