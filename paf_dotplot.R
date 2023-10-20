# Clear environment
rm(list = ls())
# Required packages
library(pafr, quietly=TRUE)
library(tidyverse, quietly = T, warn.conflicts = F)
library(ggpubr, quietly = T)
# library(gridExtra, quietly = T)
if (require(showtext, quietly = T)) {
  showtext_auto()
  if (interactive()) showtext_opts(dpi = 100) else showtext_opts(dpi = 300)
}

# Input
if (interactive()) {
  wd <- "/scratch2/kdeweese/latissima/genome_stats"
  setwd(wd)
  paf_file <- "s_lat_alignment.paf"
} else {
  line_args <- commandArgs(trailingOnly = T)
  paf_file <- line_args[1]
}

paf <- read_paf(paf_file)
long_paf <- subset(paf, alen > 1e4 & mapq > 40)
dt <- dotplot(long_paf, label_seqs = T, order_by = "qstart")
unique(long_paf$tname)
