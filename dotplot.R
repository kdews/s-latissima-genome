# Clear environment
rm(list = ls())
# Required packages
# library(pafr, quietly=TRUE)
library(tidyverse, quietly = T, warn.conflicts = F)
library(ggpubr, quietly = T)
# library(gridExtra, quietly = T)
if (require(showtext, quietly = T)) {
  showtext_auto()
  if (interactive()) showtext_opts(dpi = 100) else showtext_opts(dpi = 300)
}

# Input
psl_col <- c("matches", "misMatches", "repMatches", "nCount", "qNumInsert",
             "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName",
             "qSize", "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd",
             "blockCount", "blockSizes", "qStarts", "tStarts")
if (interactive()) {
  wd <- "/scratch2/kdeweese/latissima/genome_stats"
  setwd(wd)
  # paf_file <- "s_lat_alignment.paf"
  psl_files <- list.files(pattern = ".psl")
} else {
  line_args <- commandArgs(trailingOnly = T)
  psl_files <- line_args
  # paf_file <- line_args[1]
}

# Functions
# Reads in PSL format files with appropraite column names
read_psl <- function(fname, psl_col) {
  df <- read.table(fname)
  colnames(df) <- psl_col
  return(df)
}


# Analysis
# paf <- read_paf(paf_file)
# long_paf <- subset(paf, alen > 1e4 & mapq > 40)
# dt <- dotplot(long_paf, label_seqs = T, order_by = "qstart")
psl_list <- sapply(psl_files, read_psl, psl_col,
                   simplify = F, USE.NAMES = T)
names(psl_list) <- gsub(".*ment_|.psl", "", names(psl_list))
psl_1 <- psl_list[["Macrocystis_pyrifera_vs_Saccharina_latissima"]]
unique(psl_1$qName)
unique(psl_1$tName)

heatmap(psl_1$matches)

# Plotting
ggplot(data = psl_1, mapping = aes(x = blockCount, y = matches)) +
  geom_point(aes(color = qName))






