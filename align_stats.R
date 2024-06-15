## Initialization
rm(list = ls())
# Load required packages
suppressPackageStartupMessages(library(tidyverse, quietly = T, warn.conflicts = F))




max_match_file <- "max_matches.tsv"
df <- read.table(max_match_file, sep = "\t", header = T)

uniq_n <- length(unique(df$qNum))
uniq_len <- sum(unique(df[,c("qNum", "qSize")])[,"qSize"])*1e-6
uniq_len_perc <- (uniq_len/615)*100
uniq_match <- sum(unique(df[,c("qNum", "total_matches")])[,"total_matches"])*1e-6
uniq_match_perc <- (uniq_match/615)*100
146.5/615
