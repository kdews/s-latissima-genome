# Clear environment
rm(list = ls())
# Required packages
library(tidyverse, quietly = T, warn.conflicts = F)
library(BiocManager)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
if (require(showtext, quietly = T)) {
  showtext_auto()
  if (interactive()) showtext_opts(dpi = 100) else showtext_opts(dpi = 300)
}

wd <- "/scratch1/kdeweese/latissima/genome_stats"
setwd(wd)
coord_file <- "mummer_0_chromosome_extract_4MB_SlaSLCT1FG3_1_AssemblyScaffolds_Repeatmasked_renamed_by_size_vs_chromosome_extract_Macpyr2_AssemblyScaffolds_Repeatmasked.coords"
coords <- read.table(coord_file, skip = 4, header = F, sep = "\t")
coord_cols <- c("R_Start", "R_End",
                "Q_Start", "Q_End",
                "R_Match", "Q_Match",
                "Perc_ID",
                "R_Len", "Q_Len",
                "R_Cov", "Q_Cov",
                "R_ID", "Q_ID")
colnames(coords) <- coord_cols
coords <- coords %>%
  arrange(desc(R_Len), desc(Q_Len), desc(R_Cov), desc(Q_Cov), desc(Perc_ID))
sum_df <- coords %>%
  group_by(Q_ID, R_ID) %>%
  summarize(n=n())
sum_mat <- as.matrix(sum_df %>%
                       pivot_wider(id_cols = R_ID,
                                   names_from = Q_ID,
                                   values_from = n) %>%
                       column_to_rownames(var = "R_ID"))
heatmap(sum_mat, scale = "column", main = "Number of Alignments")
max_sum <- sum_df %>%
  filter(n == max(n)) %>%
  arrange(desc(n))

len_df <- coords %>%
  group_by(Q_ID, R_ID) %>%
  summarize(Q_Sum=sum(Q_Match))
len_mat <- as.matrix(len_df %>%
                       pivot_wider(id_cols = R_ID,
                                   names_from = Q_ID,
                                   values_from = Q_Sum) %>%
                       column_to_rownames(var = "R_ID"))
heatmap(len_mat, scale = "column", main = "Sum of Alignment Lengths")
max_match <- len_df %>%
  filter(Q_Sum == max(Q_Sum)) %>%
  arrange(desc(Q_Sum))


cov_df <- coords %>%
  group_by(Q_ID, R_ID) %>%
  summarize(Q_Cov=max(Q_Cov))
max_cov <- cov_df %>% 
  filter(Q_Cov == max(Q_Cov)) %>%
  arrange(desc(Q_Cov))
