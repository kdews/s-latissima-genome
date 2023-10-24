# Clear environment
rm(list = ls())
# Required packages
# library(pafr, quietly=TRUE)
library(tidyverse, quietly = T, warn.conflicts = F)
library(ggpubr, quietly = T)
library(RColorBrewer, quietly = T)
# library(gridExtra, quietly = T)
if (require(showtext, quietly = T)) {
  showtext_auto()
  if (interactive()) showtext_opts(dpi = 100) else showtext_opts(dpi = 300)
}

# Functions
# Convert contig IDs to factors and use to order data frame
orderPsl <- function(df_name, df_list) {
  df <- df_list[[df_name]]
  query <- unlist(strsplit(df_name, "_vs_"))[[1]]
  target <- unlist(strsplit(df_name, "_vs_"))[[2]]
  df <- df %>%
    mutate(qName=factor(qName, levels = contig_ids[[query]][["ID"]]),
           tName=factor(tName, levels = contig_ids[[target]][["ID"]])) %>%
    arrange(qName)
  return(df)
}
# Summarize PSL table by summing matches for each contig pair
sumPsl <- function(df_name, df_list) {
  df <- df_list[[df_name]]
  query <- unlist(strsplit(df_name, "_vs_"))[[1]]
  target <- unlist(strsplit(df_name, "_vs_"))[[2]]
  df_sum <- df %>%
    group_by(qName, tName) %>%
    summarize(total_matches=as.numeric(sum(matches)))
  # target_lens <- contig_ids[[target]][,c("ID", "Length")]
  df2 <- merge(df_sum, unique(df[,c("qName", "qSize")]),
               by = "qName", sort = F)
  df2 <- df2 %>%
    mutate(tPercent=total_matches/qSize*100)
  return(df2)
}
# Convert PSL summary into matrix for heatmap
matPsl <- function(df_name, df_list) {
  df <- df_list[[df_name]]
  query <- unlist(strsplit(df_name, "_vs_"))[[1]]
  target <- unlist(strsplit(df_name, "_vs_"))[[2]]
  mat <- as.matrix(df %>%
                     pivot_wider(id_cols = qName,
                                 names_from = tName,
                                 values_from = tPercent,
                                 values_fill = 0) %>%
                     column_to_rownames(var = "qName") %>%
                     select(any_of(contig_ids[[target]][["ID"]])))
  return(mat)
}

# Input
# Only take command line input if not running interactively
if (interactive()) {
  wd <- "/scratch2/kdeweese/latissima/genome_stats"
  setwd(wd)
  # paf_file <- "s_lat_alignment.paf"
  species_file <- "species.txt"
} else {
  line_args <- commandArgs(trailingOnly = T)
  # paf_file <- line_args[1]
  species_file <- line_args[1]
}
# # PAF analysis
# paf <- read_paf(paf_file)
# long_paf <- subset(paf, alen > 1e4 & mapq > 40)
# dt <- dotplot(long_paf, label_seqs = T, order_by = "qstart")
# Set column names for PSL data frames
psl_col <- c("matches", "misMatches", "repMatches", "nCount", "qNumInsert",
             "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName",
             "qSize", "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd",
             "blockCount", "blockSizes", "qStarts", "tStarts")
faidx_col <- c("ID", "Length", "Offset", "Linebases", "Linewidth")
# Import species names, contig IDs and lengths
species_tab <- read.table(species_file)
species <- species_tab$V1
sp_prefix <- basename(tools::file_path_sans_ext(species_tab$V2))
contig_files <- paste0("chromosome_extract_", sp_prefix, ".fasta.fai")
contig_ids <- sapply(contig_files, read.table, simplify = F,
                     col.names = faidx_col)
names(contig_ids) <- species
# Import PSL files into data frames
# Keep only unique species combinations
species_versus <- lapply(combn(rev(species), 2, simplify = F),
                         paste, collapse = "_vs_")
psl_files <- grep(paste(species_versus, collapse = "|"), 
                  list.files(pattern = "\\.psl"), value = T)
psl_list <- sapply(psl_files, read.table, col.names = psl_col,
                   simplify = F, USE.NAMES = T)
names(psl_list) <- gsub(".*ment_|.psl", "", names(psl_list))
psl_list <- sapply(names(psl_list), orderPsl, psl_list, simplify = F)

# Analysis
psl_sums <- sapply(names(psl_list), sumPsl, psl_list, simplify = F)
psl_mats <- sapply(names(psl_sums), matPsl, psl_sums, simplify = F)

# Plotting
heatmap(psl_mats[[1]])
heatmap(psl_mats[[2]])
heatmap(psl_mats[[3]])

subset_cols <- c("qName", "tName", "strand",
                 "qSize", "qStart", "qEnd",
                 "tSize", "tStart", "tEnd")
test <- psl_list[[1]][,subset_cols]
colnames(test)
head(test)





