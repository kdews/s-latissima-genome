# Clear environment
rm(list = ls())
# Required packages
# library(pafr, quietly=TRUE)
library(tidyverse, quietly = T, warn.conflicts = F)
library(ggpubr, quietly = T)
library(RColorBrewer, quietly = T)
library(ComplexHeatmap, quietly = T)
# library(gridExtra, quietly = T)
if (require(showtext, quietly = T)) {
  showtext_auto()
  if (interactive()) showtext_opts(dpi = 100) else showtext_opts(dpi = 300)
}

# Functions
# Gets genome name from "genome_vs_genome" character
getGen <- function(df_name, gen_type) {
  if (gen_type == "query") {
    n <- 1
  } else if (gen_type == "target") {
    n <- 2
  } else {
    stop("Error: <gen_type> must be either 'query' or 'target'.")
  }
  gen_name <- unlist(strsplit(df_name, "_vs_"))[[n]]
  return(gen_name)
}
# Fixes Undaria chromosome labels for plotting
fixUndaria <- function(contigs) {
  contigs <- gsub("\\..*", "", gsub("JABAKD0100000", "chr_", contigs))
  return(contigs)
}
# Order PSL data frame by contig size by converting contig names to factors
orderPsl <- function(df_name, df_list) {
  df <- df_list[[df_name]]
  query <- getGen(df_name, "query")
  target <- getGen(df_name, "target")
  df <- df %>%
    mutate(qName=factor(qName, levels = contig_ids[[query]][["ID"]]),
           tName=factor(tName, levels = contig_ids[[target]][["ID"]])) %>%
    arrange(qName)
  return(df)
}
# Summarize PSL table by summing matches for each contig pair
sumPsl <- function(df_name, df_list) {
  df <- df_list[[df_name]]
  query <- getGen(df_name, "query")
  target <- getGen(df_name, "target")
  df_sum <- df %>%
    group_by(qName, tName) %>%
    summarize(total_matches=as.numeric(sum(matches)))
  df2 <- merge(df_sum, unique(df[,c("qName", "qSize")]),
               by = "qName", sort = F)
  df2 <- df2 %>%
    mutate(qPercent=total_matches/qSize*100)
  return(df2)
}
# Pull out maximal matching contigs
maxMatches <- function(df_name, df_list) {
  df <- df_list[[df_name]]
  query <- getGen(df_name, "query")
  target <- getGen(df_name, "target")
  df_sum <- df %>% 
    group_by(qName) %>%
    filter(qPercent == max(qPercent)) %>%
    ungroup()
  return(df_sum)
}
# Convert PSL summary into matrix for heatmap
matPsl <- function(df_name, df_list) {
  df <- df_list[[df_name]]
  query <- getGen(df_name, "query")
  target <- getGen(df_name, "target")
  mat <- as.matrix(df %>% mutate(qName=fixUndaria(qName),
                                 tName=fixUndaria(tName)) %>%
                     pivot_wider(id_cols = qName,
                                 names_from = tName,
                                 values_from = qPercent,
                                 values_fill = 0) %>%
                     column_to_rownames(var = "qName"))
                     # select(any_of(contig_ids[[target]][["ID"]])))
  return(mat)
}
# Plot heatmap of given matrix
heatPsl <- function(mat_name, mat_list) {
  mat <- mat_list[[mat_name]]
  query <- getGen(mat_name, "query")
  target <- getGen(mat_name, "target")
  x_label <- paste("target:", gsub("_", " ", target))
  y_label <- paste("query:", gsub("_", " ", query))
  fname <- paste0(mat_name, "_heatmap.png")
  h_colors <- colorRampPalette(brewer.pal(8, "YlOrRd"))
  png(file = fname, units = "in", width = 5, height = 5, res = 150)
  heatmap(mat, main = "Query contig coverage %",
          xlab = x_label, ylab = y_label, col = h_colors(25))
  dev.off()
}
# Subset data for dotplots
dotFilt <- function(df_name, df_list, df2_list) {
  df <- df_list[[df_name]]
  df2 <- df2_list[[df_name]]
  match_list <- paste0(df2$qName, df2$tName)
  target <- getGen(df_name, "target")
  query <- getGen(df_name, "query")
  subset_cols <- c("qName", "tName", "strand",
                   "qSize", "qStart", "qEnd",
                   "tSize", "tStart", "tEnd")
  df_sub <- df[,subset_cols]
  df_sub <- df_sub %>%
    filter(paste0(qName, tName) %in% match_list)
  return(df_sub)
}
# Dotplots
dotPlot <- function(df_name, df) {
# dotPlot <- function(df_name, df_list) {
  # df <- df_list[[df_name]]
  target <- getGen(df_name, "target")
  query <- getGen(df_name, "query")
  fname <- paste0(df_name, "_dotplot.png")
  p <- ggplot(data = df,
              mapping = aes(x = qStart, y = tStart, color = strand)) +
    # geom_point() +
    geom_segment(mapping = aes(xend = qEnd, yend = tEnd),
                 linewidth = 1, lineend = "round") +
    facet_grid(cols = vars(qName), rows = vars(tName),
               switch = "both", space = "free", scale = "free", as.table = F) +
    theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank())
  showtext_opts(dpi = 300)
  ggsave(filename = fname, plot = p, width = 12, height = 6, units = "in")
  showtext_opts(dpi = 100)
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
# Summarize by contig vs. contig of each syntenic comparison
psl_sums <- sapply(names(psl_list), sumPsl, psl_list, simplify = F)
psl_match <- sapply(names(psl_sums), maxMatches, psl_sums, simplify = F)
# Matrices of syntenic blocks by contig for heatmaps
psl_mats <- sapply(names(psl_sums), matPsl, psl_sums, simplify = F)
# Heatmaps of genome vs. genome
lapply(names(psl_mats), heatPsl, psl_mats)
# Subset data for dotplots
psl_dot <- sapply(names(psl_list), dotFilt, psl_list, psl_match, simplify = F)

sub_tnames <- unique(psl_dot[[2]]$qName)[14]
df_sub <- psl_dot[[2]] %>%
  filter(tName %in% sub_tnames) %>%
  # filter(tName %in% "scaffold_1") %>%
  group_by(tName) %>%
  arrange(tStart, qStart) %>%
  mutate(qName = factor(qName, levels = unique(qName))) %>%
  ungroup()

dotPlot("Saccharina_latissima_vs_Macrocystis_pyrifera", df_sub)

scaff_7 <- psl_list[[2]] %>%
  filter(qName == "scaffold_7")
unique(scaff_7[,c("qName", "tName")])
# # Dotplots of genome vs. genome
# sapply(names(psl_dot), dotPlot, psl_dot, simplify = F)
