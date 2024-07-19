## Initialization
rm(list = ls())
# Load required packages
library(scales, quietly = T, warn.conflicts = F)
library(cooccur)
library(visNetwork)
suppressPackageStartupMessages(library(tidyverse, quietly = T, warn.conflicts = F))
# suppressPackageStartupMessages(library(ComplexHeatmap, quietly = T))


# Input
# Only take command line input if not running interactively
if (interactive()) {
  wd <- "/project/noujdine_61/kdeweese/latissima/genome_stats"
  setwd(wd)
  # Cactus seqFile
  seq_file <- "s-latissima-genome/s_lat_alignment.txt"
  # Species of interest
  spc_int <- "Saccharina_latissima"
  # Output directory
  outdir <- "s-latissima-genome/"
} else {
  line_args <- commandArgs(trailingOnly = T)
  seq_file <- line_args[1]
  spc_int <- line_args[2]
  outdir <- line_args[3]
}
match_sums_file <- "align_sums.tsv"
max_match_file <- "max_matches.tsv"
# Output files

# Prepend output directory to file name (if it exists)
if (dir.exists(outdir)) {
  match_sums_file <- paste0(outdir, match_sums_file)
}


# Functions
# Order species by decreasing relatedness
spc_order <- c("latissima","japonica", "pyrifera", "pinnatifida", "siliculosus")
# Abbreviate species genus name
abbrevSpc <- function(spc) {
  spc <- unlist(strsplit(spc, "_| "))
  let1 <- substr(spc[1], 1, 1)
  spc[1] <- paste0(let1, ".")
  spc_a <- paste(spc, collapse = " ")
  return(spc_a)
}
# Fixes chromosome labels for plotting
fixChrom <- function(contigs) {
  contigs <-
    as.character(as.numeric(str_remove_all(str_remove_all(contigs, ".*_"),
                                           "[^0-9]")))
  return(contigs)
}
# Order given ID column by size column in data frame
sizeSort <- function(df, id_col, size_col) {
  x <- df %>% select(all_of(c(id_col, size_col))) %>%
    arrange(desc(pick(all_of(size_col)))) %>%
    pull(all_of(id_col)) %>%
    unique()
  return(x)
}
# Order PSL data frame by contig size by converting contig names to factors
orderPsl <- function(df) {
  # df <- df %>%
  #   mutate(tName=paste(target, tName, sep = ";")) %>%
  #   mutate(qName=paste(query, qName, sep = ";"))
  # Order ID character vectors by size
  # query_contigs <- sizeSort(df, "qName", "qSize")
  # target_contigs <- sizeSort(df, "tName", "tSize")
  # query_nums <- sizeSort(df, "qNum", "qSize")
  # target_nums <- sizeSort(df, "tNum", "tSize")
  # Capture size rank of each target ID with index
  target_idx <- df %>% select(tName, tSize, qSize) %>%
    arrange(desc(tSize), desc(qSize)) %>% unique %>%
    # # Filter out artificial chromosomes for ranking
    # filter(fixChrom(tName) != "0") %>%
    # Use row ids for index
    rowid_to_column(var="index") %>%
    select(tName, index) %>%
    # # Force any artificial chromosomes to 0 index
    # rbind(data.frame(index = 0, tNum = "0")) %>%
    # Convert to named vector
    deframe()
  df <- df %>%
    # arrange(target, desc(tSize), desc(qSize)) %>%
    arrange(desc(tSize), desc(qSize)) %>%
    # Convert all IDs to factors ordered by size
    mutate(
      qName=factor(qName, levels = unique(qName)),
      tName=factor(tName, levels = unique(tName)),
      # qNum=factor(qNum, levels = query_nums),
      # tNum=factor(tNum, levels = target_nums),
      # Add integer index column using named vector
      # index=as.integer(target_idx[tName])
      )
    # arrange(qName)
  return(df)
}
# Pull out maximal matching contigs, where most of query contig maps onto target
maxMatches <- function(match_sums) {
  max_matches <- match_sums %>% group_by(query, target, qName) %>%
    filter(qPercent == max(qPercent)) %>% ungroup %>%
    # arrange(query, target, tNum)
    arrange(desc(tSize), desc(qSize))
  return(max_matches)
}
# Pivot PSL summary data frame for ggplot2 heatmap
pivotPsl <- function(df) {
  df_p <- df %>%
    pivot_wider(id_cols = qName,
                names_from = tName,
                values_from = qPercent,
                values_fill = 0) %>%
    pivot_longer(!qName,
                 names_to = "tName",
                 values_to = "qPercent") %>%
    rowwise %>%
    mutate(query = abbrevSpc(gsub(";.*", "", qName)),
           target = abbrevSpc(gsub(";.*", "", tName))) %>%
    ungroup %>%
    # Preserve factor order
    # mutate(tName = factor(tName, levels = levels(df$tName)))
    mutate(qName = factor(qName, levels = levels(df$qName)))
  return(df_p)
}
# Order query ID factors by max matches
plotOrder <- function(match_sums, max_matches) {
  h_order_y <- max_matches %>% arrange(desc(tSize), desc(qSize)) %>%
    # arrange(target, desc(tSize), desc(qSize)) %>%
    pull(qName) %>% as.character %>% unique
  h_order_x <- max_matches %>% arrange(desc(tSize), desc(qSize)) %>%
    pull(tName) %>% as.character %>% unique
  match_sums <- match_sums %>%
    filter(tName %in% max_matches$tName, qName %in% max_matches$qName) %>%
    mutate(qName = factor(x = qName, levels = h_order_y),
           tName = factor(x = tName, levels = h_order_x),
           qindex = as.integer(qName))
  return(match_sums)
}
old_heatPsl <- function(df) {
  # R base heatmap version
  fname <- paste0(match_name, "_heatmap.png")
  h <- heatmap(mat, scale = "row", keep.dendro = T,
               main = ttl, xlab = x_label, ylab = y_label, col = h_colors(25))
  hclust_mat <- as.hclust(h$Rowv)
  png(file = fname, units = "in", width = 5, height = 5, res = 300)
  heatmap(mat, main = ttl, xlab = x_label, ylab = y_label, col = h_colors(25),
          scale = "row")
  dev.off()
}
# Plot heatmap of given matrix
heatPsl <- function(match_sums_piv) {
  leg_name <- paste("Synteny", "coverage", sep = "\n")
  match_sums_piv <- match_sums_piv %>%
    # Filter out artificial chromosomes
    filter(fixChrom(tName) != "0", fixChrom(qName) != "0")
  # x_label <- abbrevSpc(spc_int)
  # y_label <- abbrevSpc("Saccharina_japonica")
  # ttl <- paste(x_label, "vs", y_label)
  # if (x_label == abbrevSpc(spc_int)) {
  #   p <- ggplot(data = match_sums_piv,
  #               mapping = aes(fill = qPercent, x = tName, y = qName)) +
  #     scale_y_continuous(expand = c(0, 0), breaks = breaks_width(100))
  # } else {
  #   p <- ggplot(data = match_sums_piv,
  #               mapping = aes(fill = qPercent, x = tName, y = qName))
  # }
  p <- ggplot(data = match_sums_piv,
              mapping = aes(fill = qPercent, x = tName, y = qName)) +
    geom_tile() +
    scale_fill_viridis_c(option = "turbo", labels = label_percent(),
                         name = leg_name) +
    # facet_grid(rows = vars(query), cols = vars(target),
    #            scales = "free", space = "free", switch = "both") +
    # labs(x = x_label, y = y_label) +
    theme(axis.title = element_text(face = "italic"),
          axis.text = element_blank(),
          legend.position = "left")

  return(p)
}
# Function to use ggsave on plots
plotSave <- function(plot_name, plot_type, plot_list, outdir = NULL,
                     width, height) {
  p <- plot_list[[plot_name]]
  fname <- paste0(plot_type, "_", plot_name, ".png")
  if (dir.exists(outdir)) fname <- paste0(outdir, fname)
  ggsave(filename = fname, plot = p, bg = "white",
         width = width, height = height, units = "in")
  print(fname)
  return(fname)
}

# Read in PSL summary files
match_sums <- read.table(match_sums_file, header = T, sep = "\t")
match_sums <- match_sums %>% rowwise %>%
  # Make unique column ID for each chr Num
  mutate(tID=paste(target, tNum, sep = "."),
         qID=paste(query, qNum, sep = "."),
         # Format species names
         query=abbrevSpc(query),
         target=abbrevSpc(target)) %>% ungroup
# Order species by decreasing relatedness to species of interest
spc_lvls <- unname(sapply(spc_order, grep,
                          unique(c(match_sums$target, match_sums$query)),
                          value = T))
match_sums <- match_sums %>%
  mutate(query=factor(query, levels = spc_lvls),
         target=factor(target, levels = spc_lvls))
# Keep only species of interest
match_sums <- match_sums %>% filter(query == abbrevSpc(spc_int))
# Get maximal match information per query scaffold
match_sums <- match_sums %>% group_by(target, qName) %>%
  mutate(max_match=case_when(qPercent == max(qPercent) ~ 1,
                             .default = 0)) %>% ungroup
max_matches <- match_sums %>%
  group_by(target, qName) %>% filter(qPercent == max(qPercent)) %>% ungroup

perc_heat <- max_matches %>% filter(tNum != 0) %>% 
  select(qNum, qPercent, target) %>%
  pivot_wider(id_cols = qNum,
              names_from = target,
              values_from = qPercent,
              values_fill = 0) %>%
  column_to_rownames(var = "qNum") %>%
  as.matrix
heatmap(perc_heat, scale = "none")

scaf_groups <- max_matches %>%
  filter(tNum != 0) %>%
  select(qNum, tNum, target) %>%
  pivot_wider(id_cols = qNum,
              names_from = target,
              values_from = tNum)


scaf_by_chrom <- match_sums %>%
  filter(tNum != 0) %>%
  select(qID, tID, max_match) %>%
  pivot_wider(id_cols = tID,
              names_from = qID,
              values_from = max_match,
              values_fill = 0) %>%
  column_to_rownames(var = "tID") %>%
  as.matrix
# hc <- heatmap(crossprod(scaf_by_chrom), symm = T, keep.dendro = T)
hc <- hclust(dist(crossprod(scaf_by_chrom)))
my_clust <- cutree(hc, k = 32)
match_sums <- match_sums %>% mutate(clust=my_clust[qID])

df <- match_sums %>% select(qNum, qSize, qPercent, clust) %>% unique %>%
  mutate(clust=as.character(clust))
match_sums <- match_sums %>% filter(tNum != 0) %>%
  mutate(clust=as.character(clust))
ggplot(data = df,
       mapping = aes(x = qSize, y = qPercent,
                     color = clust)) +
  scale_x_log10() +
  scale_y_log10() +
  geom_point()

# match_sums <- orderPsl(match_sums)
# max_matches <- maxMatches(match_sums)
# match_sums_pivs <- pivotPsl(match_sums)
# match_sums_pivs_ord <- plotOrder(match_sums_pivs, max_matches)
q_ord <- sizeSort(match_sums, "qName", "qSize")
t_ord <- sizeSort(match_sums, "tName", "tSize")
match_sums <- match_sums %>%
  mutate(qName=factor(qName, levels = q_ord, ordered = T),
         tName=factor(tName, levels = t_ord, ordered = T),
         qIndex=as.numeric(qName))
h_order <- max_matches %>%
  arrange(tName) %>% pull(qName) %>% as.character %>% unique

match_sums_pivs <- match_sums %>%
  pivot_wider(id_cols = qName,
              names_from = tName,
              values_from = qPercent,
              values_fill = 0) %>%
  pivot_longer(!c(qName, qIndex),
               names_to = "tName",
               values_to = "qPercent") %>%
  rowwise %>%
  mutate(query = abbrevSpc(gsub(";.*", "", qName)),
         target = abbrevSpc(gsub(";.*", "", tName))) %>%
  ungroup %>%
  # Preserve factor order
  mutate(tName = factor(tName, levels = levels(match_sums$tName)))
match_sums_pivs <- match_sums_pivs %>%
  filter(
    tName %in% max_matches$tName,
    qName %in% max_matches$qName
  ) %>%
  mutate(qName = factor(x = qName, levels = h_order, ordered = T))
heatPsl(match_sums_pivs)
 
# match_sums <- read.table(match_sums_file, header = T, sep = "\t")
# max_matches <- maxMatches(match_sums)
# # Convert scaffold names to size-ordered factors
# match_sums_ord <- orderPsl(match_sums)
# max_matches_ord <- orderPsl(max_matches)
# # Pivot PSL summary data frame for ggplot2 heatmap
# match_sums_piv <- pivotPsl(match_sums_ord)
# match_sums_piv <- plotOrder(match_sums_piv, max_matches_ord)

# # Plots
# (psl_heat <- heatPsl(match_sums_piv))
# # Heatmaps of genome vs. genome synteny
# print("Saving heatmaps...")
# h_fnames <- unlist(sapply(names(psl_heats), plotSave, "heatmap", psl_heats,
#                           outdir, 9, 9, simplify = F))

