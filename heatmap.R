## Initialization
rm(list = ls())
# Load required packages
library(scales, quietly = T, warn.conflicts = F)
library(cooccur)
library(visNetwork)
library(dendextend)
suppressPackageStartupMessages(library(ComplexHeatmap, quietly = T))
suppressPackageStartupMessages(library(tidyverse, quietly = T, warn.conflicts = F))
library(ggpubr)
library(paletteer)

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
# heatmap(perc_heat, scale = "none")

scaf_groups <- max_matches %>%
  filter(tNum != 0) %>%
  select(qNum, tNum, target) %>%
  pivot_wider(id_cols = qNum,
              names_from = target,
              values_from = tNum)

# Query scaffolds correlation
# scafs_perc <- match_sums %>%
#   # filter(tNum != 0) %>%
#   select(qID, tID, qPercent) %>%
#   pivot_wider(id_cols = qID,
#               names_from = tID,
#               values_from = qPercent,
#               values_fill = 0) %>%
#   column_to_rownames(var = "qID") %>%
#   as.matrix
scafs_max <- match_sums %>%
  filter(tNum != 0) %>%
  select(qID, tID, max_match) %>%
  pivot_wider(id_cols = qID,
              names_from = tID,
              values_from = max_match,
              values_fill = 0) %>%
  column_to_rownames(var = "qID") %>% as.matrix %>% tcrossprod
scafs_max_0 <- match_sums %>%
  # filter(tNum != 0) %>%
  select(qID, tID, max_match) %>%
  pivot_wider(id_cols = qID,
              names_from = tID,
              values_from = max_match,
              values_fill = 0) %>%
  column_to_rownames(var = "qID") %>% as.matrix %>% tcrossprod
# Heatmap with k splits showing clustering
heatClust <- function(mat, k) {
  h_labs <- sapply(rownames(mat), fixChrom)
  # hc <- hclust(dist(mat))
  # h <- heatmap(mat, scale = "none", keep.dendro = T,
  #              labRow = h_labs, labCol = h_labs, revC = F)
  # hclust_mat <- as.hclust(h$Rowv)
  # return(hclust_mat)
  pal <- circlize::colorRamp2(
    breaks = 0:4,
    reverse = T,
    # hcl_palette = "Grays",
    hcl_palette = "Inferno",
  )
  linecol <- "grey"
  row_dend <- hclust(dist(mat)) # row clustering
  col_dend <- hclust(dist(t(mat))) # column clustering
  h2 <- Heatmap(mat, name = "my_heat",
                show_row_names = F, show_column_names = F,
                border = linecol, col = pal,
                row_split = k, column_split = k,
                row_gap = unit(0, "mm"), column_gap = unit(0, "mm"),
                heatmap_legend_param = list(title = "Co-ocurrence"),
                row_title_gp = gpar(fontsize = 10),
                column_title_gp = gpar(fontsize = 10),
                # column_title_gp = gpar(fontsize = 10, col = rainbow(k)),
                # column_title_gp = gpar(fontsize = 10, col = rainbow(k)),
                # cluster_rows = color_branches(row_dend, k = k),
                # cluster_columns = color_branches(col_dend, k = k),
                )
  # for (i in 1:k) {
  #   decorate_heatmap_body("my_heat", row_slice = i, column_slice = i, {
  #     grid.rect(gp = gpar(col = rainbow(32)[i], fill = NA))
  #   })
  # }
  return(h2)
}
h <- heatClust(scafs_max, 32)
h_0 <- heatClust(scafs_max_0, 32)

getKH <- function(mat) {
  hc <- hclust(dist(mat))
  dend <- as.dendrogram(hc)
  kh <- heights_per_k.dendrogram(dend)
  kh <- tibble(k=names(kh), height=unname(kh)) %>%
    mutate(across(everything(), .fns = as.numeric)) %>%
    filter(height >= 0)
  return(kh)
}
kh_list <- list(no_chr0=getKH(scafs_max), all=getKH(scafs_max_0))
elbowK <- function(kh_list) {
  df <- bind_rows(kh_list, .id = "filter")
  p <- ggplot(data = df, mapping = aes(x = k, y = height, col = filter)) +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = scales::breaks_width(5), expand = c(0, 0),
                       limits = c(0, 200)) +
    # scale_x_log10() +
    # scale_y_log10() +
    geom_vline(xintercept = 39, col = "grey", lty = 2) +
    geom_vline(xintercept = 27, col = "grey", lty = 2) +
    theme_classic()
  return(p)
  # ggdendrogram(dend, leaf_labels = F, labels = F)
}
elbowK(kh_list)

hc1 <- hclust(dist(scafs_max))
clust <- cutree(hc1, k = 32)

# p1 <- ggdendro::ggdendrogram(dend, leaf_labels = F, labels = F)
# p2 <- ggdendro::ggdendrogram(dend_0, leaf_labels = F, labels = F)
# as.data.frame(heights_per_k.dendrogram(dend))
# plot(heights_per_k.dendrogram(dend_0))
# ggarrange(p1, p2)
# 
# plot(hclust_mat, labels = F)
# abline(h = 25)
clust <- cutree(
  hclust_mat,
  k = 32,
  # h = 25,
)

prepDf <- function(df) {
  df <- df %>% filter(!grepl("contig", tName), tNum <= 36)
  ord_tID <- df %>% 
    group_by(tID, clust) %>% summarize(n = n(), .groups = "drop") %>%
    pivot_wider(id_cols = tID,
                names_from = clust,
                values_from = n,
                values_fill = 0) %>%
    gather(key = "clust", value = "n", !tID) %>%
    group_by(tID) %>%
    mutate(likely_assignment=clust[which.max(n)],
           assignment_prob=max(n)) %>%
    arrange(as.integer(likely_assignment), desc(assignment_prob)) %>% 
    ungroup
  # ord_tID_1 <- ord_tID %>% rowwise %>% filter(fixChrom(tID) != 0) %>%
  #   pull(tID) %>% unique
  # ord_tID_0 <- ord_tID %>% rowwise %>% filter(fixChrom(tID) == 0) %>%
  #   pull(tID) %>% unique
  # ord_tID <- c(ord_tID_1, ord_tID_0)
  ord_tID <- ord_tID %>% pull(tID) %>% unique
  df <- df %>%
    mutate(tID=factor(tID, levels = ord_tID)) %>%
    arrange(clust, tID) %>%
    mutate(qID=factor(qID, levels = unique(qID)))
  return(df)
}
pivDf <- function(df) {
  df_p <- df %>%
    pivot_wider(id_cols = qID,
                names_from = tID,
                # values_from = clust) %>%
                values_from = qPercent,
                values_fill = 0) %>%
    pivot_longer(!qID,
                 names_to = "tID",
                 # values_to = "clust") %>%
                 values_to = "qPercent") %>%
    rowwise %>%
    mutate(query=abbrevSpc(gsub("\\..*", "", qID)),
           target=abbrevSpc(gsub("\\..*", "", tID)),
           target=factor(target, levels = spc_lvls)) %>%
    ungroup %>%
    # Preserve factor ordering
    mutate(tID=factor(tID, levels = levels(df$tID)),
           qID=factor(qID, levels = rev(levels(df$qID))))
  return(df_p)
}
barClust <- function(df, pal) {
  k <- length(unique(df$clust))
  leg_name <- paste("k = ", k)
  pal <- colorRampPalette(paletteer_d(pal))
  df <- df %>% summarize(n=n(), .by = c(tID, clust, query, target))
  # Convert cluster from integer to factor
  df <- df %>% mutate(clust=as.factor(clust))
  p <- ggplot(data = df, mapping = aes(x = tID, y = n, fill = clust)) +
    geom_col(width = 1, position = "fill") +
    scale_x_discrete(labels = as_labeller(fixChrom)) +
    scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
    scale_fill_manual(values = pal(k), name = leg_name) +
    facet_grid(rows = vars(query), cols = vars(target),
               scales = "free", switch = "both") +
    theme_minimal() +
    theme(strip.text = element_text(face = "italic"),
          strip.placement = "outside",
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 90))
  return(p)
}
heatClust <- function(df_p, pal) {
  leg_name <- paste("Synteny", "coverage", sep = "\n")
  # pal <- colorRampPalette(paletteer_d(pal))
  p <- ggplot(data = df_p,
              mapping = aes(x = tID, y = qID,
                            fill = qPercent)) +
                            # fill = clust)) +
    geom_tile() +
    scale_x_discrete(labels = as_labeller(fixChrom)) +
    scale_y_discrete(labels = as_labeller(fixChrom)) +
    # scale_fill_paletteer_c(name = leg_name, palette = pal) +
    # scale_fill_paletteer_d(name = leg_name, palette = pal) +
    # scale_fill_gradient(name = leg_name, low = "black", high = "white") +
    scale_fill_viridis_c(name = leg_name, option = "turbo",
                         labels = scales::label_percent()) +
    facet_grid(rows = vars(query), cols = vars(target),
               scales = "free", switch = "both") +
    theme_minimal() +
    theme(strip.text = element_text(face = "italic"),
          strip.placement = "outside",
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 90))
  return(p)
}
makePlots <- function(df, df_p, pal_b, pal_h) {
  p_bar <- barClust(df, pal_b)
  p_heat <- heatClust(df_p, pal_h)
  p <- ggarrange(p_bar, p_heat, nrow = 2, labels = "AUTO", align = "hv")
  return(p)
}

match_sums <- match_sums %>%
  mutate(clust=clust[qID]) %>%
  filter(!is.na(clust))
  # mutate(clust=as.factor(clust))
df <- match_sums
# df <- df %>% filter(tNum != 0)
# df <- df %>% filter(max_match == 1)
df <- prepDf(df)
df_p <- pivDf(df)
makePlots(df, df_p, "khroma::smoothrainbow", "grDevices::Oslo")
# pal <- "dichromat::BluetoGray_8"
# pal <- "colorBlindness::Blue2Orange8Steps"
# pal <- "grDevices::Grays"
# pal <- "ggprism::black_and_white"


# clust_cross <- match_sums %>%
#   filter(!is.na(clust),
#          max_match == 1) %>%
#   select(tNum, clust) %>%
#   pivot_wider(id_cols = tNum,
#               names_from = clust,
#               )
  # as.matrix %>% crossprod



# match_sums <- orderPsl(match_sums)
# max_matches <- maxMatches(match_sums)
# match_sums_pivs <- pivotPsl(match_sums)
# match_sums_pivs_ord <- plotOrder(match_sums_pivs, max_matches)
# q_ord <- sizeSort(match_sums, "qName", "qSize")
# t_ord <- sizeSort(match_sums, "tName", "tSize")
# match_sums <- match_sums %>%
#   mutate(qName=factor(qName, levels = q_ord, ordered = T),
#          tName=factor(tName, levels = t_ord, ordered = T),
#          qIndex=as.numeric(qName))
# h_order <- max_matches %>%
#   arrange(tName) %>% pull(qName) %>% as.character %>% unique
# 
# match_sums_pivs <- match_sums %>%
#   pivot_wider(id_cols = qName,
#               names_from = tName,
#               values_from = qPercent,
#               values_fill = 0) %>%
#   pivot_longer(!c(qName, qIndex),
#                names_to = "tName",
#                values_to = "qPercent") %>%
#   rowwise %>%
#   mutate(query = abbrevSpc(gsub(";.*", "", qName)),
#          target = abbrevSpc(gsub(";.*", "", tName))) %>%
#   ungroup %>%
#   # Preserve factor order
#   mutate(tName = factor(tName, levels = levels(match_sums$tName)))
# match_sums_pivs <- match_sums_pivs %>%
#   filter(
#     tName %in% max_matches$tName,
#     qName %in% max_matches$qName
#   ) %>%
#   mutate(qName = factor(x = qName, levels = h_order, ordered = T))
# heatPsl(match_sums_pivs)
 
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

