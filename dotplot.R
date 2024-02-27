# Clear environment
rm(list = ls())
# Required packages
suppressPackageStartupMessages(library(tidyverse,
                                       quietly = T, warn.conflicts = F))
# library(caret)
library(ggpubr, quietly = T)
library(RColorBrewer, quietly = T)
suppressPackageStartupMessages(library(ComplexHeatmap, quietly = T))
# library(gridExtra, quietly = T)
if (require(showtext, quietly = T)) {
  showtext_auto()
  if (interactive()) showtext_opts(dpi = 100) else showtext_opts(dpi = 300)
}

# Variables and functions
# Set column names for PSL data frames
psl_col <- c("matches", "misMatches", "repMatches", "nCount", "qNumInsert",
             "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName",
             "qSize", "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd",
             "blockCount", "blockSizes", "qStarts", "tStarts")
# Gets genome name from "genome_vs_genome" character named: query_vs_target
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
  contigs <- gsub("\\..*", "", gsub("JABAKD01", "scaffold", contigs))
  return(contigs)
}
# Order PSL data frame by contig size by converting contig names to factors
orderPsl <- function(df_name, df_list) {
  df <- df_list[[df_name]] %>%
    mutate(qName=fixUndaria(qName),
           tName=fixUndaria(tName))
  query <- getGen(df_name, "query")
  target <- getGen(df_name, "target")
  # Convert scaffold names to factors ordered by size
  query_contigs <- df %>%
    select(qName, qSize) %>%
    arrange(desc(qSize)) %>%
    pull(qName) %>%
    unique()
  target_contigs <- df %>%
    select(tName, tSize) %>%
    arrange(desc(tSize)) %>%
    pull(tName) %>%
    unique()
  df <- df %>%
    mutate(qName=factor(qName, levels = query_contigs),
           tName=factor(tName, levels = target_contigs)) %>%
    # Filter out extremely small contigs
    filter(qSize > 1e6 & tSize > 1e6) %>%
    # Filter out extremely large contigs
    filter(qSize < 30e6 & tSize < 30e6) %>%
    arrange(qName)
  return(df)
}
# Calculate metrics from PSL table
metricsPsl <- function(df_name, df_list) {
  df <- df_list[[df_name]]
  query <- getGen(df_name, "query")
  target <- getGen(df_name, "target")
  df <- df %>%
    # Calculate synteny within each alignment range
    mutate(synt=matches/abs(qStart-qEnd)) %>%
    # Keep only highest synteny for each scaffold-scaffold pair
    group_by(qName, tName) %>%
    filter(synt == max(synt)) %>%
    ungroup()
  # %>%
  return(df)
    # # Calculate fragmentedness of alignment
    # mutate(p_blockCount=blockCount/abs(qStart-qEnd)) %>%
    # # Calculate maximum block size (bp) of each alignment
    # mutate(max_block=max(as.numeric(unlist(strsplit(blockSizes, ","))))) %>%
    # # Calculate each alignment size
    # mutate(aln_size=abs(qStart-qEnd))
  df_cov <- df %>%
    group_by(qName, tName) %>%
    summarize(total_aln = sum(aln_size),
              total_p_blockCount = sum(p_blockCount)) %>%
    ungroup()
  df <- merge(df, df_cov)
  df_scale <- df %>%
    # Calculate coverage of each alignment
    mutate(cov = total_aln/qSize) %>%
    # Scale metrics [0,1] for each query scaffold
    group_by(qName) %>%
    mutate(
      synt_scaled=
        (synt-min(synt))/(max(synt)-min(synt)),
      cov_scaled=
        (cov-min(cov))/(max(cov)-min(cov)),
      max_block_scaled=
        (max_block-min(max_block))/(max(max_block)-min(max_block)),
      blockCount_scaled=
        (blockCount-min(blockCount))/(max(blockCount)-min(blockCount)),
      p_blockCount_scaled=
             (p_blockCount-min(p_blockCount))/(max(p_blockCount)-min(p_blockCount)),
      total_p_blockCount_scaled=
        (total_p_blockCount-min(total_p_blockCount))/(max(total_p_blockCount)-min(total_p_blockCount))
      )
  # For query scaffolds with 1 match, scale to 1
  df_scale[is.na.data.frame(df_scale)] <- 1
  return(df_scale)
}
# Plot 
allSyntPlot <- function(df_name, df_list) {
  df <- df_list[[df_name]]
  query <- getGen(df_name, "query")
  target <- getGen(df_name, "target")
  df <- df %>%
    filter(qName %in% unique(df$qName)[1:50])
  p <- ggplot(data = df,
              # Plot alignment blocks along query scaffold
              mapping = aes(x = qStart, xend = qEnd,
                            # Label scaffolds on y-axis
                            y = tStart, yend = tEnd,
                            # # Shade line segments by degree of synteny
                            # alpha = synt,
                            # Shade line segments by relative coverage
                            alpha = cov_scaled,
                            # Color line segments by degree of fragmentedness
                            col = total_p_blockCount_scaled)) +
    # facet_grid(tName~qName) +
    # scale_alpha_continuous(trans = "reverse") +
    facet_grid(cols = vars(qName), rows = vars(tName),
               switch = "both",
               space = "free",
               scale = "free",
               as.table = F) +
    coord_cartesian(clip="off") +
    scale_color_continuous(type = "viridis") +
    geom_point(mapping = aes(x = qSize, y = tName),
               alpha = 0) +
    geom_segment(linewidth = 3, lineend = "round") +
    labs(x = query, y = target, title = "Genome synteny by scaffold") +
    theme_classic() +
    theme(axis.title = element_text(face = "italic"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.1),
          panel.spacing = unit(0, "lines"),
          # strip.clip = "off",
          strip.text.x.bottom = element_text(size = rel(0.75), angle = 90),
          strip.text.y.left = element_text(size = rel(0.75), angle = 0),
          strip.placement = "outside",
          # strip.text.y.left = element_text(angle = 0),
          strip.background.x = element_rect(color = NA,  fill=NA),
          strip.background.y = element_rect(color = NA,  fill=NA))
  return(p)
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
    mutate(qPercent=total_matches/qSize) %>%
    arrange(qName, desc(qPercent))
  return(df2)
}
# Pull out maximal matching contigs,
# where most of query contig maps onto target
maxMatches <- function(df_name, df_list) {
  df <- df_list[[df_name]]
  query <- getGen(df_name, "query")
  target <- getGen(df_name, "target")
  df_match <- df %>% 
    group_by(qName) %>%
    filter(qPercent == max(qPercent)) %>%
    ungroup()
  return(df_match)
}


# Convert PSL summary into matrix for heatmap
matPsl <- function(df_name, df_list) {
  df <- df_list[[df_name]]
  query <- getGen(df_name, "query")
  target <- getGen(df_name, "target")
  mat <- as.matrix(df %>%
                     pivot_wider(id_cols = qName,
                                 names_from = tName,
                                 values_from = qPercent,
                                 # values_from = total_matches,
                                 values_fill = 0) %>%
                     column_to_rownames(var = "qName"))
  return(mat)
}
# Plot heatmap of given matrix
heatPsl <- function(mat_name, mat_list, df_list) {
  mat <- mat_list[[mat_name]]
  df <- df_list[[mat_name]]
  query <- getGen(mat_name, "query")
  target <- getGen(mat_name, "target")
  x_label <- paste("target:", gsub("_", " ", target))
  y_label <- paste("query:", gsub("_", " ", query))
  ttl <- gsub("_", " ", mat_name)
  fname <- paste0(mat_name, "_heatmap.png")
  h_colors <- colorRampPalette(brewer.pal(8, "YlOrRd"))
  # h <- heatmap(mat, scale = "none", keep.dendro = T)
  h <- heatmap(mat, scale = "row", keep.dendro = T,
               main = ttl, xlab = x_label, ylab = y_label, col = h_colors(25))
  hclust_mat <- as.hclust(h$Rowv)
  # df$BUSCO <- factor(x = df$BUSCO,
  #                    levels = c(rownames(mat)[hclust_order], "Genes"),
  #                    ordered = T)
  df <- df %>%
    mutate(qName = factor(x = qName,
                      # levels = c(hclust_mat$labels[hclust_mat$order], "Genes"),
                      levels = hclust_mat$labels[hclust_mat$order],
                      ordered = T))
  p <- ggplot(mapping = aes(x = tName, y = qName)) +
    geom_tile(data = df,
              mapping = aes(fill = qPercent)) +
    # scale_fill_manual(values = h_colors) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45),
          axis.text.y = element_text(angle = 45),
          axis.line = element_line(color = "black"),
          legend.position = "left") +
    ggtitle(ttl)
  return(p)
  # png(file = fname, units = "in", width = 5, height = 5, res = 300)
  # heatmap(mat, main = ttl, xlab = x_label, ylab = y_label, col = h_colors(25),
  #         scale = "row")
  # dev.off()
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
  df_sub <- df_sub %>%
    group_by(tName) %>%
    arrange(tStart, qStart) %>%
    mutate(qName = factor(qName, levels = unique(qName))) %>%
    ungroup()
  return(df_sub)
}
# Dotplots
dotPlot <- function(df_name, df_list, filter_ids = NULL) {
  df <- df_list[[df_name]]
  if (missing(filter_ids)) {
  } else {
    if (all(filter_ids %in% df$tName)) {
      df <- df %>%
        filter(tName %in% filter_ids)
    } else {
      print(filter_ids %in% df$tName)
      print(filter_ids)
      print(class(filter_ids))
      print(unique(df$tName))
      print(class(df$tName))
      stop("Error: ensure all 'filter_ids' are present in 'tName' column.")
    }
  } 
  df <- df %>%
    mutate(Zeros = 0,
           qSize = qSize/1e6,
           qStart = qStart/1e6,
           qEnd = qEnd/1e6,
           tSize = tSize/1e6,
           tStart = tStart/1e6,
           tEnd = tEnd/1e6)
  target <- getGen(df_name, "target")
  target_nice <- gsub("_", " ", target)
  query <- getGen(df_name, "query")
  query_nice <- gsub("_", " ", query)
  fname <- paste0(df_name, "_dotplot.png")
  p <- ggplot(data = df,
              mapping = aes(x = qStart, y = tStart, color = strand)) +
    geom_segment(mapping = aes(xend = qEnd, yend = tEnd),
                 # linewidth = 2,
                 lineend = "round") +
    # geom_point(mapping = aes(x = qSize, y = tSize),
    #            alpha = 0) +
    # geom_point(mapping = aes(x = Zeros, y = Zeros),
    #            alpha = 0) +
    facet_grid(cols = vars(qName), rows = vars(tName),
               switch = "both",
               space = "free",
               scale = "free",
               as.table = F) +
    coord_cartesian(clip="off") +
    theme_classic() +
    theme(axis.title = element_text(face = "italic"),
          # strip.clip = "off",
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          strip.background.x = element_rect(color = NA,  fill=NA),
          strip.background.y = element_rect(color = NA,  fill=NA)) +
    xlab("") +
    ylab("")
  if (missing(filter_ids)) {
    p <- p +
      xlab(query_nice) +
      ylab(target_nice) +
      labs(title = "Genome synteny by scaffold") +
      theme(axis.ticks = element_blank(),
            axis.text = element_blank(),
            # axis.title = element_text(size = rel(3)),
            # legend.title = element_text(size = rel(3)),
            # legend.text = element_text(size = rel(3)),
            # plot.title = element_text(size = rel(4)),
            panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.1),
            panel.spacing = unit(0, "lines"),
            strip.text.x.bottom = element_text(size = rel(0.75), angle = 90),
            strip.text.y.left = element_text(size = rel(0.75), angle = 0))
            # strip.text.x.bottom = element_text(size = rel(2), angle = 45),
            # strip.text.y.left = element_text(size = rel(2), angle = 0))
  }
  # showtext_opts(dpi = 300)
  # ggsave(filename = fname, plot = p, width = 30, height = 40, units = "in")
  # showtext_opts(dpi = 100)
  return(p)
}
# Segmented synteny plots
sepDots <-  function(df_name, df_list) {
  df <- df_list[[df_name]]
  tIds <- levels(df$tName)[levels(df$tName) %in% unique(df$tName)]
  sep_list <- list()
  for (i in 1:length(tIds)) {
    p <- dotPlot(df_name, df_list, tIds[i])
    sep_list[[i]] <- p
  }
  return(sep_list)
}
# Plot list of ggplots in one figure
plotList <- function(plot_name, plot_list) {
  plots <- plot_list[[plot_name]]
  ttl <- gsub("_", " ", plot_name)
  p <- ggarrange(plotlist = plots, common.legend = T)
  p <- annotate_figure(p, fig.lab = ttl)
  return(p)
}
# Function to use ggsave on plots
plotSave <- function(plot_name, plot_list, width, height) {
  p <- plot_list[[plot_name]]
  fname <- paste0(plot_name, "_dotplot.png")
  showtext_opts(dpi = 300)
  ggsave(filename = fname, plot = p, 
         width = width, height = height, units = "in")
  showtext_opts(dpi = 100)
  return(fname)
}


# Input
# Only take command line input if not running interactively
if (interactive()) {
  wd <- "/scratch1/kdeweese/latissima/genome_stats"
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
# Import dataframe of species and associated file names
species_tab <- read.table(seq_file, sep = "\t", skip = 1)
species <- species_tab$V1
# Keep only unique species combinations
species_versus <- lapply(combn(rev(species), 2, simplify = F),
                         paste, collapse = "_vs_")
# # Designate species of interest as query genome in all combinations
# species_versus <- paste(spc_int,
#                         grep(spc_int, species, invert = T, value = T),
#                         sep = "_vs_")
# Create list of PSL filenames
psl_files <- grep(paste(species_versus, collapse = "|"),
                  list.files(full.names = T, pattern = "\\.psl"), value = T)
# Read PSL files into list of dataframes
psl_list <- sapply(psl_files, read.table, col.names = psl_col,
                   simplify = F, USE.NAMES = T)
# Remove file extensions from list names
names(psl_list) <- gsub(".*ment_|.psl", "", names(psl_list))
# Convert scaffold names to size-ordered factors
psl_list <- sapply(names(psl_list), orderPsl, psl_list, simplify = F)

# Analysis
# Change to output directory for plots (if it exists)
if (dir.exists(outdir)) setwd(outdir)
# Summarize by contig vs. contig of each syntenic comparison
psl_sums <- sapply(names(psl_list), sumPsl, psl_list, simplify = F)
# Select maximal matching contigs
psl_match <- sapply(names(psl_sums), maxMatches, psl_sums, simplify = F)
# Subset data for dotplots
psl_dot <- sapply(names(psl_list), dotFilt, psl_list, psl_match, simplify = F)
# Summarize subsetted data
psl_dot_sums <- sapply(names(psl_dot), sumPsl, psl_dot, simplify = F)
# Matrices of syntenic blocks by contig for heatmaps
psl_mats1 <- sapply(names(psl_sums), matPsl, psl_sums, simplify = F)
psl_mats2 <- sapply(names(psl_dot_sums), matPsl, psl_dot_sums, simplify = F)
# Heatmaps of genome vs. genome
psl_heats1 <- sapply(names(psl_mats1), heatPsl, psl_mats2, psl_sums, simplify = F)
psl_heats2 <- sapply(names(psl_mats2), heatPsl, psl_mats1, psl_dot_sums, simplify = F)
psl_heats[["Saccharina_latissima_vs_Macrocystis_pyrifera"]]
# gplots::heatmap.2(psl_mats[["Undaria_pinnatifida_vs_Macrocystis_pyrifera"]])


# # Calculate scaffold-scaffold alignment metrics in each PSL dataframe
# psl_metrics <- sapply(names(psl_list), metricsPsl, psl_list, simplify = F)
# test <- psl_metrics[["Saccharina_latissima_vs_Macrocystis_pyrifera"]]
# ggplot(data = test,
#        mapping = aes(x = total_p_blockCount_scaled,
#                      y = cov_scaled,
#                      col = qName)) +
#   geom_point() +
#   theme(legend.position = "none")
# test$euc_dist <- (1-test$cov_scaled)^2+(0+test$total_p_blockCount_scaled)^2
# test2 <- test %>%
#   group_by(qName) %>%
#   filter(euc_dist == min(euc_dist)) %>%
#   filter(qName == "scaffold_1")
# ggplot(data = test,
#        mapping = aes(x = total_p_blockCount_scaled,
#                      y = cov_scaled,
#                      col = qName)) +
#   geom_point() +
#   facet_wrap(~qName) +
#   theme_classic() +
#   theme(legend.position = "none")
# psl_plots <- sapply(names(psl_metrics), allSyntPlot, psl_metrics, simplify = F)
# allSyntPlot("Saccharina_latissima_vs_Macrocystis_pyrifera", psl_metrics)




# # Summarize by contig vs. contig of each syntenic comparison
# psl_sums <- sapply(names(psl_list), sumPsl, psl_list, simplify = F)
# # Select maximal matching contigs
# psl_match <- sapply(names(psl_sums), maxMatches, psl_sums, simplify = F)
# # Matrices of syntenic blocks by contig for heatmaps
# psl_mats <- sapply(names(psl_sums), matPsl, psl_sums, simplify = F)
# # Subset data for dotplots
# psl_dot <- sapply(names(psl_list), dotFilt, psl_list, psl_match, simplify = F)
# 
# # Plots
# # Heatmaps of genome vs. genome
# lapply(names(psl_mats), heatPsl, psl_mats)
# # Dotplots of genome vs. genome
# p_list <- sapply(names(psl_dot), dotPlot, psl_dot, simplify = F)
# fnames <- sapply(names(p_list), plotSave, p_list, 15, 10, simplify = F)
# print(fnames)
# # Separate out certain syntenic regions
# # seps <- sapply(names(psl_dot), sepDots, psl_dot, simplify = F)
# # sep_list <- sapply(names(seps), plotList, seps, simplify = F)
# # sapply(names(sep_list), plotSave, sep_list, 7, 7, simplify = F)
