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
           qSize = qSize/1000000,
           qStart = qStart/1000000,
           qEnd = qEnd/1000000,
           tSize = tSize/1000000,
           tStart = tStart/1000000,
           tEnd = tEnd/1000000)
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
    geom_point(mapping = aes(x = qSize, y = tSize),
               alpha = 0) +
    geom_point(mapping = aes(x = Zeros, y = Zeros),
               alpha = 0) +
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
# lapply(names(psl_mats), heatPsl, psl_mats)
# Subset data for dotplots
psl_dot <- sapply(names(psl_list), dotFilt, psl_list, psl_match, simplify = F)

# Dotplots of genome vs. genome
p_list <- sapply(names(psl_dot), dotPlot, psl_dot, simplify = F)
seps <- sapply(names(psl_dot), sepDots, psl_dot, simplify = F)
sep_list <- sapply(names(seps), plotList, seps, simplify = F)
sapply(names(p_list), plotSave, p_list, 15, 10, simplify = F)
# sapply(names(sep_list), plotSave, sep_list, 7, 7, simplify = F)


# test <- psl_dot[[1]] %>%
#   filter(tName == "JABAKD010000001.1")
# ggplot(data = test) +
#   facet_grid(cols = vars(tName), rows = vars(qName)) +
#   geom_segment(aes(x = 0, y = 0, xend = tSize, yend = 0)) +
#   geom_segment(aes(x = 0, y = 1, xend = qSize, yend = 1)) +
#   geom_segment(aes(x = tStart, y = 0, xend = qStart, yend = 1, color = strand),
#                alpha = 0.2)

# showtext_opts(dpi = 300)
# ggsave(filename = fname, plot = p, width = 30, height = 40, units = "in")
# showtext_opts(dpi = 100)