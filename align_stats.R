## Initialization
rm(list = ls())
# Load required packages
# Required packages
# library(Hmisc, quietly = T, warn.conflicts = F)
library(scales, quietly = T, warn.conflicts = F)
suppressPackageStartupMessages(library(tidyverse, quietly = T,
                                       warn.conflicts = F))
library(VennDiagram, quietly = T)
library(ggVennDiagram, quietly = T)
library(ggpubr, quietly = T, warn.conflicts = F)
suppressPackageStartupMessages(library(ggpmisc, quietly = T,
                                       warn.conflicts = F))
library(RColorBrewer, quietly = T)
library(BiocManager, quietly = T)
if (require(showtext, quietly = T)) {
  showtext_auto()
  if (interactive()) showtext_opts(dpi = 100) else showtext_opts(dpi = 300)
}

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
lens_file <- "lens_by_chrom.tsv"
align_report_file <- "alignment_report.tsv"
# Output
venn_file <- "align_venn_diagram.png"
align_plot_file <- "align_length_and_n.png"
# Prepend output directory to file names (if it exists)
if (dir.exists(outdir)) {
  match_sums_file <- paste0(outdir, match_sums_file)
  max_match_file <- paste0(outdir, max_match_file)
  lens_file <- paste0(outdir, lens_file)
  align_report_file <- paste0(outdir, align_report_file)
  venn_file <- paste0(outdir, venn_file)
  align_plot_file <- paste0(outdir, align_plot_file)
}

# df <- read.table(max_match_file, sep = "\t", header = T)
# 
# uniq_n <- length(unique(df$qNum))
# uniq_len <- sum(unique(df[,c("qNum", "qSize")])[,"qSize"])*1e-6
# uniq_len_perc <- (uniq_len/615)*100
# uniq_match <- sum(unique(df[,c("qNum", "total_matches")])[,"total_matches"])*1e-6
# uniq_match_perc <- (uniq_match/615)*100
# 146.5/615

# Functions
# Total genome length of species of interest
spc_int_len <- 615545555
# Abbreviate species genus name
abbrevSpc <- function(spc) {
  spc <- unlist(strsplit(spc, "_| "))
  let1 <- substr(spc[1], 1, 1)
  spc[1] <- paste0(let1, ".")
  spc_a <- paste(spc, collapse = " ")
  return(spc_a)
}
# Clean up data frame for plotting
cleanDf <- function(df) {
  # Order "Species" column by decreasing relatedness
  spc_order <- c("japonica", "pyrifera", "pinnatifida", "siliculosus")
  lvls <- unname(sapply(spc_order, grep,
                        unique(df$Species), value = T))
  df_clean <- df %>% mutate(Species=factor(Species, levels = lvls),
                            # Convert bp to Mb
                            `Reference chromosome length (Mb)`=tSize*1e-6,
                            `Homologous scaffolds (Mb)`=sum_homolog*1e-6,
                            `Exact matches (Mb)`=sum_match*1e-6,
                            # Change Species variable name for legend
                            Reference=Species)
  # Catch spaces converted to periods in read.table
  colnames(df_clean) <- gsub("\\.", " ", colnames(df_clean))
  return(df_clean)
}
# Summarize alignment statistics between species of interest and all references
homoOverlap <- function(match_sums) {
  df_filt <- match_sums %>% filter(query == spc_int) %>%
    rowwise() %>% mutate(Species=abbrevSpc(target))
  df_venn <- df_filt %>% group_by(Species) %>%
    summarize(uniq_qName=list(unique(qName)), n=length(unique(qName)),
              .groups = "drop")
  un_n <- df_filt %>% select(qName, qSize) %>% unique %>% pull(qName) %>% length
  un_len <- df_filt %>% select(qName, qSize) %>% unique %>% pull(qSize) %>% sum
  un_len_perc <- un_len/spc_int_len*100
  ttl <- abbrevSpc(spc_int)
  subttl <- paste0(un_n, " unique homologous scaffolds (", round(un_len*1e-6, 2), " Mb)")
  v1 <- venn.diagram(df_venn$uniq_qName, category.names = df_venn$Species,
                     fill = rainbow(length(df_venn$Species)),
                     main = ttl, sub = subttl, print.mode = "raw",
                     main.cex = 2,
                     # Italicize species names
                     main.fontface = "italic", cat.fontface = "italic",
                     filename = NULL, disable.logging = T)
  png(filename = venn_file, res = showtext_opts()$dpi,
      width = 7, height = 5, units = "in")
  grid.newpage()
  grid.draw(v1)
  dev.off()
}
# Plot length of scaffold matches vs. chromosome length
plotLens <- function(df) {
  # Clean up data frame for plotting
  df <- cleanDf(df)
  # Filter out artificial chromosomes
  df <- df %>% filter(tNum != 0)
  # Function for individual plots
  lenPlot <- function(my_var, lab_pos) {
    y_label <- gsub("(\\(.+\\))",
                    paste0("in \\*", abbrevSpc(spc_int), "\\* \\1"),
                    my_var)
    p <- ggplot(data = df,
                mapping = aes(x = `Reference chromosome length (Mb)`,
                              y = .data[[my_var]], group = Reference,
                              col = Reference, fill = Reference)) +
      geom_point(alpha = 0.5) +
      stat_poly_line(formula = y~x+0, alpha = 0.1) +
      stat_poly_eq(formula = y~x+0, mapping = use_label(labels = c("R2", "eq")),
                   label.x = lab_pos[["x"]], label.y = lab_pos[["y"]]) +
      scale_x_continuous(breaks = pretty_breaks()) +
      scale_y_continuous(breaks = pretty_breaks()) +
      ylab(y_label) +
      theme_bw() +
      theme(legend.text = element_text(face = "italic"),
            axis.title.y = ggtext::element_markdown())
    return(p)
  }
  p1 <- lenPlot("Homologous scaffolds (Mb)", c(x = "right", y = "bottom"))
  p2 <- lenPlot("Exact matches (Mb)", c(x = "left", y = "top"))
  p <- ggarrange(p1, p2, nrow = 2, labels = "AUTO", common.legend = T,
                 legend = "left", align = "hv")
  return(p)
}
# Boxplots of n aligned per species
plotN <- function(df) {
  # Clean up data frame for plotting
  df <- cleanDf(df)
  y_label <- paste0("*", abbrevSpc(spc_int), "* scaffolds aligned per chromosome")
  p <- ggplot(data = df, mapping = aes(x = Reference, y = `n homologous scaffolds`,
                                       col = Reference, fill = Reference)) +
    geom_point(show.legend = F) +
    geom_boxplot(alpha = 0.5, width = 0.3,
                 outlier.alpha = 1, outlier.shape = 24, outlier.size = 2,
                 show.legend = F) +
    scale_y_continuous(breaks = breaks_width(20)) +
    ylab(y_label) +
    theme_bw() +
    theme(axis.text.x = element_text(face = "italic"),
          axis.title.y = ggtext::element_markdown())
  # If highest point is outlier, separate plot into high and low range
  rng <- diff(range(df$`n homologous scaffolds`))
  max_y <- sort(df$`n homologous scaffolds`, decreasing = T)[1]
  max_y_2 <- sort(df$`n homologous scaffolds`, decreasing = T)[2]
  if ((max_y - max_y_2)/rng > 0.5) {
    low_rng <- p + coord_cartesian(ylim = c(0, max_y_2))
    hi_rng <- p + coord_cartesian(ylim = c(max_y-10, max_y+10)) +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
    p <- ggarrange(hi_rng, low_rng, nrow = 2, heights = c(1, 5), align = "v")
  }
  # p <- annotate_figure(p, fig.lab = "C", fig.lab.face = "bold")
  return(p)
}

# Read in PSL summary files
match_sums <- read.table(match_sums_file, header = T, sep = "\t")
max_matches <- read.table(max_match_file, header = T, sep = "\t")
max_match_lens_sum <- read.table(lens_file, header = T, sep = "\t")

# Plots
# Venn diagram of species of interest homologous scaffolds versus references
homoOverlap(match_sums)
# Length of alignment by chromosome and species
lens_plot <- plotLens(max_match_lens_sum)
# n aligned by chromosome and species
n_box <- plotN(max_match_lens_sum)
n_box
# fig <- ggarrange(lens_plot, n_box, ncol = 2, labels = c("", "C"),
#                  widths = c(1, 0.8))
# # Save plots
# showtext_opts(dpi = 300)
# ggsave(filename = align_plot_file, plot = fig, bg = "white",
#        width = 10, height = 8, units = "in")
# showtext_opts(dpi = 100)
