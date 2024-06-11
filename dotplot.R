# Clear environment
rm(list = ls())
# Required packages
library(Hmisc)
library(scales)
suppressPackageStartupMessages(library(tidyverse, quietly = T, warn.conflicts = F))
library(ggpubr, quietly = T)
library(ggpmisc)
library(RColorBrewer, quietly = T)
library(BiocManager)
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
# Order PSL data frame by contig size by converting contig names to factors
orderPsl <- function(df_name, df_list) {
  df <- df_list[[df_name]] %>%
    mutate(qNum=fixChrom(qName),
           tNum=fixChrom(tName))
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
  query_nums <- df %>%
    select(qNum, qSize) %>%
    arrange(desc(qSize)) %>%
    pull(qNum) %>%
    unique()
  target_nums <- df %>%
    select(tNum, tSize) %>%
    arrange(desc(tSize)) %>%
    pull(tNum) %>%
    unique()
  target_index_df <- df %>%
    select(tNum, tSize) %>%
    arrange(desc(tSize)) %>%
    unique() %>%
    filter(tNum != "0") %>%
    rowid_to_column(var="index") %>%
    select(index, tNum)
  target_index <- c(target_index_df$index, 0)
  names(target_index) <- c(target_index_df$tNum, "0")
  df <- df %>%
    mutate(index=target_index[tNum],
           qName=factor(qName, levels = query_contigs),
           tName=factor(tName, levels = target_contigs),
           qNum=factor(qNum, levels = query_nums),
           tNum=factor(tNum, levels = target_nums)) %>%
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
# Summarize PSL table by summing matches for each contig pair
sumPsl <- function(df_name, df_list) {
  df <- df_list[[df_name]]
  query <- getGen(df_name, "query")
  target <- getGen(df_name, "target")
  df_sum <- df %>%
    group_by(index, qNum, tNum, qSize, tSize) %>%
    summarize(total_matches=as.numeric(sum(matches)), .groups = "keep")
  # df2 <- merge(df_sum, unique(df[,c("qNum", "qSize")]),
  #              by = "qNum", sort = F)
  # df2 <- df2 %>%
  df_sum <- df_sum %>% mutate(qPercent=total_matches/qSize) %>%
    arrange(qNum, desc(qPercent))
  return(df_sum)
}
# Pull out maximal matching contigs, where most of query contig maps onto target
maxMatches <- function(df_name, df_list) {
  df <- df_list[[df_name]]
  query <- getGen(df_name, "query")
  target <- getGen(df_name, "target")
  df_match <- df %>% 
    group_by(qNum) %>%
    filter(qPercent == max(qPercent)) %>%
    ungroup() %>%
    arrange(tNum)
  return(df_match)
}
# Order query ID factors by max matches
plotOrder <- function(match_name, match_list, df_list) {
  match <- match_list[[match_name]]
  df <- df_list[[match_name]]
  query <- getGen(match_name, "query")
  target <- getGen(match_name, "target")
  h_order <- match %>%
    arrange(tNum) %>%
    pull(qNum) %>%
    as.character()
  df <- df %>%
    mutate(qNum = factor(x = qNum, levels = h_order, ordered = T))
  return(df)
}
# Correlate alignments between all species
groupScafs <- function(match_list, spc_int) {
  m_names <- grep(spc_int, names(match_list), value = T)
  grp_list <- list()
  for (i in 1:length(m_names)) {
    m_name <- m_names[[i]]
    match <- match_list[[m_name]]
    query <- getGen(m_name, "query")
    target <- getGen(m_name, "target")
    grp <- match %>%
      select(qNum, tNum) %>%
      pivot_wider(names_from = tNum,
                  values_from = qNum)
    #   group_by(tNum) %>%
    #   group_rows()
    grp_list[[m_name]] <- grp
  }
  return(grp_list)
}
# Pivot PSL summary dataframe for ggplot2 heatmap
pivotPsl <- function(df_name, df_list) {
  df <- df_list[[df_name]]
  query <- getGen(df_name, "query")
  target <- getGen(df_name, "target")
  df_p <- df %>%
    pivot_wider(id_cols = qNum,
                names_from = tNum,
                values_from = qPercent,
                values_fill = 0) %>%
    pivot_longer(!qNum,
                 names_to = "tNum",
                 values_to = "qPercent") %>%
    mutate(tNum = factor(tNum, levels = levels(df$tNum)))
  return(df_p)
}
# # Plot number & length of scaffold matches vs. chromosome length (qcount)
# plotLens <- function(qcount_name, qcount_list) {
#   qcount <- qcount_list[[qcount_name]]
#   query <- getGen(qcount_name, "query")
#   target <- getGen(qcount_name, "target")
#   x_label <- "Reference chromosome length (Mb)"
#   y_label1 <- "Length of homologs (Mb)"
#   y_label2 <- "n homologs"
#   ttl <- paste(abbrevSpc(query), "alignment to", abbrevSpc(target))
#   # Convert bp to Mb
#   qcount <- qcount %>% mutate(tSize=tSize*1e-6) %>% mutate(n_sum=n_sum*1e-6)
#   # Limits for all scatter plots, scaled to Mb
#   x_max <- max(unlist(lapply(qcount_list, function(x) pull(x, tSize))))*1e-6
#   y_max <- max(unlist(lapply(qcount_list, function(x) pull(x, n_sum))))*1e-6
#   common_lims <- coord_cartesian(xlim = c(0, x_max + 1), ylim = c(0, y_max + 1))
#   # # Add zero point to force fit to x=0 y=0 intercept
#   # levels(qcount$tNum) <- c("0", levels(qcount$tNum))
#   # fakezero <- c(tNum=factor("0", levels = levels(qcount$tNum)),
#   #               tSize=as.integer(0), n=NA, n_sum=as.integer(0))
#   # qcount <- rbind(fakezero, qcount)
#   barwidth <- (resolution(qcount$tSize)*1000)/max(qcount$tSize)
#   p1 <- ggplot(data = qcount,
#               mapping = aes(x = tSize, y = n_sum)) +
#     geom_point() +
#     stat_poly_line(formula = y~x+0, se = F, fullrange = T) +
#     stat_poly_eq(mapping = use_label(c("R2", "eq")), formula = y~x+0) +
#     # common_lims +
#     labs(y = y_label1) +
#     theme_classic() +
#     theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
#           axis.line.x = element_blank(), axis.ticks.x = element_blank())
#   p2 <- ggplot(data = qcount,
#                mapping = aes(x = tSize, y = n)) +
#     geom_col(alpha = 0.2, position = "identity", width = barwidth) +
#     # coord_cartesian(xlim = c(0, x_max + 1)) +
#     labs(x = x_label, y = y_label2) +
#     theme_classic()
#   p <- ggarrange(p1, p2, nrow = 2, align = "v", heights = c(2, 1))
#   p <- annotate_figure(
#     p,
#     bottom = text_grob(abbrevSpc(target), face = "italic"),
#     left = text_grob(abbrevSpc(query), face = "italic", rot = 90)
#   )
#   return(p)
# }
# Collapse list of maxMatch data frames into one data frame labeling target species
collapseMatch <- function(match_list) {
  # Subset species of interest
  match_list <- match_list[grep(spc_int, names(match_list), value = T)]
  # Rename data frames by target species
  names(match_list) <- sapply(names(match_list), getGen, "target", USE.NAMES = F)
  # Format species names
  names(match_list) <- sapply(names(match_list), abbrevSpc, USE.NAMES = F)
  # Collapse list of data frames
  match_lens <- bind_rows(match_list, .id = "Species")
  # Convert species column to ordered factor sorted by species relatedness
  spc_order <- c("japonica", "pyrifera", "pinnatifida", "siliculosus")
  lvls <- unname(sapply(spc_order, grep, unique(match_lens$Species), value = T))
  match_lens <- match_lens %>%
    mutate(qNum=as.character(qNum),
           tNum=as.character(tNum),
           Species=factor(Species, levels = lvls)) %>%
    # Filter out artificial chromosomes
    filter(tNum != "0")
    # # Pivot dataframe longer to plot metrics on same graph
    # pivot_longer(!c(tNum, tSize, Species),
    #              names_to = "metric", values_to = "value")
  return(match_lens)
}
# Plot number and length of scaffold matches vs. chromosome length (match_lens)
plotLens <- function(match_lens) {
  match_lens_sum <- match_lens %>% group_by(tNum, tSize, Species) %>%
    summarize(sum_homolog=sum(qSize), sum_match=sum(total_matches),
              `n homologs`=n(), .groups = "keep") %>%
    # Convert bp to Mb
    mutate(`Reference chromosome length (Mb)`=tSize*1e-6,
           `summed homolog lengths (Mb)`=sum_homolog*1e-6,
           `summed match lengths (Mb)`=sum_match*1e-6,
           Reference=Species) %>%
    arrange(Species, desc(tSize)) %>% ungroup()
  lenPlot <- function(my_var) {
    p <- ggplot(data = match_lens_sum,
                mapping = aes(x = `Reference chromosome length (Mb)`,
                              group = Reference,
                              col = Reference, fill = Reference)) +
      geom_point(mapping = aes(y = .data[[my_var]]), alpha = 0.5) +
      stat_poly_line(mapping = aes(y = .data[[my_var]]), formula = y~x+0) +
      stat_poly_eq(formula = y~x+0, 
                   mapping = aes(y = .data[[my_var]],
                                 label = paste(after_stat(rr.label),
                                               after_stat(eq.label),
                                               sep = "*\", \"*"))) +
                   # geom = "label_npc", col = "black", alpha = 0.5) +
      scale_x_continuous(breaks = pretty_breaks()) +
      scale_y_continuous(breaks = pretty_breaks()) +
      theme_bw() +
      theme(legend.text = element_text(face = "italic"))
    return(p)
  }
  addSpcInt <- function(p) {
    p <- annotate_figure(
      p,
      left = text_grob(abbrevSpc(spc_int), face = "italic", rot = 90)
    )
    return(p)
  }
  p1 <- lenPlot("summed homolog lengths (Mb)")
  leg1 <- get_legend(p1)
  p1 <- addSpcInt(p1 + theme(legend.position = "none"))
  p2 <- lenPlot("summed match lengths (Mb)")
  p2 <- addSpcInt(p2 + theme(legend.position = "none"))
  p3 <- ggplot(data = match_lens %>% filter(index < 30),
               mapping = aes(x = index, y = total_matches,
                             fill = qNum, group = Species)) +
    # geom_col(col = "white") +
    geom_col(mapping = aes(col = qNum)) +
    stat_summary(mapping = aes(y = tSize), geom = "point",
                 fun = mean, col = "black") +
    # scale_x_continuous(breaks = pretty_breaks(),
    #                    labels = label_number(scale = 1e-6)) +
    scale_y_continuous(breaks = pretty_breaks(),
                       labels = label_number(scale = 1e-6)) +
    facet_grid(cols = vars(Species), switch = "both") +
    theme_bw() +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_text(face = "italic"),
          strip.placement = "outside")
  return(p3)
  p <- ggarrange(p1, p2, align = "hv", legend.grob = leg1, legend = "right")
  return(p)
}
# Plot heatmap of given matrix
heatPsl <- function(match_name, match_list, df_list) {
  match <- match_list[[match_name]]
  df <- df_list[[match_name]]
  query <- getGen(match_name, "query")
  target <- getGen(match_name, "target")
  x_label <- abbrevSpc(target)
  y_label <- abbrevSpc(query)
  ttl <- gsub("_", " ", match_name)
  # Filter out artificial chromosomes
  match <- match %>% filter(tNum != "0", qNum != "0")
    # # Filter out small scaffolds/contigs from query (remove <1Mb)
    # filter(qSize > 1e6)
  # Order of query contigs while their respective homologs (tNum) are size-ordered
  h_order <- match %>% arrange(tNum) %>% pull(qNum) %>% as.character()
  # Apply filtering from match data frame
  df <- df %>% filter(tNum %in% match$tNum, qNum %in% h_order) %>%
    # Use "match" data frame to order query contig factors in heatmap
    mutate(qNum = factor(x = qNum, levels = h_order, ordered = T))
  p <- ggplot(mapping = aes(x = tNum, y = qNum)) +
    geom_tile(data = df,
              mapping = aes(fill = qPercent)) +
    scale_fill_viridis_c(option = "turbo") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black"),
          legend.position = "left") +
    labs(x = x_label, y = y_label) +
    theme(axis.title = element_text(face="italic"))
  # # R base heatmap version
  # fname <- paste0(match_name, "_heatmap.png")
  # h <- heatmap(mat, scale = "row", keep.dendro = T,
  #              main = ttl, xlab = x_label, ylab = y_label, col = h_colors(25))
  # hclust_mat <- as.hclust(h$Rowv)
  # png(file = fname, units = "in", width = 5, height = 5, res = 300)
  # heatmap(mat, main = ttl, xlab = x_label, ylab = y_label, col = h_colors(25),
  #         scale = "row")
  # dev.off()
  return(p)
}
# Dotplots
dotPlot <- function(df_name, df_list, filter_ids = NULL) {
  df <- df_list[[df_name]]
  if (missing(filter_ids)) {
  } else {
    if (all(filter_ids %in% df$tNum)) {
      df <- df %>%
        filter(tNum %in% filter_ids)
    } else {
      print(filter_ids %in% df$tNum)
      print(filter_ids)
      print(class(filter_ids))
      print(unique(df$tNum))
      print(class(df$tNum))
      stop("Error: ensure all 'filter_ids' are present in 'tNum' column.")
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
              mapping = aes(x = tStart, y = qStart, color = strand)) +
    geom_segment(mapping = aes(xend = tEnd, yend = qEnd),
                 lineend = "round") +
    geom_point(mapping = aes(x = tSize, y = qSize),
               alpha = 0) +
    geom_point(mapping = aes(x = Zeros, y = Zeros),
               alpha = 0) +
    facet_grid(cols = vars(tNum), rows = vars(qNum),
               switch = "both",
               space = "free",
               scale = "free",
               as.table = F) +
    coord_cartesian(clip="off") +
    theme_classic() +
    theme(axis.title = element_text(face = "italic"),
          # strip.clip = "off",
          strip.placement = "outside",
          strip.background.x = element_rect(color = NA,  fill=NA),
          strip.background.y = element_rect(color = NA,  fill=NA)) +
    xlab("") +
    ylab("")
  if (missing(filter_ids)) {
    p <- p +
      labs(title = "Genome synteny by scaffold",
           x = target_nice, y = query_nice) +
      theme(axis.ticks = element_blank(),
            axis.text = element_blank(),
            panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.1),
            panel.spacing = unit(0, "lines"),
            strip.text.x.bottom = element_text(size = rel(0.75), angle = 0),
            strip.text.y.left = element_text(size = rel(0.75), angle = 0))
  }
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
plotSave <- function(plot_name, plot_type, plot_list, outdir = NULL,
                     width, height) {
  p <- plot_list[[plot_name]]
  fname <- paste0(plot_type, "_", plot_name, ".png")
  if (dir.exists(outdir)) {
      fname <- paste0(outdir, fname)
  }
  showtext_opts(dpi = 300)
  ggsave(filename = fname, plot = p,
         width = width, height = height, units = "in")
  showtext_opts(dpi = 100)
  return(fname)
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
# Import dataframe of species and associated file names
species_tab <- read.table(seq_file, sep = "\t", skip = 1)
species <- species_tab$V1
# Keep all combinations of remaining species
other_spc <- gtools::permutations(v = species[!species %in% spc_int],
                                  n = length(species[!species %in% spc_int]),
                                  r = 2)
other_spc <- paste(other_spc[,1], other_spc[,2], sep = "_vs_")
# Designate species of interest as query genome in all combinations
species_versus <- c(paste(spc_int,
                          grep(spc_int, species, invert = T, value = T),
                          sep = "_vs_"),
                    other_spc)
    # # Keep only unique combinations of remaining species
    # unlist(lapply(combn(species[!species %in% spc_int], 2, simplify = F),
    #               paste, collapse = "_vs_")))
# Create list of PSL filenames
psl_files <- sapply(species_versus, grep, list.files(pattern = "\\.psl",
                                                     path = ".",
                                                     full.names = T), value = T)
psl_files <- grep("new", psl_files, value = T, invert = T)
# Read PSL files into list of dataframes
psl_list <- sapply(psl_files, read.table, col.names = psl_col,
                   simplify = F, USE.NAMES = T)
# Remove file extensions from list names
names(psl_list) <- gsub(".*ment_|.psl", "", names(psl_list))
# Convert scaffold names to size-ordered factors
psl_list <- sapply(names(psl_list), orderPsl, psl_list, simplify = F)

# Analysis
# Summarize by contig vs. contig of each syntenic comparison
psl_sums <- sapply(names(psl_list), sumPsl, psl_list, simplify = F)
# Select maximal matching contigs
psl_match <- sapply(names(psl_sums), maxMatches, psl_sums, simplify = F)
# Compare matches
match_lens <- collapseMatch(psl_match)
plotLens(match_lens)
# match_lens <- sapply(grep("latissima", names(psl_match), value = T), lenvsMatch,
#                     psl_match, psl_list, simplify = F)
# match_lens <- match_lens[sapply(spc_order, grep, names(match_lens), value = T)]
# len_plots <- sapply(names(match_lens), plotLens, match_lens, simplify = F)
# ((len_fig <- ggarrange(plotlist = len_plots)))
# showtext_opts(dpi = 300)
# ggsave(plot = len_fig, filename = "homolog_lengths.png", path = outdir, bg = "white",
#        width = 12, height = 9, units = "in")
# showtext_opts(dpi = 100)

# # # Rearrange scaffold ID factors by matches for plotting
# # psl_list <- sapply(names(psl_match), plotOrder, psl_match, psl_list,
# #                    simplify = F)
# # Pivot summarized data for heatmaps
# psl_pivs <- sapply(names(psl_sums), pivotPsl, psl_sums, simplify = F)
# # # Group scaffolds
# # test <- groupScafs(psl_match, spc_int)
# 
# # Plots
# # Heatmaps of genome vs. genome synteny
# psl_heats <- sapply(names(psl_match), heatPsl, psl_match, psl_pivs,
#                     simplify = F)
# # p_heat <- ggarrange(plotlist = psl_heats, legend = "right")
# ((p2 <- ggarrange(psl_heats[["Saccharina_latissima_vs_Macrocystis_pyrifera"]],
#                 psl_heats[["Saccharina_japonica_vs_Macrocystis_pyrifera"]],
#                 psl_heats[["Saccharina_latissima_vs_Saccharina_japonica"]],
#                 psl_heats[["Macrocystis_pyrifera_vs_Saccharina_japonica"]],
#                 legend = "right")))
# h_fnames <- unlist(sapply(names(psl_heats), plotSave, "heatmap", psl_heats,
#                           outdir, 7, 10, simplify = F))
# print(unname(h_fnames))
# 
# # Dotplots of genome vs. genome
# psl_dots <- sapply(names(psl_list), dotPlot, psl_list, simplify = F)
# d_fnames <- unlist(sapply(names(psl_dots), plotSave, "dotplot", psl_dots, outdir,
#                    15, 10, simplify = F))
# print(unname(d_fnames))


# # Separate out certain syntenic regions
# seps <- sapply(names(psl_list), sepDots, psl_list, simplify = F)
# sep_list <- sapply(names(seps), plotList, seps, simplify = F)
# sapply(names(sep_list), plotSave, sep_list, 7, 7, simplify = F)


# # Calculate scaffold-scaffold alignment metrics in each PSL dataframe
# psl_metrics <- sapply(names(psl_list), metricsPsl, psl_list, simplify = F)
# test <- psl_metrics[["Saccharina_latissima_vs_Macrocystis_pyrifera"]]
# test$euc_dist <- (1-test$cov_scaled)^2+(0+test$total_p_blockCount_scaled)^2



