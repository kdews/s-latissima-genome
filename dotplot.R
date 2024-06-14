# Clear environment
rm(list = ls())
# Required packages
library(Hmisc, quietly = T, warn.conflicts = F)
library(scales, quietly = T, warn.conflicts = F)
suppressPackageStartupMessages(library(tidyverse, quietly = T, warn.conflicts = F))
library(ggpubr, quietly = T, warn.conflicts = F)
suppressPackageStartupMessages(library(ggpmisc, quietly = T, warn.conflicts = F))
library(RColorBrewer, quietly = T)
library(BiocManager, quietly = T)
suppressPackageStartupMessages(library(ComplexHeatmap, quietly = T))
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
# Output files
max_match_file <- "max_matches.tsv"
align_report_file <- "alignment_report.tsv"
# Prepend output directory to file name (if it exists)
if (dir.exists(outdir)) align_report_file <- paste0(outdir, align_report_file)


# Functions
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
# Order given ID column by size column in data frame
sizeSort <- function(df, id_col, size_col) {
  x <- df %>% select(all_of(c(id_col, size_col))) %>%
    arrange(desc(pick(all_of(size_col)))) %>%
    pull(all_of(id_col)) %>%
    unique()
  return(x)
}
# Order PSL data frame by contig size by converting contig names to factors
orderPsl <- function(df_name, df_list) {
  df <- df_list[[df_name]]
  query <- getGen(df_name, "query")
  target <- getGen(df_name, "target")
  # Extract numeric information from IDs
  df <- df %>% mutate(qNum=fixChrom(qName), tNum=fixChrom(tName))
  # Order ID character vectors by size
  query_contigs <- sizeSort(df, "qName", "qSize")
  target_contigs <- sizeSort(df, "tName", "tSize")
  query_nums <- sizeSort(df, "qNum", "qSize")
  target_nums <- sizeSort(df, "tNum", "tSize")
  # Capture size rank of each target ID with index
  target_idx <- df %>% select(tNum, tSize) %>%
    arrange(desc(tSize)) %>%
    unique() %>%
    # Filter out artificial chromosomes for ranking
    filter(tNum != "0") %>%
    # Use row ids for index
    rowid_to_column(var="index") %>%
    select(tNum, index) %>%
    # Force any artificial chromosomes to 0 index
    rbind(data.frame(index = 0, tNum = "0")) %>%
    # Convert to named vector
    deframe()
  df <- df %>% mutate(
    # Add integer index column using named vector
    index=as.integer(target_idx[tNum]),
    # Convert all IDs to factors ordered by size
    qName=factor(qName, levels = query_contigs),
    tName=factor(tName, levels = target_contigs),
    qNum=factor(qNum, levels = query_nums),
    tNum=factor(tNum, levels = target_nums)
  ) %>%
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
    filter(tNum != "0") %>%
    arrange(Species)
  write.table(match_lens, file = max_match_file,
              quote = F, row.names = F, sep = "\t")
  print(paste("Table of maximal matches in:", max_match_file))
  return(match_lens)
}
# Summarize collapsed match_lens data frame
sumMatch <- function(match_lens) {
  match_lens_sum <- match_lens %>% group_by(index, tNum, tSize, Species) %>%
    # Sums per ID (i.e., per contig)
    summarize(sum_homolog=sum(qSize),
              sum_match=sum(total_matches),
              `n homologs`=n(),
              .groups = "keep") %>% 
    ungroup() %>% arrange(Species, desc(tSize)) 
  return(match_lens_sum)
}
# Report alignment statistics in table
alnReport <- function(match_lens_sum) {
  align_report <- match_lens_sum %>% 
    # Convert bp to Mb
    mutate_at(grep("sum|size", colnames(.), ignore.case = T, value = T),
              ~.*1e-6) %>%
    # Per species statistics
    group_by(Species) %>%
    summarize(`Average homologs mapped per chromosome`=
                paste(round(smean.sd(`n homologs`), 2), collapse = " ± "),
              `Maximum homologs mapped per chromosome`=max(`n homologs`),
              `Minimum homologs mapped per chromosome`=min(`n homologs`),
              `Total homologs mapped`=sum(`n homologs`),
              `Average exact matches (Mb) per chromosome`=
                paste(round(smean.sd(sum_match), 2), collapse = " ± "),
              `Average exact matches (%) per chromosome`=
                paste(round(smean.sd(sum_match/tSize*100), 2), collapse = " ± "),
              `Maximum exact matches (Mb) per chromosome`=max(sum_match),
              `Minimum exact matches (Mb) per chromosome`=min(sum_match),
              `Total exact matches (Mb)`=sum(sum_match),
              .groups = "keep") %>%
    # Coerce all values to characters
    ungroup() %>% mutate(across(everything(), as.character)) %>%
    # Transpose table
    pivot_longer(cols = !Species,
                 names_to = " ", values_to = "Value") %>%
    pivot_wider(names_from = "Species", values_from = "Value")
  write.table(align_report, file = align_report_file,
              quote = F, row.names = F, sep = "\t")
  print(paste("Table of alignment statistics in:", align_report_file))
  return(align_report)
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
    mutate(qNum = factor(x = qNum, levels = h_order, ordered = T),
           qindex = as.integer(qNum))
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
# Plot number and length of scaffold matches vs. chromosome length (match_lens)
plotLens <- function(match_lens_sum) {
  # Clean up data frame for plotting
  match_lens_sum <- match_lens_sum %>%
    # Convert bp to Mb
    mutate(`Reference chromosome length (Mb)`=tSize*1e-6,
           `summed homolog lengths (Mb)`=sum_homolog*1e-6,
           `summed match lengths (Mb)`=sum_match*1e-6,
           # Change Species variable name for legend
           Reference=Species)
  # Function for individual plots
  lenPlot <- function(my_var, lab_pos) {
    p <- ggplot(data = match_lens_sum,
                mapping = aes(x = `Reference chromosome length (Mb)`,
                              y = .data[[my_var]], group = Reference,
                              col = Reference, fill = Reference)) +
      geom_point(alpha = 0.5) +
      stat_poly_line(formula = y~x+0, alpha = 0.1) +
      stat_poly_eq(formula = y~x+0, mapping = use_label(labels = c("R2", "eq")),
                   label.x = lab_pos[["x"]], label.y = lab_pos[["y"]]) +
      scale_x_continuous(breaks = pretty_breaks()) +
      scale_y_continuous(breaks = pretty_breaks()) +
      theme_bw() +
      theme(legend.text = element_text(face = "italic"))
    return(p)
  }
  # Function to annotate left of plot with italicized name of species of interest
  annotSpcInt <- function(p) {
    p <- annotate_figure(
      p,
      left = text_grob(abbrevSpc(spc_int), face = "italic", rot = 90)
    )
    return(p)
  }
  p1 <- lenPlot("summed homolog lengths (Mb)", c(x = "right", y = "bottom"))
  leg1 <- get_legend(p1)
  p1 <- annotSpcInt(p1 + theme(legend.position = "none"))
  p2 <- lenPlot("summed match lengths (Mb)", c(x = "right", y = "top"))
  p2 <- annotSpcInt(p2 + theme(legend.position = "none"))
  p <- ggarrange(p1, p2, nrow = 2, align = "hv",
                 legend.grob = leg1, legend = "right", labels = "AUTO")
  p <- setNames(list(p), spc_int)
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
  leg_name <- paste("Synteny", "coverage", sep = "\n")
  # Filter out artificial chromosomes
  match <- match %>% filter(tNum != "0", qNum != "0")
    # # Filter out small scaffolds/contigs from query (remove <1Mb)
    # filter(qSize > 1e6)
  # Order of query contigs while their respective homologs (tNum) are size-ordered
  h_order <- match %>% arrange(tNum) %>% pull(qNum) %>% as.character()
  # Apply filtering from match data frame
  df <- df %>% filter(tNum %in% match$tNum, qNum %in% match$qNum) %>%
    # Use "match" data frame to order query contig factors in heatmap
    mutate(qNum = factor(x = qNum, levels = h_order, ordered = T),
           qindex = as.integer(qNum))
  if (query == spc_int) {
    p <- ggplot(data = df,
                mapping = aes(fill = qPercent, x = tNum, y = qindex)) +
      scale_y_continuous(expand = c(0, 0), breaks = breaks_width(100))
  } else {
    p <- ggplot(data = df, mapping = aes(fill = qPercent, x = tNum, y = qNum))
    }
  p <- p + geom_tile() +
    scale_fill_viridis_c(option = "turbo", labels = label_percent(),
                         name = leg_name) +
    labs(x = x_label, y = y_label) +
    theme(axis.title = element_text(face="italic"),
          legend.position = "left")
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
    df <- df %>% filter(qSize > 1e6, qNum != "0", tNum != "0")
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
  query <- getGen(df_name, "query")
  x_label <- abbrevSpc(target)
  y_label <- abbrevSpc(query)
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
    labs(x = "", y = "")
  if (missing(filter_ids)) {
    p <- p +
      labs(title = "Genome synteny by scaffold", x = x_label, y = y_label) +
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
  if (dir.exists(outdir)) fname <- paste0(outdir, fname)
  ggsave(filename = fname, plot = p, bg = "white",
         width = width, height = height, units = "in")
  print(fname)
  return(fname)
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
# Select out PSLs with new filtering
psl_files <- grep("new", psl_files, value = T, invert = T)
# Read PSL files into list of dataframes
psl_list_raw <- sapply(psl_files, read.table, col.names = psl_col, simplify = F)
# Remove file extensions from list names
names(psl_list_raw) <- gsub(".*ment_|.psl", "", names(psl_list_raw))
# Convert scaffold names to size-ordered factors
psl_list <- sapply(names(psl_list_raw), orderPsl, psl_list_raw, simplify = F)

# Analysis
# Summarize by contig vs. contig of each syntenic comparison
psl_sums <- sapply(names(psl_list), sumPsl, psl_list, simplify = F)
# Select maximal matching contigs
psl_match <- sapply(names(psl_sums), maxMatches, psl_sums, simplify = F)
# Compare alignment statistics between genomes and species of interest
match_lens <- collapseMatch(psl_match)
# Summarize alignment statistics per reference contig
match_lens_sum <- sumMatch(match_lens)
# Save report of alignment statistics by reference species used
align_report <- alnReport(match_lens_sum)
# Pivot summarized data for heatmaps
psl_pivs <- sapply(names(psl_sums), pivotPsl, psl_sums, simplify = F)

# Plots
# Alignment statistics against species of interest
len_fig <- plotLens(match_lens_sum)
showtext_opts(dpi = 300)
homolog_lens_plot <- plotSave(spc_int, "homolog_lengths", len_fig, outdir,
                              6.5, 9.5)
# Heatmaps of genome vs. genome synteny
psl_heats <- sapply(names(psl_match), heatPsl, psl_match, psl_pivs,
                    simplify = F)
print("Saving heatmaps...")
h_fnames <- unlist(sapply(names(psl_heats), plotSave, "heatmap", psl_heats,
                          outdir, 9, 9, simplify = F))
# Dotplots of genome vs. genome
# Rearrange scaffold ID factors by matches for plotting
psl_list <- sapply(names(psl_match), plotOrder, psl_match, psl_list,
                   simplify = F)
dotPlot(names(psl_list)[[4]], psl_list)
print("Saving dotplots...")
psl_dots <- sapply(names(psl_list), dotPlot, psl_list, simplify = F)
d_fnames <- unlist(sapply(names(psl_dots), plotSave, "dotplot", psl_dots, outdir,
                   10, 10, simplify = F))










# # Separate out certain syntenic regions
# seps <- sapply(names(psl_list), sepDots, psl_list, simplify = F)
# sep_list <- sapply(names(seps), plotList, seps, simplify = F)
# sapply(names(sep_list), plotSave, sep_list, 7, 7, simplify = F)


# # Calculate scaffold-scaffold alignment metrics in each PSL dataframe
# psl_metrics <- sapply(names(psl_list), metricsPsl, psl_list, simplify = F)
# test <- psl_metrics[["Saccharina_latissima_vs_Macrocystis_pyrifera"]]
# test$euc_dist <- (1-test$cov_scaled)^2+(0+test$total_p_blockCount_scaled)^2



