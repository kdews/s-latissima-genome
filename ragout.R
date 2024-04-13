# Clear environment
rm(list = ls())
# Required packages
suppressPackageStartupMessages(library(tidyverse, quietly = T, warn.conflicts = F))
# library(ggrepel)
suppressPackageStartupMessages(library(ggpmisc, quietly = T, warn.conflicts = F))
library(ape, quietly = T, warn.conflicts = F)
library(viridisLite, quietly = T, warn.conflicts = F)
library(RColorBrewer, quietly = T, warn.conflicts = F)
library(gridExtra, quietly = T, warn.conflicts = F)
library(ggpubr, quietly = T, warn.conflicts = F)
if (require(showtext, quietly = T, warn.conflicts = F)) {
  showtext_auto()
  if (interactive()) showtext_opts(dpi = 100) else showtext_opts(dpi = 300)
}

# Input
if (interactive()) {
  wd <- "/scratch1/kdeweese/latissima/genome_stats"
  setwd(wd)
  # Table of species in analysis
  seqFile <- "s-latissima-genome/s_lat_alignment.txt"
  # Ragout output directory
  ragout_dir <- "ragout-out-all-solid-refine-update"
  # Species of interest
  spc_int <- "Saccharina_latissima"
  # Output directory
  outdir <- "s-latissima-genome/"
} else {
  line_args <- commandArgs(trailingOnly = T)
  # Table of species in analysis
  seqFile <- line_args[1]
  # Ragout output directory
  ragout_dir <- line_args[2]
  # Species of interest
  spc_int <- line_args[3]
  # Output directory
  outdir <- line_args[4]
}

# Output plot file names
outfiles <- list(
  len_plot = "ragout_length_dist.png",
  map_plot = "ragout_pseudo_mapping.png",
  comp_len_plot = "ragout_comp_length_dist.png",
  comp_map_plot = "ragout_comp_pseudo_mapping.png"
)
# If exists, prepend output directory to output file names
redirect <- function(filename, outdir) {
  if (dir.exists(outdir)) filename <- paste0(outdir, filename)
  return(filename)
}
outfiles <- sapply(outfiles, redirect, outdir, simplify = F)

# Global variables
types <- c("input", "pseudochromosomes", "chimera", "remainder", "excluded")
label_names <- c("input", "gap")

# Functions
# Checks existence of file or directory and errors if FALSE
checkPath <- function(test_path) {
  if (file.exists(test_path)) {
    return(TRUE)
  } else {
    stop(paste0("Cannot follow path (", test_path, ")."))
  }
}
# Fixes chromosome labels for plotting
fixChrom <- function(contigs) {
  contigs <-
    as.character(as.numeric(str_remove_all(str_remove_all(contigs, ".*_"),
                                           "[^0-9]")))
  return(contigs)
}
# Function to extract title from directory name
extTtl <- function(ragout_dir) {
  ttl_dict <- list(all="No filtering,",
                   `no-sac`="Short contigs filtered out (except S. latissima),",
                   filt="Short contigs filtered out,",
                   solid="unbroken scaffolds",
                   chimera="broken scaffolds")
  ttl <- gsub("ragout-out-", "", ragout_dir)
  vals <- list()
  for (val in names(ttl_dict)) if (grepl(val, ttl)) vals <- append(vals, val)
  for (val in vals) ttl <- gsub(val, ttl_dict[[val]], ttl)
  ttl <- gsub("-refine", "", ttl)
  ttl <- gsub("-", " ", ttl)
  return(ttl)
}
# Read Cactus-formatted seqFile
readSpecies <- function(seqFile) {
  seq_cols <- c("Species", "Assembly")
  checkPath(seqFile)
  seqs <- read.table(seqFile, sep = "\t", fill = NA, header = F,
                     # Omit phylogenetic tree line
                     comment.char = "(",
                     col.names = seq_cols)
  colnames(seqs) <- seq_cols[1:length(seqs)]
  seqs <- seqs[,na.omit(colnames(seqs))]
  return(seqs)
}
# Import NCBI AGP file (v2.0)
readAgp <- function(ragout_dir) {
  base_cols <- c("seqid_scaf", "start_scaf", "end_scaf",
                 "component_number", "component_type")
  agp_cols <- c(base_cols,
                "seqid_comp", "start_comp", "end_comp", "orientation")
  gap_cols <- c(base_cols,
                "gap_length", "gap_type", "linkage", "evidence")
  names(agp_cols) <- gap_cols
  # AGP file describing scaffolding
  checkPath(ragout_dir)
  agp_file <- list.files(path = ragout_dir, pattern = ".*agp", full.names = T)
  agp_raw <- read.table(agp_file, col.names = agp_cols)
  agp <- agp_raw %>%
    filter(component_type == "W") %>%
    select(!contains("component")) %>%
    mutate(label = factor("input", levels = label_names))
  gaps_df <- agp_raw %>%
    rename(all_of(agp_cols)) %>%
    filter(component_type == "N") %>%
    select(!contains(c("component", "linkage", "type"))) %>%
    # Quantify number of species used as evidence for gap
    rowwise() %>%
    mutate(n_evidence = length(unlist(strsplit(evidence, ",")))) %>%
    ungroup() %>%
    mutate(label = factor("gap", levels = label_names))
  agp_combo <- merge(agp, gaps_df, all = T)
  agp_combo <- agp_combo %>%
    # Convert genomic position columns to numeric type and bp to Mb
    mutate_at(.vars = vars(grep("start|end|length", colnames(.), value = T)),
              .funs = ~ as.numeric(.x)*10^-6) %>%
    group_by(seqid_scaf) %>%
    # Sum N content of each scaffold
    mutate(N_scaf = sum(gap_length, na.rm = T)) %>%
    ungroup()
  return(agp_combo)
}
# Filters out gaps from AGP dataframe
filtAgp <- function(agp) {
  # Calculate N content and total length before filtering
  agp_filt <- agp %>%
    mutate(N_count = sum(as.numeric(grep("scaffold_", seqid_comp,
                                         value = T, invert = T)), na.rm = T),
           total_len = sum(abs(as.numeric(start_scaf)-as.numeric(end_scaf)+1)),
           perc_N = round(N_count/total_len, digits = 4)) %>%
    filter(component_type == "W") %>%
    select(!contains("component")) %>%
    # Convert genomic position columns to numeric type
    mutate_at(.vars = vars(grep("start|end", colnames(.), value = T)),
              .funs = as.numeric)
  return(agp_filt)
}
# Import FASTA index file (.fasta.fai)
readFai <- function(idx_file, suff = NULL) {
  checkPath(idx_file)
  idx <- read.table(idx_file)
  idx <- idx[,1:2]
  colnames(idx) <- c("seqid", "length")
  idx <- idx %>%
    arrange(desc(length)) %>%
    mutate(seqid = factor(seqid, levels = seqid))
  if (!missing(suff)) {
    colnames(idx) <- paste0(colnames(idx), suff)
  }
  return(idx)
}
# Generate annotated index of pre- and post-Ragout scaffolds
genIdx <- function(ragout_dir, agp_list, seqs) {
  agp <- agp_list[[ragout_dir]] %>% filter(label == "input")
  # Pre-Ragout scaffold index (later filtered for plotting)
  pre_file <- paste0(pull(filter(seqs, grepl(spc_int, Species)),
                           Assembly), ".fai")
  pre_idx <- readFai(pre_file)
  # Ragout rescaffolded scaffold index
  scaf_file <- list.files(path = ragout_dir,
                          pattern = "_scaffolds\\.fasta\\.fai",
                          full.names = T)
  scaf <- readFai(scaf_file)
  # Unplaced scaffold index
  unplc_file <- list.files(path = ragout_dir,
                           pattern = "_unplaced\\.fasta\\.fai",
                           full.names = T)
  unplc <- readFai(unplc_file)
  # Create vector of chimeric scaffolds
  chims <- unique(gsub("\\[.*\\]", "", grep("\\[", unplc$seqid, value = T)))
  # Create vector of rescaffolded and unplaced seqids
  rag_seqids <- c(agp$seqid_comp, gsub("\\[.*\\]", "", unplc$seqid))
  # Label scaffolds with "type" column before combining
  # Type: rescaffolded
  scaf <- scaf %>% mutate(type = "pseudochromosomes")
  # Type: remainder
  unplc <- unplc %>% mutate(type = "remainder") %>%
    # Filter out chimeric scaffolds from remainder
    rowwise() %>%
    filter(!grepl("\\[", seqid)) %>%
    ungroup()
  pre_idx <- pre_idx %>%
    # Type: chimera
    mutate(type = case_when(seqid %in% chims ~ "chimera",
                            # Type: input
                            (seqid %in% rag_seqids &
                               !seqid %in% chims) ~ "input",
                            # Type: excluded
                            .default = "excluded"))
  # Combine all scaffold types
  idx <- rbind(scaf, unplc, pre_idx) %>%
    mutate(type = factor(type, levels = types),
           seqid = factor(seqid, levels = unique(as.character(seqid))),
           length = log10(length)) %>%
    arrange(desc(length))
  return(idx)
}
# Annotate AGP dataframe with component synteny to chr0
annotChr0 <- function(idx_agp) {
  comp_chr0 <- read.table("comp_chr0.txt", header = T, sep = "\t")
  comp_chr0 <- comp_chr0 %>%
    mutate(qName = paste0("scaffold_", qNum))
  dict_chr0 <- comp_chr0$qPercent
  names(dict_chr0) <- comp_chr0$qName
  idx_agp <- idx_agp %>%
    mutate(chr0 = case_when(seqid_comp %in% names(dict_chr0) ~
                              dict_chr0[as.character(seqid_comp)],
                            .default = 0))
  return(idx_agp)
}
# Index AGP dataframe with scaffold/contig lengths
idxAgp <- function(ragout_dir, agp_list, idx_list) {
  agp <- agp_list[[ragout_dir]]
  idx <- idx_list[[ragout_dir]] %>% filter(type == "pseudochromosomes")
  pre_idx <- idx_list[[ragout_dir]] %>%
    filter(type %in% c("input", "chimera")) %>%
    arrange(desc(length)) %>%
    mutate(seqid = factor(seqid, levels = unique(as.character(seqid))))
  idx_dict <- idx$length
  names(idx_dict) <- idx$seqid
  pre_idx_dict <- pre_idx$length
  names(pre_idx_dict) <- pre_idx$seqid
  idx_agp <- agp %>%
    mutate(
      # Order component seqid factors by length (for plotting)
      seqid_comp = factor(seqid_comp, levels = levels(pre_idx$seqid)),
      # Create length columns
      length_scaf = as.numeric(idx_dict[seqid_scaf]),
      length_comp = as.numeric(pre_idx_dict[seqid_comp]),
      # Convert length out of log10(bp) scale to Mb scale
      length_scaf = (10^length_scaf)*10^-6,
      length_comp = (10^length_comp)*10^-6
           ) %>%
    arrange(desc(length_scaf), start_scaf)
  idx_agp <- annotChr0(idx_agp)
  return(idx_agp)
}
# Takes factor vector and returns factor-named vector of colors (i.e., dict)
getColors <- function(vec, palFun, pal = NULL) {
  if (!is.factor(vec)) stop("Error: vector class is not factor.")
  if (identical(palFun, brewer.pal)) {
    if (missing(pal)) pal = "Paired"
    vec_colors <- palFun(length(levels(vec)), pal)
    names(vec_colors) <- levels(vec)
    # Color chimera red using "Paired" palette
    if (any(names(vec_colors) %in% c("chimera"))) {
      new_ord <- c("input", "pseudochromosomes", "remainder", "excluded", "chimera")
      if (any(names(vec_colors) %in% c("total"))) new_ord <- c(new_ord, "total")
      names(vec_colors) <- new_ord
    }
  } else if (identical(palFun, viridis)) {
    if (missing(pal)) pal = "turbo"
    vec_colors <- viridis(n = length(levels(vec)), option = pal)
    names(vec_colors) <- levels(vec)
  }
  return(vec_colors)
}
# Get common legend for variable in list of dataframes using a plotting function
getCommon <- function(plotFun, var_name, df_list, ...) {
  var_count <- sapply(df_list,
                      function(x) length(unique(pull(x, matches(var_name)))),
                      simplify = F)
  legend_name <- names(which.max(var_count))
  common_legend <- get_legend(plotFun(legend_name, df_list, ...))
  return(common_legend)
}
# Barplot of scaffolds indexed by length
idxPlot <- function(ragout_dir, idx_list, ttls = NULL, pal = "Paired",
                    exclude = NULL) {
  x_max <- max(unlist(sapply(idx_list,
                             function(x) length(unique(pull(x, seqid))),
                             simplify = F)), na.rm = T)
  y_max <- max(unlist(sapply(idx_list, pull, length, simplify = F)), na.rm = T)
  common_lims <- lims(x = c(0, x_max + 1), y = c(0, y_max + 1))
  idx <- idx_list[[ragout_dir]]
  ttl <- ttls[[ragout_dir]]
  if(!missing(exclude)) idx <- idx %>% filter(!type %in% exclude)
  # Filter out input scaffolds that are chimera, remainder, or excluded
  if (any(idx$type %in% c("input"))) {
    filt_seqids <- idx %>% filter(type != "input") %>% pull(seqid) %>%
      gsub("\\[.*\\]", "", .)
    idx <- idx %>% filter(!(type == "input" & seqid %in% filt_seqids))
  }
  # Colors for plot
  vec_colors <- getColors(idx$type, brewer.pal, pal)
  idx_sum <- idx %>%
    group_by(type) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    # Add total n row
    mutate(type = factor(type, levels = c(types, "total"))) %>%
    rbind(list(factor("total", levels = c(types, "total")), sum(.$n)))
  # Colors for table
  tbl_colors <- getColors(idx$type, brewer.pal, pal)
  tbl_colors[["total"]] <- "white"
  fnt_colors <- rep("black", length(tbl_colors))
  names(fnt_colors) <- names(tbl_colors)
  fnt_colors[names(fnt_colors) %in% c("pseudochromosomes", "excluded")] <- "white"
  tbl_colors <- tbl_colors[as.character(unique(idx_sum$type))]
  fnt_colors <- fnt_colors[as.character(unique(idx_sum$type))]
  annot <- annotate(geom = "table", x = x_max, y = y_max,
                    label = idx_sum, table.colnames = F,
                    table.theme =
                      ttheme_default(padding = unit(c(2, 2), "mm"),
                                     base_size = 10,
                                     core = list(bg_params =
                                                   list(fill = tbl_colors),
                                                 fg_params =
                                                   list(col = fnt_colors))))
  idx <- idx %>% arrange(type)
  p <- ggplot(data = idx, mapping = aes(fill = type)) +
    geom_col(mapping = aes(x = as.numeric(row.names(idx)), y = length),
             width = 1) +
    common_lims +
    scale_fill_manual(values = vec_colors) +
    labs(title = str_wrap(ttl, width = 35),
         x = "Scaffold index", y = "Scaffold length (log10 bp)") +
    theme_minimal() +
    annot +
    theme(legend.position = "top",
          legend.title = element_blank())
  return(p)
}
# Combine pre- and post-Ragout index plots
combPlot <- function(ragout_dir, idx_list, ttls = NULL, pal = "Paired",
                     legend = F) {
  p1 <- idxPlot(ragout_dir, idx_list, ttls = ttls, pal = pal,
                exclude = c("pseudochromosomes")) +
    labs(subtitle = "Before", x = "", y = "")
  p2 <- idxPlot(ragout_dir, idx_list, pal = pal,
                exclude = c("input")) +
    labs(subtitle = "After", x = "", y = "")
  p_list <- list(p1, p2)
  if (legend) {
    common_legend <- getCommon(idxPlot, "type", idx_list)
    p <- ggarrange(plotlist = p_list, nrow = length(p_list),
                   align = "hv", legend.grob = common_legend, legend = "top")
    p <- annotate_figure(p,
                         bottom = "Scaffold index",
                         left = "Scaffold length (log10 bp)")
  } else {
    p <- ggarrange(plotlist = p_list, nrow = length(p_list),
                   align = "hv", legend = "none")
  }
  return(p)
}
# Plot of ragout rescaffolding
ragoutPlot <- function(ragout_dir, idx_agp_list, ttls = NULL, pal = "Paired",
                       labels = T) {
  ttl <- ttls[[ragout_dir]]
  idx_agp <- idx_agp_list[[ragout_dir]]
  # Common x-axis limits for aligning all plots
  x_max <- max(unlist(sapply(idx_agp_list, pull, length_scaf)), na.rm = T)*1.1
  common_x_lims <- xlim(c(0, x_max + 1))
  # # Color each seqid with consistent colors in all graphs
  # all_seqids <- lapply(idx_agp_list, select, c(seqid_comp, length_comp)) %>%
  #   bind_rows() %>%
  #   unique() %>%
  #   arrange(desc(length_comp)) %>%
  #   mutate(seqid_comp = as.numeric(fixChrom(seqid_comp))) %>%
  #   pull(seqid_comp)
  # all_seqids <- factor(all_seqids, levels = all_seqids)
  # myColors <- getColors(all_seqids, viridis, "turbo")
  # Factor seqids by scaffold length and start position (for plotting)
  order_scaf <- unique(pull(idx_agp, seqid_scaf))
  order_comp <- fixChrom(levels(pull(idx_agp, seqid_comp)))
  # Arrange data by scaffold length, then start position
  idx_agp <- idx_agp %>%
    mutate(seqid_comp = fixChrom(seqid_comp),
           seqid_comp = factor(seqid_comp, levels = order_comp),
           seqid_scaf = as.factor(as.numeric(factor(seqid_scaf,
                                                    levels = order_scaf))),
           perc_N = (N_scaf/length_scaf)*100)
  # Summary of data
  # Positions
  xpos <- x_max*0.8
  ypos <- length(levels(idx_agp$seqid_scaf))
  sum_df <- idx_agp %>%
    select(seqid_scaf, N_scaf, length_scaf, perc_N) %>% unique()
  tot_len <- sum_df %>% pull("length_scaf") %>% sum()
  tot_N <- sum_df %>% pull("N_scaf") %>% sum()
  tot_perc_N <- (tot_N/tot_len)*100
  sum_df <- sum_df %>%
    mutate(perc_N = paste0(round(perc_N, digits = 2), "%"))
  # Color by labels with consistent colors in all graphs
  myColors <- getColors(idx_agp$label, brewer.pal, pal)
  # myColors <- viridis(4, option = "rocket")
  p <- ggplot(data = idx_agp) +
    geom_segment(mapping = aes(x = start_scaf, xend = end_scaf,
                               y = seqid_scaf, yend = seqid_scaf,
                               # col = perc_N,
                               col = label,
                               linewidth = label)) +
    annotate(geom = "text", x = sum_df$length_scaf, y = sum_df$seqid_scaf,
             label = sum_df$perc_N, hjust = -0.3) +
    annotate(geom = "text", x = xpos, y = levels(idx_agp$seqid_scaf)[ypos],
             label = paste(round(tot_len, digits = 2), "Mb")) +
    annotate(geom = "text", x = xpos, y = levels(idx_agp$seqid_scaf)[ypos-1],
             label = paste0(round(tot_perc_N, digits = 2), "% N's")) +
    labs(title = str_wrap(ttl, width = 35),
         subtitle = "Scaffold mapping onto pseudochromosomes",
         x = "Scaffold length (Mb)", y = "Pseudochromosome index") +
    theme_classic() +
    scale_color_manual(values = myColors) +
    # scale_color_gradientn(colors = myColors,
    #                       # values = c(0, 0.1, 0.25, 1),
    #                       limits = c(0, 100)) +
    # scale_color_gradient2(low = myColors[1],
    #                       mid = myColors[2],
    #                       high = myColors[3],
    #                       midpoint = 5, limits = c(0, 100)) +
    scale_linewidth_discrete(range = c(3, 1)) +
    theme(legend.position = "top") +
    common_x_lims
  if (!labels) {
    p <- p + labs(x = "", y = "")
  }
  return(p)
}
# Run all functions on data
runAnalysis <- function(ragout_dirs, seqs, plot1, plot2) {
  # Extract titles for plots
  ttls <- sapply(ragout_dirs, extTtl)
  # Wrangle data
  # Import AGP
  agp_list <- sapply(ragout_dirs, readAgp, simplify = F)
  # Import FASTA indices
  idx_list <- sapply(ragout_dirs, genIdx, agp_list, seqs, simplify = F)
  # Index AGP
  idx_agp_list <- sapply(ragout_dirs, idxAgp, agp_list, idx_list, simplify = F)
  # Plots
  # Bar plot length distribution before and after rescaffolding
  p_list <- sapply(ragout_dirs, combPlot, idx_list,
                   ttls = ttls,
                   simplify = F)
  common_legend <- getCommon(idxPlot, "type", idx_list)
  all_bar <- ggarrange(plotlist = p_list,
                       legend.grob = common_legend, legend = "top",
                       align = "hv", ncol = length(p_list))
  all_bar <- annotate_figure(all_bar,
                             bottom = "Scaffold index",
                             left = "Scaffold length (log10 bp)")
  # Save plot
  ht <- 7
  wd <- 7
  if (length(p_list) > 1) wd <- 7*length(p_list)*0.6
  ggsave(plot1, all_bar,
         height = ht,
         width = wd,
         bg = "white")
  print(paste("Saved plot:", plot1))
  # Line graph mapping of original scaffolds onto pseudochromosomes
  dot_list <- sapply(ragout_dirs, ragoutPlot, idx_agp_list, labels = F,
                     ttls = ttls,
                     simplify = F)
  common_legend <- getCommon(ragoutPlot, "label", idx_agp_list)
  all_dot <- ggarrange(plotlist = dot_list,
                       legend.grob = common_legend, legend = "top",
                       align = "hv", ncol = length(dot_list))
  all_dot <- annotate_figure(all_dot,
                             top =
                               text_grob("Reference-based scaffold ordering",
                                         face = "bold", size = 14),
                             bottom = "Scaffold length (Mb)",
                             left = "Scaffold index")
  # Save plot2
  ht <- 7
  wd <- 7
  if (length(dot_list) > 1) wd <- 7*length(dot_list)*0.6
  ggsave(plot2, all_dot,
         height = ht,
         width = wd,
         bg = "white")
  print(paste("Saved plot:", plot2))
  # # Combine all plots into figure
  # comp_rag <- ggarrange(all_bar, all_dot, align = "hv", nrow = 2)
  # comp_rag <- annotate_figure(comp_rag,
  #                             top =
  #                               text_grob("Reference-based scaffold ordering",
  #                                         face = "bold", size = 14))
  
}


# Import data
seqs <- readSpecies(seqFile)
ragout_dirs <- list.files(pattern = "ragout-out-")
ragout_dirs <- grep("refine|filt", ragout_dirs, value = T)
suppressWarnings(
  runAnalysis(c(ragout_dir), seqs, outfiles$len_plot, outfiles$map_plot)
)
suppressWarnings(
  runAnalysis(ragout_dirs, seqs, outfiles$comp_len_plot, outfiles$comp_map_plot)
)

