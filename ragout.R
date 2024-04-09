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
  assembly_file <- "s-latissima-genome/species_table.txt"
  # Ragout output directory
  ragout_dir <- "ragout-out-new-cactus"
  # Species of interest
  spc_int <- "Saccharina_latissima"
  # Output directory
  outdir <- "s-latissima-genome/"
} else {
  line_args <- commandArgs(trailingOnly = T)
  # Table of species in analysis
  assembly_file <- line_args[1]
  # Ragout output directory
  ragout_dir <- line_args[2]
  # Species of interest
  spc_int <- line_args[3]
  # Output directory
  outdir <- line_args[4]
}
# Output plot filenames
rescaf_plot <- "rescaffolding.png"
# Prepend output directory to plot name (if it exists)
if (dir.exists(outdir)) {
  rescaf_plot <- paste0(outdir, rescaf_plot)
}

# Global variables
spc_int <- gsub("_"," ", spc_int)
label_names <- c("input", "rescaffolded", "chimera", "remainder", "excluded")

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
# Read table of species names and associated files
readSpecies <- function(assembly_file) {
  assembly_file_cols <- c("Species", "Assembly", "Annotation", "Proteins")
  checkPath(assembly_file)
  species_tab <- read.table(assembly_file, sep = "\t", fill = NA, header = F)
  colnames(species_tab) <- assembly_file_cols[1:length(species_tab)]
  species_tab <- species_tab[,na.omit(colnames(species_tab))]
  return(species_tab)
}
# Import NCBI AGP file (v2.0)
readAgp <- function(ragout_dir) {
  agp_cols <- c("seqid_scaf", "start_scaf", "end_scaf",
                "component_number", "component_type",
                "seqid_comp", "start_comp", "end_comp", "orientation")
  # AGP file describing scaffolding
  checkPath(ragout_dir)
  agp_file <- list.files(path = ragout_dir, pattern = ".*agp", full.names = T)
  agp <- read.table(agp_file, col.names = agp_cols)
  return(agp)
}
# Filters out gaps from AGP dataframe
filtAgp <- function(agp) {
  evidence <- agp %>% filter((component_type == "N")) %>% pull("orientation")
  # Calculate N content and total length before filtering
  agp_filt <- agp %>%
    mutate(N_count = sum(as.numeric(grep("scaffold_", seqid_comp,
                                         value = T, invert = T)), na.rm = T),
           total_len = sum(abs(as.numeric(start_scaf)-as.numeric(end_scaf)+1)),
           perc_N = round(N_count/total_len, digits = 4)) %>%
    filter(component_type == "W") %>%
    select(!contains("component")) %>%
    # Convert numeric columns
    mutate_at(vars(grep("start|end", colnames(.), value = T)), as.numeric)
  # agp_filt <- agp_filt %>%
  #   mutate(evidence = evidence)
  return(agp_filt)
}
# Keep gaps from AGP
gapsAgp <- function(agp) {
  gap_df <- agp %>%
    filter(component_type == "N") %>%
    rowwise() %>%
    mutate(evidence = length(unlist(strsplit(orientation, ","))),
           position = mean(c(start_scaf, end_scaf)))
  return(gap_df)
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
genIdx <- function(ragout_dir, agp_list, species_tab) {
  agp <- agp_list[[ragout_dir]]
  # Pre-Ragout scaffold index (later filtered for plotting)
  pre_file <- paste0(pull(filter(species_tab, grepl(spc_int, Species)),
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
  scaf <- scaf %>% mutate(type = "rescaffolded")
  # Type: remainder
  unplc <- unplc %>% mutate(type = "remainder") %>%
    rowwise() %>%
    # Filter out chimeric scaffolds from remainder
    filter(!grepl("\\[", seqid))
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
    mutate(type = factor(type, levels = label_names),
           seqid = factor(seqid, levels = unique(as.character(seqid))),
           length = log10(length)) %>%
    arrange(desc(length))
  return(idx)
}
# Index AGP dataframe with scaffold/contig lengths
idxAgp <- function(ragout_dir, filt_agp_list, idx_list) {
  filt_agp <- filt_agp_list[[ragout_dir]]
  idx <- idx_list[[ragout_dir]] %>% filter(type == "rescaffolded")
  pre_idx <- idx_list[[ragout_dir]] %>%
    filter(type %in% c("input", "chimera")) %>%
    arrange(desc(length)) %>%
    mutate(seqid = factor(seqid, levels = unique(as.character(seqid))))
  idx_dict <- idx$length
  names(idx_dict) <- idx$seqid
  pre_idx_dict <- pre_idx$length
  names(pre_idx_dict) <- pre_idx$seqid
  idx_agp <- filt_agp %>%
    # Convert length out of log10 scale
    mutate(length_scaf = 10^(as.numeric(idx_dict[seqid_scaf])),
           length_comp = 10^(as.numeric(pre_idx_dict[seqid_comp])),
           # Order component seqid factors by length (for plotting)
           seqid_comp = factor(seqid_comp, levels = levels(pre_idx$seqid))) %>%
    arrange(desc(length_scaf), start_scaf)
  return(idx_agp)
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
# Takes factor vector and returns factor-named vector of colors (i.e., dict)
getColors <- function(vec, palFun, pal = NULL) {
  if (!is.factor(vec)) stop("Error: vector class is not factor.")
  if (identical(palFun, brewer.pal)) {
    if (missing(pal)) pal = "Paired"
    vec_colors <- palFun(length(levels(vec)), pal)
    names(vec_colors) <- levels(vec)
    # Color chimera red using "Paired" palette
    if (any(names(vec_colors) %in% c("chimera"))) {
      new_ord <- c("input", "rescaffolded", "remainder", "excluded", "chimera")
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
                             simplify = F)))
  y_max <- max(unlist(sapply(idx_list, pull, length, simplify = F)))
  common_lims <- lims(x = c(0, x_max), y = c(0, y_max))
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
    mutate(type = factor(type, levels = c(label_names, "total"))) %>%
    rbind(list(factor("total", levels = c(label_names, "total")), sum(.$n)))
  # Colors for table
  tbl_colors <- getColors(idx$type, brewer.pal, pal)
  tbl_colors[["total"]] <- "lightgrey"
  fnt_colors <- rep("black", length(tbl_colors))
  names(fnt_colors) <- names(tbl_colors)
  fnt_colors[names(fnt_colors) %in% c("rescaffolded", "excluded")] <- "white"
  tbl_colors <- tbl_colors[as.character(unique(idx_sum$type))]
  fnt_colors <- fnt_colors[as.character(unique(idx_sum$type))]
  annot <- annotate(geom = "table", x = x_max, y = y_max,
                    label = idx_sum, table.colnames = F,
                    table.theme =
                      ttheme_default(base_size = 7,
                                     padding = unit(c(2, 2), "mm"),
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
    labs(title = ttl, x = "", y = "length [log10(bp)]") +
    theme_minimal() +
    annot +
    theme(legend.position = "left")
  return(p)
}
# Combine pre- and post-Ragout index plots
combPlot <- function(ragout_dir, idx_list, ttls = NULL, pal = "Paired",
                     legend = F) {
  p1 <- idxPlot(ragout_dir, idx_list, ttls = ttls, pal = pal,
                exclude = c("rescaffolded"))
  p2 <- idxPlot(ragout_dir, idx_list, pal = pal,
                exclude = c("input"))
  p_list <- list(p1, p2)
  if (legend) {
    common_legend <- getCommon(idxPlot, "type", idx_list)
    p <- ggarrange(plotlist = p_list, nrow = length(p_list),
                   align = "hv", legend = common_legend)
  } else {
    p <- ggarrange(plotlist = p_list, nrow = length(p_list),
                   align = "hv", legend = "none")
  }
  return(p)
}
# Plot of ragout rescaffolding
ragoutPlot <- function(ragout_dir, gaps_list, idx_agp_list, ttls = NULL) {
  ttl <- ttls[[ragout_dir]]
  gap_df <- gaps_list[[ragout_dir]]
  idx_agp <- idx_agp_list[[ragout_dir]]
  # Common x-axis limits for aligning all plots
  x_max <- max(unlist(sapply(idx_agp_list, pull, length_scaf)))
  common_x_lims <- xlim(c(0, x_max))
  # Color each seqid with consistent color in all graphs
  all_seqids <- lapply(idx_agp_list, select, c(seqid_comp, length_comp)) %>%
    bind_rows() %>%
    unique() %>%
    arrange(desc(length_comp)) %>%
    mutate(seqid_comp = as.numeric(fixChrom(seqid_comp))) %>%
    pull(seqid_comp)
  all_seqids <- factor(all_seqids, levels = all_seqids)
  # Factor seqids by scaffold length and start position (for plotting)
  order_scaf <- unique(pull(idx_agp, seqid_scaf))
  order_comp <- fixChrom(levels(pull(idx_agp, seqid_comp)))
  # Arrange data by scaffold length, then start position
  idx_agp <- idx_agp %>%
    mutate(seqid_comp = fixChrom(seqid_comp),
           seqid_comp = factor(seqid_comp, levels = order_comp),
           seqid_scaf = factor(seqid_scaf, levels = order_scaf))
  myColors <- getColors(all_seqids, viridis, "turbo")
  ypos <- dim(idx_agp)[1]
  xpos <- max(idx_agp$length_scaf)*0.8
  p <- ggplot(data = idx_agp) +
    geom_segment(mapping = aes(x = start_scaf, xend = end_scaf,
                               y = seqid_scaf, yend = seqid_scaf),
                               # col = chr0),
                               # col = length_comp),
                               # col = seqid_comp),
                 linewidth = 2) +
    annotate(geom = "text", x = xpos,
             y = idx_agp$seqid_scaf[ypos],
             label = paste(round(unique(idx_agp$total_len/10^6),
                                 digits = 2), "Mb")) +
    annotate(geom = "text", x = xpos,
             y = idx_agp$seqid_scaf[ypos-2],
             label = paste(round(unique(idx_agp$perc_N*100),
                                 digits = 2), "% N's")) +
    geom_segment(data = gap_df,
                 mapping = aes(x = start_scaf, xend = end_scaf,
                               y = seqid_scaf, yend = seqid_scaf,
                               col = evidence)) +
    # geom_text(data = gap_df,
    #            mapping = aes(x = position, y = seqid_scaf,
    #                          label = evidence)) +
    labs(title = ttl, x = "Length (bp)", y = "rescaffolded") +
    theme_classic() +
    scale_color_gradient(low = "yellow", high = "blue") +
    # scale_color_manual(values = myColors) +
    theme(legend.position = "right",
          # axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    common_x_lims
  return(p)
}

# Import data
species_tab <- readSpecies(assembly_file)
ragout_dirs <- c(ragout_dir)
ragout_dirs <- c("ragout-out-all-solid-refine")
# ragout_dirs <- c("ragout-out-filt-chimera",
#                  "ragout-out-all-chimera",
#                  "ragout-out-all-chimera-refine",
#                  "ragout-out-all-solid",
#                  "ragout-out-all-solid-refine",
#                  "ragout-out-no-sac-solid",
#                  "ragout-out-no-sac-solid-refine",
#                  "ragout-out-no-sac-chimera-refine")
# ttls <- c("Size filtration applied to all (chimeric)",
#           "No size filtration (chimeric)",
#           "No size filtration (chimeric, refined)",
#           "No size filtration (solid)",
#           "No size filtration (solid, refined)",
#           "Size filtration except S. latissima (solid)",
#           "Size filtration except S. latissima (solid, refined)",
#           "Size filtration except S. latissima (chimeric, refined)")
# names(ttls) <- ragout_dirs
agp_list <- sapply(ragout_dirs, readAgp, simplify = F)
filt_agp_list <- sapply(agp_list, filtAgp, simplify = F)
gaps_list <- sapply(agp_list, gapsAgp, simplify = F)
idx_list <- sapply(ragout_dirs, genIdx, filt_agp_list, species_tab,
                   simplify = F)
# Index AGP
idx_agp_list <- sapply(ragout_dirs, idxAgp, filt_agp_list, idx_list,
                         simplify = F)
# Annotate synteny to S. japonica chr0
idx_agp_list <- sapply(idx_agp_list, annotChr0, simplify = F)
 
# Plots
# Barplots of rescaffolded length distributions
p_list <- sapply(ragout_dirs, combPlot, idx_list,
                 # ttls = ttls,
                 simplify = F)
common_legend <- getCommon(idxPlot, "type", idx_list)
all_bar <- ggarrange(plotlist = p_list,  align = "hv",
                     legend.grob = common_legend, legend = "left")

# Plots of rescaffolding mapping
dot_list <- sapply(ragout_dirs, ragoutPlot, gaps_list, idx_agp_list,
                   simplify = F)
common_legend <- getCommon(ragoutPlot, "evidence", gaps_list, idx_agp_list)
all_dot <- ggarrange(plotlist = dot_list, align = "hv",
                     legend.grob = common_legend, legend = "right")
# Combine all plots into figure
comp_rag <- ggarrange(all_bar, all_dot, align = "hv")
comp_rag <- annotate_figure(comp_rag,
                            top = text_grob("Rescaffolding",
                                            face = "bold", size = 14))
showtext_opts(dpi = 300)
ggsave(rescaf_plot, comp_rag, bg = "white",
       height = 7, width = 12, units = "in")

