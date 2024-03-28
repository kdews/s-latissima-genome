# Clear environment
rm(list = ls())
# Required packages
suppressPackageStartupMessages(library(tidyverse, quietly = T, warn.conflicts = F))
library(ggrepel)
library(ggpmisc)
library(viridisLite)
library(ape, quietly = T, warn.conflicts = F)
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
  assembly_file <- "s-latissima-genome/new_filt_species_table.txt"
  # Ragout output directory
  ragout_dir <- "ragout-out"
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
# # Output plot filenames
# filt_plot <- paste0(spc_int, "_filtering.png")
# vio_plot <- "scaffold_sizes_violin.png"
# # Append output directory to plot name (if it exists)
# if (dir.exists(outdir)) {
#   # filt_plot <- paste0(outdir, filt_plot)
#   # vio_plot <- paste0(outdir, vio_plot)
# }

# Global variables
spc_int <- gsub("_"," ", spc_int)
label_names <- c("input", "rescaffolded", "remainder", "excluded")

# Functions
# Checks existence of file or directory and errors if FALSE
checkPath <- function(test_path) {
  if (file.exists(test_path)) {
    return(TRUE)
  } else {
    stop(paste0("Error: cannot follow path (", test_path, ")."))
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
genIdx <- function(ragout_dir, agp_list, species_tab) {
  agp <- agp_list[[ragout_dir]]
  scaf_file <- list.files(path = ragout_dir,
                          pattern = "_scaffolds\\.fasta\\.fai",
                          full.names = T)
  unplc_file <- list.files(path = ragout_dir,
                           pattern = "_unplaced\\.fasta\\.fai",
                           full.names = T)
  pre_file <- paste0(pull(filter(species_tab, grepl(spc_int, Species)),
                           Assembly), ".fai")
  scaf <- readFai(scaf_file) %>% mutate(type = "rescaffolded")
  unplc <- readFai(unplc_file) %>% mutate(type = "remainder")
  rag_seqids <- c(agp$seqid_comp, gsub("\\[.*\\]", "", unplc$seqid))
  pre_idx <- readFai(pre_file) %>% mutate(type = "input")
    # mutate(type = case_when(seqid %in% agp$seqid_comp ~ "input",
    #                         .default = type))
  excl <- pre_idx %>% mutate(type = "excluded") %>%
    filter(!seqid %in% rag_seqids)
  idx <- rbind(scaf, unplc, excl, pre_idx) %>%
    mutate(type = factor(type, levels = label_names),
           length = log10(length)) %>%
    arrange(desc(length))
  return(idx)
}
# Index AGP dataframe with scaffold/contig lengths
idxAgp <- function(ragout_dir, filt_agp_list, idx_list) {
  filt_agp <- filt_agp_list[[ragout_dir]]
  idx <- idx_list[[ragout_dir]] %>% filter(type != "input")
  pre_idx <- idx_list[[ragout_dir]] %>% filter(type == "input")
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
getColors <- function(vec, pal = "Paired") {
  if (is.factor(vec)) {
    vec_colors <- brewer.pal(length(levels(vec)), pal)
    names(vec_colors) <- levels(vec)
  } else {
    vec_colors <- viridis(n = length(unique(vec)), option = pal)
    names(vec_colors) <- unique(vec)
  }
  return(vec_colors)
}
# Get common legend
getCommon <- function(idx_list) {
  type_count <- sapply(idx_list, function(x) length(unique(pull(x, type))),
                       simplify = F)
  legend_name <- names(which.max(type_count))
  common_legend <- get_legend(idxPlot(legend_name, idx_list))
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
  idx_sum <- idx %>%
    group_by(type) %>%
    summarise(n = n())
  if(!missing(exclude)) idx <- idx %>% filter(!type %in% exclude)
  if (any(idx$type %in% c("input"))) {
    filt_seqids <- idx %>% filter(type != "input") %>% pull(seqid)
    idx <- idx %>%
      filter(!(seqid %in% filt_seqids & type == "input"))
    annot <- NULL
  } else {
    tbl_colors <- getColors(idx$type, pal)
    fnt_colors <- rep("black", length(tbl_colors))
    names(fnt_colors) <- names(tbl_colors)
    fnt_colors[names(fnt_colors) %in% c("rescaffolded", "excluded")] <- "white"
    # tbl_colors <- tbl_colors[unique(idx$type)]
    # fnt_colors <- fnt_colors[unique(idx$type)]
    annot <- annotate(geom = "table", x = x_max, y = y_max,
                      label = list(idx_sum), table.colnames = F,
                      table.theme = ttheme_default(
                        core = list(bg_params = list(fill = tbl_colors),
                                    fg_params = list(col = fnt_colors))))
  }
  idx <- idx %>% arrange(type)
  vec_colors <- getColors(idx$type, pal)
  p <- ggplot(data = idx, mapping = aes(fill = type)) +
    geom_col(mapping = aes(x = as.numeric(row.names(idx)), y = length),
             width = 1) +
    common_lims +
    scale_fill_manual(values = vec_colors) +
    labs(title = ttl, x = "", y = "length [log10(bp)]") +
    theme_minimal() +
    # coord_cartesian(clip = "off") +
    annot +
    theme(legend.position = "top")
  return(p)
}
# Combine pre- and post-Ragout index plots
combPlot <- function(ragout_dir, idx_list, ttls = NULL, pal = "Paired",
                     legend = F) {
  p1 <- idxPlot(ragout_dir, idx_list, ttls = ttls, pal = pal,
                # exclude = c("rescaffolded", "remainder", "excluded"))
                exclude = c("rescaffolded"))
  p2 <- idxPlot(ragout_dir, idx_list, pal = pal,
                exclude = c("input"))
  if (legend) {
    common_legend <- getCommon(idx_list)
    p <- ggarrange(p1, p2, nrow = 2, align = "hv", legend.grob = common_legend)
  } else {
    p <- ggarrange(p1, p2, nrow = 2, align = "hv", legend = "none")
  }
  
  return(p)
}
# Plot of ragout rescaffolding
ragoutPlot <- function(ragout_dir, idx_agp_list, ttls, common_x_lims = NULL,
                      leg_help) {
  ttl <- ttls[[ragout_dir]]
  idx_agp <- idx_agp_list[[ragout_dir]]
  # Factor seqids by scaffold length and start position (for plotting)
  order_scaf <- unique(pull(idx_agp, seqid_scaf))
  order_comp <- fixChrom(levels(pull(idx_agp, seqid_comp)))
  # Arrange data by scaffold length, then start position
  idx_agp <- idx_agp %>%
    mutate(seqid_comp = fixChrom(seqid_comp),
           seqid_comp = factor(seqid_comp, levels = order_comp),
           seqid_scaf = factor(seqid_scaf, levels = order_scaf),
           temp = start_scaf,
           start_scaf = case_when(orientation == "-" ~ end_scaf,
                                  .default = start_scaf),
           end_scaf = case_when(orientation == "-" ~ temp,
                                .default = end_scaf))
  # myColors <- viridis(n = length(levels(idx_agp$seqid_comp)))
  # names(myColors) <- levels(idx_agp$seqid_comp)
  myColors <- viridis(n = length(levels(leg_help)), option = "turbo")
  names(myColors) <- levels(leg_help)
  ypos <- dim(idx_agp)[1]
  xpos <- max(idx_agp$length_scaf)*0.8
  p <- ggplot(data = idx_agp) +
    geom_segment(mapping = aes(x = min(start_scaf), xend = length_scaf,
                               y = seqid_scaf, yend = seqid_scaf),
                 col = "grey") +
    geom_segment(mapping = aes(x = start_scaf, xend = end_scaf,
                               y = seqid_scaf, yend = seqid_scaf,
                               # col = chr0),
                               # col = seqid_comp),
                               col = length_comp),
                 linewidth = 2) +
    annotate(geom = "text", x = xpos,
             y = idx_agp$seqid_scaf[ypos],
             label = paste(round(unique(idx_agp$total_len/10^6),
                                 digits = 2), "Mb")) +
    annotate(geom = "text", x = xpos,
             y = idx_agp$seqid_scaf[ypos-2],
             label = paste(round(unique(idx_agp$perc_N*100),
                                 digits = 2), "% N's")) +
    labs(title = ttl, x = "assembly v1.0", y = "rescaffolded") +
    theme_classic() +
    # scale_color_manual(values = myColors) +
    scale_color_viridis_c(option = "plasma") +
    theme(legend.position = "none", axis.title.y = element_blank()) +
    common_x_lims
  return(p)
}

# Import data
species_tab <- readSpecies(assembly_file)
ragout_dirs <- c("ragout-out-filt",
                 "ragout-out-chimera",
                 "ragout-out-solid",
                 "ragout-out",
                 "ragout-out-refine",
                 "ragout-out-refine-chimera")
ttls <- c("Size filtration applied to all (chimeric)",
          "No size filtration (chimeric)",
          "No size filtration (solid)",
          "Size filtration except S. latissima (solid)",
          "Size filtration except S. latissima (solid, refined)",
          "Size filtration except S. latissima (chimeric, refined)")
names(ttls) <- ragout_dirs
agp_list <- sapply(ragout_dirs, readAgp, simplify = F)
filt_agp_list <- sapply(agp_list, filtAgp, simplify = F)
idx_list <- sapply(ragout_dirs, genIdx, filt_agp_list, species_tab, simplify = F)
 
# test <- function(ragout_dir, filt_agp_list, idx_list) {
#   rescaf <- unique(filt_agp_list[[ragout_dir]]$seqid_comp)
#   remain <- unique(gsub("\\[.*\\]", "",
#                         as.character(idx_list[[ragout_dir]] %>%
#                                        filter(type == "remainder") %>%
#                                        pull(seqid))))
#   remain <- remain[!remain %in% rescaf]
#   print(sum(c(length(rescaf), length(remain))))
# }
# 
# sapply(ragout_dirs, test, filt_agp_list, idx_list)

# comp_file <- paste0(pull(filter(species_tab, grepl(spc_int, Species)),
#                          Assembly), ".fai")
# pre_idx <- readFai(comp_file) %>%
#   mutate(type = factor("input", levels = label_names), length = log10(length))
 
# # Index AGP
# idx_agp_list <- sapply(ragout_dirs, idxAgp, filt_agp_list, idx_list,
#                          simplify = F)
# # Annotate synteny to S. japonica chr0
# idx_agp_list <- sapply(idx_agp_list, annotChr0, simplify = F)
# 
# # Plots
# # Barplots of rescaffolded length distributions
# # x_max <- max(unlist(sapply(idx_list, dim, simplify = F)))
# # y_max <- max(unlist(sapply(idx_list, pull, length, simplify = F)))
# # common_lims <- lims(x = c(0, x_max), y = c(0, y_max))
# p_list <- sapply(ragout_dirs, idxPlot, idx_list, ttls = ttls,
#                  exclude = c("input"), simplify = F)
p_list <- sapply(ragout_dirs, combPlot, idx_list, ttls = ttls, simplify = F)
# pre_list <- sapply(ragout_dirs, idxPlot, idx_list, ttls = ttls,
#                    exclude = c("rescaffolded", "remainder", "excluded"), simplify = F)

# p_list2 <- list(`pre-ragout`= idxPlot(1, list(pre_idx),
#                                       list("Before rescaffolding"),
#                                       common_lims = ggplot_build(p_list[[1]])$layout$panel_scales_y[[1]]$range$range))
# p_list2 <- append(p_list2, p_list)
# # common_legend <- get_legend(p_list[[legend_name]])
common_legend <- getCommon(idx_list)
(comp_rag <- ggarrange(plotlist = p_list, align = "hv",
                       # ncol = length(p_list),
                       legend.grob = common_legend,
                       legend = "top"))
# # Plots of rescaffolding mapping
# common_x_lims <- xlim(c(0, max(unlist(sapply(idx_agp_list, pull, length_scaf)))))
# leg_help <- fixChrom(as.character(sort(unique(unlist(sapply(idx_agp_list, pull, seqid_comp))))))
# leg_help <- factor(leg_help, levels = leg_help)
# dot_list <- sapply(ragout_dirs, ragoutPlot, idx_agp_list, ttls, common_x_lims,
#                    leg_help, simplify = F)
# (all_dot <- ggarrange(plotlist = dot_list,
#                      nrow = 3, ncol = 2,
#                      align = "hv"))

