# Clear environment
rm(list = ls())
# Required packages
suppressPackageStartupMessages(library(tidyverse, quietly = T, warn.conflicts = F))
library(ggrepel)
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
  agp <- read.table(agp_file)
  colnames(agp) <- agp_cols
  return(agp)
}
# Filters out gaps from AGP dataframe
filtAgp <- function(agp) {
  agp_filt <- agp %>%
    mutate(N_count = sum(as.numeric(seqid_comp), na.rm = T),
           total_len = sum(abs(start_scaf-end_scaf+1)),
           perc_N = round(N_count/total_len, digits = 4)) %>%
    filter(component_type == "W") %>%
    select(!contains("component"))
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
  comp_file <- paste0(pull(filter(species_tab, grepl(spc_int, Species)),
                           Assembly), ".fai")
  scaf <- readFai(scaf_file) %>% mutate(type = "rescaffolded")
  unplc <- readFai(unplc_file) %>% mutate(type = "remainder")
  comp <- readFai(comp_file) %>%
    mutate(type = case_when(seqid %in% agp$seqid_comp ~ "input",
                            as.numeric(row.names(.)) <= 155 ~ "remainder",
                            .default = "excluded")) %>%
    filter(!seqid %in% unique(gsub("\\[.*\\]", "", unplc$seqid)))
  idx <- rbind(scaf, unplc, comp) %>%
    mutate(type = factor(type, levels = label_names),
           # length = length) %>%
           length = log10(length)) %>%
    arrange(desc(length))
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
# Barplot of scaffolds indexed by length
idxPlot <- function(ragout_dir, idx_list, ttls = NULL, pal = "Paired",
                    common_lims = NULL, include = NULL, exclude = NULL) {
  idx <- idx_list[[ragout_dir]]
  ttl <- ttls[[ragout_dir]]
  if(!missing(include)) idx <- idx %>% filter(type %in% include)
  if(!missing(exclude)) idx <- idx %>% filter(!type %in% exclude)
  idx <- idx %>% arrange(type)
  vec_colors <- getColors(idx$type, pal)
  p <- ggplot(data = idx, mapping = aes(fill = type)) +
    geom_col(mapping = aes(x = as.numeric(row.names(idx)), y = length),
             width = 1) +
    common_lims +
    scale_fill_manual(values = vec_colors) +
    labs(title = ttl, x = "", y = "length [log10(bp)]") +
    theme_minimal()
  return(p)
}
# Combine pre- and post-Ragout index plots
combPlot <- function(idx, pal = "Paired", common_lims = NULL) {
  common_legend <- get_legend(idxPlot(idx, pal))
  p1 <- idxPlot(idx, pal, exclude = c("rescaffolded")) + common_lims
  p2 <- idxPlot(idx, pal, exclude = c("input")) + common_lims
  p <- ggarrange(p1, p2, nrow = 2,
                 # legend.grob = common_legend,
                 legend = "none",
                 align = "hv")
  return(p)
}
# Merge agp and lengths
mergeAgp <- function(ragout_dir, filt_agp_list, idx_list, pre_idx) {
  filt_agp <- filt_agp_list[[ragout_dir]]
  idx <- idx_list[[ragout_dir]]
  # !!!!! merge with pre_idx
  merge_agp <- merge(filt_agp, idx, by.x = "seqid_scaf", by.y = "seqid") %>%
    mutate(length_scaf = 10^length,
           seqid_comp = factor(seqid_comp, levels = levels(pre_idx$seqid))) %>%
    mutate_at(vars(grep("start|end", colnames(.), value = T)), as.numeric) %>%
    arrange(desc(length_scaf), start_scaf)
  return(merge_agp)
}
# Annotate synteny to chr0
annotChr0 <- function(merge_agp) {
  comp_chr0 <- read.table("comp_chr0.txt", header = T, sep = "\t")
  comp_chr0 <- comp_chr0 %>%
    mutate(qName = paste0("scaffold_", qNum))
  dict_chr0 <- pull(comp_chr0, qPercent)
  names(dict_chr0) <- pull(comp_chr0, qName)
  merge_agp <- merge_agp %>%
    mutate(chr0 = case_when(seqid_comp %in% names(dict_chr0) ~
                              dict_chr0[as.character(seqid_comp)],
                            .default = 0))
  return(merge_agp)
}
# Dotplot of ragout rescaffolding
ragoutDot <- function(ragout_dir, merge_agp_list, ttls, common_x_lims = NULL,
                      leg_help) {
  ttl <- ttls[[ragout_dir]]
  merge_agp <- merge_agp_list[[ragout_dir]]
  # Factor seqids by scaffold length and start position (for plotting)
  order_scaf <- unique(pull(merge_agp, seqid_scaf))
  order_comp <- fixChrom(levels(pull(merge_agp, seqid_comp)))
  # Arrange dataframe by scaffold length, then start position
  merge_agp <- merge_agp %>%
    mutate(seqid_comp = fixChrom(seqid_comp),
           seqid_comp = factor(seqid_comp, levels = order_comp),
           seqid_scaf = factor(seqid_scaf, levels = order_scaf),
           temp = start_scaf,
           start_scaf = case_when(orientation == "-" ~ end_scaf,
                                  .default = start_scaf),
           end_scaf = case_when(orientation == "-" ~ temp,
                                .default = end_scaf))
  # myColors <- viridis(n = length(levels(merge_agp$seqid_comp)))
  # names(myColors) <- levels(merge_agp$seqid_comp)
  myColors <- viridis(n = length(levels(leg_help)), option = "turbo")
  names(myColors) <- levels(leg_help)
  ypos <- dim(merge_agp)[1]
  xpos <- max(merge_agp$length_scaf)*0.8
  p <- ggplot(data = merge_agp) +
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
             y = merge_agp$seqid_scaf[ypos],
             label = paste(round(unique(merge_agp$total_len/10^6),
                                 digits = 2), "Mb")) +
    annotate(geom = "text", x = xpos,
             y = merge_agp$seqid_scaf[ypos-2],
             label = paste(round(unique(merge_agp$perc_N*100),
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
comp_file <- paste0(pull(filter(species_tab, grepl(spc_int, Species)),
                         Assembly), ".fai")
pre_idx <- readFai(comp_file) %>%
  mutate(type = factor("input", levels = label_names),
         length = log10(length))
legend_name <- names(which.max(sapply(idx_list,
                                      function(x) length(unique(pull(x, type))),
                                      simplify = F)))
common_legend <- get_legend(idxPlot(legend_name, idx_list))
x_max <- max(unlist(sapply(idx_list, dim, simplify = F)))
y_max <- max(unlist(sapply(idx_list, pull, length, simplify = F)))
common_lims <- lims(x = c(0, x_max), y = c(0, y_max))
p_list <- sapply(ragout_dirs, idxPlot, idx_list, ttls = ttls, common_lims = common_lims, 
               exclude = c("input"), simplify = F)
p_list2 <- list(`pre-ragout`= idxPlot(1, list(pre_idx),
                                     list("Before rescaffolding"),
                                     common_lims = common_lims))
p_list2 <- append(p_list2, p_list)
# (comp_rag <- ggarrange(plotlist = p_list2,
#                       nrow = length(p_list2),
#                       legend.grob = common_legend, legend = "right",
#                       align = "hv"))
merge_agp_list <- sapply(ragout_dirs, mergeAgp, filt_agp_list, idx_list, pre_idx, simplify = F)
merge_agp_list <- sapply(merge_agp_list, annotChr0, simplify = F)
common_x_lims <- xlim(c(0, max(unlist(sapply(merge_agp_list, pull, length_scaf)))))
leg_help <- fixChrom(as.character(sort(unique(unlist(sapply(merge_agp_list, pull, seqid_comp))))))
leg_help <- factor(leg_help, levels = leg_help)
dot_list <- sapply(ragout_dirs, ragoutDot, merge_agp_list, ttls, common_x_lims,
                   leg_help, simplify = F)
(all_dot <- ggarrange(plotlist = dot_list,
                     nrow = 3, ncol = 2,
                     align = "hv"))

