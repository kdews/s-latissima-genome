# Clear environment
rm(list = ls())
# Required packages
suppressPackageStartupMessages(library(tidyverse, quietly = T, warn.conflicts = F))
library(ape)
library(scales)
library(gridExtra, quietly = T, warn.conflicts = F)
library(ggpubr, quietly = T)
if (require(showtext, quietly = T)) {
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
  colnames(species_tab) <- assembly_file_cols
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
  agp_filt <- agp %>%
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
    mutate(seqid = factor(seqid, levels = seqid),
           rank=as.numeric(rownames(.)))
  if (!missing(suff)) {
    colnames(idx) <- paste0(colnames(idx), suff)
  }
  return(idx)
}
# Dotplot of ragout rescaffolding
ragoutDot <- function(gff_update) {
  # Arrange dataframe by scaffold length, then start position
  gff_update <- gff_update %>%
    mutate(seqid_comp = fixChrom(seqid_comp)) %>%
    arrange(desc(length_scaf), start_scaf)
  # Factor seqids by scaffold length and start position (for plotting)
  order_scaf <- unique(pull(gff_update, seqid_scaf))
  order_comp <- unique(pull(gff_update, seqid_comp))
  gff_update <- gff_update %>%
    mutate(seqid_scaf = factor(seqid_scaf, levels = order_scaf),
           seqid_comp = factor(seqid_comp, levels = order_comp),
           temp = start_scaf,
           start_scaf = case_when(orientation == "-" ~ end_scaf,
                                  .default = start_scaf),
           end_scaf = case_when(orientation == "-" ~ temp,
                                .default = end_scaf))
  p <- ggplot(data = gff_update) +
    geom_segment(mapping = aes(x = start_scaf, xend = end_scaf,
                               y = start_comp, yend = end_comp),
                 col = "black") +
    facet_grid(cols = vars(seqid_scaf), rows = vars(seqid_comp),
               switch = "both", scales = "free", as.table = F) +
    labs(title = "Rescaffolding with Ragout", x = "Ragout assembly",
         y = "Assembly v1.0") +
    theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank(),
          strip.placement = "outside",
          strip.background.x = element_rect(color = NA,  fill=NA),
          strip.background.y = element_rect(color = NA,  fill=NA),
          strip.text.x.bottom = element_text(size = rel(0.75), angle = 0),
          strip.text.y.left = element_text(size = rel(0.75), angle = 0),
          strip.clip = "off",
          legend.position = "none")
  return(p)
}
idxPlot <- function(idx, y_axis_max, pal) {
  # idx <- idx %>%
  #   filter(type != "excluded")
  p <- ggplot(data = idx, mapping = aes(fill = type)) +
    geom_col(mapping = aes(x = rank, y = length), width = 1) +
    ylim(0, y_axis_max) +
    scale_fill_manual(values = pal) +
    # scale_fill_brewer(palette = pal, direction = -1) +
    # scale_fill_viridis_d(option = pal, direction = -1) +
    theme_bw()
  return(p)
}


# Import data
species_tab <- readSpecies(assembly_file)
agp <- readAgp(ragout_dir)
idx_file_comp <- paste0(pull(filter(species_tab, grepl(spc_int, Species)),
                             Assembly), ".fai")
idx_file_scaf <- list.files(path = ragout_dir,
                            pattern = "_scaffolds\\.fasta\\.fai",
                            full.names = T)
idx_file_unscaf <- list.files(path = ragout_dir,
                            pattern = "_unplaced\\.fasta\\.fai",
                            full.names = T)
idx_comp <- readFai(idx_file_comp) %>%
  mutate(type = case_when(rank <= 155 ~ "Ragout input", .default = "excluded"),
         type = factor(type, levels = c("Ragout input", "Ragout assembly",
                                        "remainder", "excluded")))
idx_scaf <- readFai(idx_file_scaf) %>%
  mutate(type = "Ragout assembly")
idx_unscaf <- readFai(idx_file_unscaf) %>%
  mutate(type = "remainder",
         rank = rank + max(idx_scaf$rank))
idx_ragout <- rbind(idx_scaf, idx_unscaf,
                    idx_comp %>% filter(type == "excluded")) %>% 
  mutate(rank = as.numeric(row.names(.)),
         type = factor(type, levels = c("Ragout input", "Ragout assembly",
                                        "remainder", "excluded")))
y_axis_max <- max(c(idx_comp$length, idx_ragout$length))

vec.colors <- brewer.pal(4, "Paired")
names(vec.colors) <- levels(idx_comp$type)

# Plots
p1 <- idxPlot(idx_comp, y_axis_max, vec.colors)
p2 <- idxPlot(idx_ragout, y_axis_max, vec.colors)
p3 <- idxPlot(idx_comp %>% filter(type != "excluded"), y_axis_max, vec.colors)
p4 <- idxPlot(idx_ragout %>% filter(type != "excluded"), y_axis_max, vec.colors)
ggarrange(p1, p3, p2, p4, ncol = 2, nrow = 2, align = "hv")

