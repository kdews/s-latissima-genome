# Clear environment
rm(list = ls())
# Required packages
suppressPackageStartupMessages(library(tidyverse, quietly = T, warn.conflicts = F))
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
  assembly_file <- "s-latissima-genome/species_table.txt"
  ortho_dir <- "ortho_prots/OrthoFinder/Results_Mar08/"
  # Species of interest
  spc_int <- "Saccharina_latissima"
  # Outgroup species
  out_spc <- "Ectocarpus_siliculosus"
  # Output directory
  outdir <- "s-latissima-genome/"
} else {
  line_args <- commandArgs(trailingOnly = T)
  assembly_file <- line_args[1]
  ortho_dir <- line_args[2]
  # Species of interest
  spc_int <- line_args[3]
  # Outgroup species
  out_spc <- line_args[4]
  # Output directory
  outdir <- line_args[5]
}
# # Output plot filenames
# filt_plot <- paste0(spc_int, "_filtering.png")
# vio_plot <- "scaffold_sizes_violin.png"
# # Append output directory to plot name (if it exists)
# if (dir.exists(outdir)) {
#   # filt_plot <- paste0(outdir, filt_plot)
#   # vio_plot <- paste0(outdir, vio_plot)
# }

# Functions
# readOrthoTsv()

# Analysis
spc_int <- gsub("_"," ", spc_int)
species_tab <- read.table(assembly_file, sep = "\t", fill = NA, header = F)
colnames(species_tab) <- c("Species", "Assembly", "Annotation", "Proteins")
species_tab <- species_tab[,na.omit(colnames(species_tab))]
ortho_prefix <- basename(tools::file_path_sans_ext(species_tab$Proteins))
ortho_int <- ortho_prefix[grep(spc_int, species_tab$Species)]

ptn <- paste0(ortho_int, "\\.tsv")
list.files(path = "orthofinder/Results_Mar06/",
           pattern = ptn,
           full.names = T)



# read.table()

