# Clear environment
rm(list = ls())
# Required packages
suppressPackageStartupMessages(library(tidyverse, quietly = T, warn.conflicts = F))
library(ape)
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
  # Directory containing input protein FASTAs for OrthoFinder
  ortho_dir <- "ortho_prots/"
  # Log file from OrthoFinder run
  log_file <- "orthofinder.log"
  # Species of interest
  spc_int <- "Saccharina_latissima"
  # Outgroup species
  out_spc <- "Ectocarpus_siliculosus"
  # Output directory
  outdir <- "s-latissima-genome/"
} else {
  line_args <- commandArgs(trailingOnly = T)
  assembly_file <- line_args[1]
  # Directory containing input protein FASTAs for OrthoFinder
  ortho_dir <- line_args[2]
  # Log file from OrthoFinder run
  log_file <- line_args[3]
  # Species of interest
  spc_int <- line_args[4]
  # Outgroup species
  out_spc <- line_args[5]
  # Output directory
  outdir <- line_args[6]
}
# # Output plot filenames
# filt_plot <- paste0(spc_int, "_filtering.png")
# vio_plot <- "scaffold_sizes_violin.png"
# # Append output directory to plot name (if it exists)
# if (dir.exists(outdir)) {
#   # filt_plot <- paste0(outdir, filt_plot)
#   # vio_plot <- paste0(outdir, vio_plot)
# }
spc_int <- gsub("_"," ", spc_int)
out_spc <- gsub("_"," ", out_spc)

# Global variables
c_names <- c("Species", "Assembly", "Annotation", "Proteins")

# Functions
# Checks existence of file or directory and errors if FALSE
checkPath <- function(test_path) {
  if (file.exists(test_path)) {
    return(TRUE)
  } else {
    stop(paste0("Error: cannot follow path (", test_path, ")."))
  }
}
# Parse JGI protein IDs
parseJGI <- function(x) {
  if (grepl("\\|", x)) {
    protein_ID <- unlist(strsplit(x, "\\|"))[3]
  } else {
    protein_ID <- x
  }
  return(protein_ID)
}
# Abbreviate species names
abbrevSpc <- function(spc) {
  spc <- unlist(strsplit(spc, " "))[1:2]
  let1 <- substr(spc[1], 1, 1)
  let2 <- substr(spc[2], 1, 1)
  spc_a <- toupper(paste0(let1, let2))
  return(spc_a)
}
# Read table of species names and associated files
readSpecies <- function(assembly_file) {
  checkPath(assembly_file)
  species_tab <- read.table(assembly_file, sep = "\t", fill = NA, header = F)
  colnames(species_tab) <- c_names
  species_tab <- species_tab[,na.omit(colnames(species_tab))]
  return(species_tab)
}
# Find OrthoFinder output directory path from given log file
getOrthoDir <- function(log_file) {
  checkPath(log_file)
  ortho_log <- readLines("orthofinder.log")
  res_dir <- gsub(" ", "", ortho_log[grep("Results:", ortho_log) + 1])
  return(res_dir)
}
# Correlate basename of OrthoFinder result files with species names
dictOrtho <- function(species_tab) {
  ortho_prefix <- basename(tools::file_path_sans_ext(species_tab$Proteins))
  ortho_dict <- species_tab$Species
  names(ortho_dict) <- ortho_prefix
  return(ortho_dict)
}
# Import TSV of orthologs for species as dataframe
readOrthoTsv <- function(spc, ortho_dict, res_dir) {
  spc <- grep(spc, ortho_dict, value = T)
  ortho_int <- names(ortho_dict[ortho_dict == spc])
  ptn <- paste0(ortho_int, "\\.tsv")
  logues_path <- paste0(res_dir, "Orthologues")
  checkPath(logues_path)
  ortho_tsv <- list.files(path = logues_path, pattern = ptn, full.names = T)
  df <- read.table(ortho_tsv, sep = "\t", header = T, comment.char = "")
  # Filter for single-copy orthologs
  df <- df %>%
    rowwise() %>%
    filter(all(!grepl(",", pick(colnames(df))))) %>%
    mutate_at(vars(ortho_int), parseJGI) %>%
    mutate(Species = ortho_dict[[Species]],
           Orthologs = parseJGI(Orthologs))
  tempcol <- grep(ortho_int, colnames(df))
  colnames(df)[tempcol] <- spc
  return(df)
}
# Read GFF3 file into dataframe
readGff3 <- function(spc, species_tab) {
  gff_file <- species_tab %>%
    filter(Species == spc) %>%
    pull(Annotation)
  gff <- read.gff(gff_file)
  gff <- gff %>%
    filter(type == "gene") %>%
    mutate(protein_ID = gsub(".*proteinId=|;.*", "", attributes)) %>%
    select(protein_ID, seqid, start, end)
  return(gff)
}
# Annotate orthologs with gene location information
idxOrtho <- function(spc, df, gff_list) {
  # Suffixes for merged dataframe
  suffs <- sapply(c(spc, spc_int), abbrevSpc, USE.NAMES = F)
  suffs <- paste0("_", suffs)
  # Given species
  gff1 <- gff_list[[spc]]
  # Species of interest
  if (any(grepl(spc_int, colnames(df)))) {
    spc_int1 <- grep(spc_int, colnames(df), value = T)
  } else {
    stop(paste0("Error: species of interest (", spc_int,
                ") not found in dataframe."))
  }
  gff2 <- gff_list[[spc_int1]]
  df1 <- df %>%
    filter(Species == spc) %>%
    select(Orthogroup, matches(spc_int1), Orthologs)
  df2 <- merge(df1, gff1, by.x = "Orthologs", by.y = "protein_ID")
  colnames(df2)[colnames(df2) == "Orthologs"] <- spc
  df3 <- merge(df2, gff2, by.x = spc_int1, by.y = "protein_ID",
               suffixes = suffs)
  cols1 <- c("Orthogroup", colnames(df3)[colnames(df3) != "Orthogroup"])
  df3 <- df3 %>%
    arrange(Orthogroup) %>%
    select(matches(cols1))
  return(df3)
}

# Import data
# Species table
species_tab <- readSpecies(assembly_file)
# OrthoFinder results directory
res_dir <- getOrthoDir(log_file)
# Dictionary of OrthoFinder file prefixes and species names
ortho_dict <- dictOrtho(species_tab)
# Results for species of interest
df <- readOrthoTsv(spc_int, ortho_dict, res_dir)
# Import GFF3s as dataframes
gff_list <- sapply(unname(ortho_dict), readGff3, species_tab, simplify = F)

df2 <- idxOrtho(unname(ortho_dict)[1], df, gff_list)
head(df2)


# Analysis
# df %>%
#   group_by(Species) %>%
#   summarise(n = n())
# ptn <- paste0("OrthologuesStats_one-to-one", "\\.tsv")
# stats_path <- paste0(res_dir, "Comparative_Genomics_Statistics")
# sc_stats_file <- list.files(path = stats_path, pattern = ptn, full.names = T)
# test <- read.table(sc_stats_file, sep = "\t", header = T, row.names = 1)
# heatmap(as.matrix(test))





