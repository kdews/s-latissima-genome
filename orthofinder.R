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
  # Path to results from OrthoFinder run
  res_dir <- "ortho_prots/OrthoFinder/Results_Mar08/"
  # Ragout output directory
  ragout_dir <- "ragout-out"
  # Species of interest
  spc_int <- "Saccharina_latissima"
  # Output directory
  outdir <- "s-latissima-genome/"
} else {
  line_args <- commandArgs(trailingOnly = T)
  assembly_file <- line_args[1]
  # Directory containing input protein FASTAs for OrthoFinder
  ortho_dir <- line_args[2]
  # Path to results from OrthoFinder run
  res_dir <- line_args[3]
  # Ragout output directory
  ragout_dir <- line_args[4]
  # Species of interest
  spc_int <- line_args[5]
  # Output directory
  outdir <- line_args[6]
}
# # Output plot filenames
# filt_plot <- paste0(spc_int, "_filtering.png")
# vio_plot <- "scaffold_sizes_violin.png"
# # Prepend output directory to plot name (if it exists)
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
    stop(paste0("Cannot follow path (", test_path, ")."))
  }
}
# Parse JGI protein IDs
parseJGI <- function(x) {
  if (grepl("\\|", x)) {
    protein_ID <- unlist(strsplit(x, "\\|"))[3]
  } else {
    protein_ID <- gsub("\\..*", "", x)
  }
  return(protein_ID)
}
# Fixes chromosome labels for plotting
fixChrom <- function(contigs) {
  contigs <-
    as.character(as.numeric(str_remove_all(str_remove_all(contigs, ".*_"),
                                           "[^0-9]")))
  return(contigs)
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
  assembly_file_cols <- c("Species", "Assembly", "Annotation", "Proteins")
  checkPath(assembly_file)
  species_tab <- read.table(assembly_file, sep = "\t", fill = NA, header = F,
                            col.names = assembly_file_cols)
  species_tab <- species_tab[,na.omit(colnames(species_tab))]
  return(species_tab)
}
# # Find OrthoFinder output directory path from given log file
# getOrthoDir <- function(log_file) {
#   checkPath(log_file)
#   ortho_log <- readLines("orthofinder.log")
#   res_dir <- gsub(" ", "", ortho_log[grep("Results:", ortho_log) + 1])
#   checkPath(res_dir)
#   return(res_dir)
# }
# Correlate basename of OrthoFinder result files with species names
dictOrtho <- function(species_tab) {
  ortho_prefix <- basename(tools::file_path_sans_ext(species_tab$Proteins))
  ortho_dict <- species_tab$Species
  names(ortho_dict) <- ortho_prefix
  return(ortho_dict)
}
# Import TSV of orthologs for species as dataframe
readOrthoTsv <- function(spc, ortho_dict, res_dir) {
  spc <- unname(grep(spc, ortho_dict, value = T))
  ortho_int <- as.character(names(ortho_dict[ortho_dict == spc]))
  ptn <- paste0(ortho_int, "\\.tsv")
  logues_path <- paste0(res_dir, "Orthologues")
  checkPath(logues_path)
  ortho_tsv <- list.files(path = logues_path, pattern = ptn, full.names = T)
  df <- read.table(ortho_tsv, sep = "\t", header = T, comment.char = "")
  df <- df %>%
    rowwise() %>%
    # # Filter for single-copy orthologs
    # filter(all(!grepl(",", pick(colnames(df))))) %>%
    # mutate_at(vars(ortho_int), parseJGI) %>%
    mutate(Species = as.character(ortho_dict[[Species]]))
  # Orthologs = parseJGI(Orthologs))
  tempcol <- grep(ortho_int, colnames(df))
  colnames(df)[tempcol] <- spc
  return(df)
}
# Read GFF3 file into dataframe
readGff <- function(spc, species_tab) {
  gff_file <- species_tab %>%
    filter(Species == spc) %>%
    pull(Annotation)
  gff <- read.gff(gff_file)
  gff <- gff %>%
    filter(type == "gene") %>%
    rowwise() %>%
    mutate(protein_ID = parseJGI(gsub(".*Name=|;.*", "", attributes)),
           length = abs(end-start)) %>%
    select(protein_ID, seqid, start, end, length)
  return(gff)
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
  if (!missing(suff)) {
    colnames(idx) <- paste0(colnames(idx), suff)
  }
  return(idx)
}
# Index AGP dataframe with component and scaffolded FASTA indices
idxAgp <- function(agp, idx_comp, idx_scaf) {
  agp_idx <- merge(agp, idx_comp, by = "seqid_comp", sort = F)
  agp_idx <- merge(agp_idx, idx_scaf, by = "seqid_scaf", sort = F)
  return(agp_idx)
}
# Recalibration of GFF3 for scaffolded genome
calGff <- function(agp_idx, gff_list) {
  spc <- grep(spc_int, names(gff_list))
  gff <- gff_list[[spc]]
  gff_update <- merge(agp_idx, gff,
                      by.x = "seqid_comp", by.y = "seqid", sort = F)
  gff_update <- gff_update %>%
    mutate_at(vars(contains("start"), contains("end"), contains("length")),
              as.numeric) %>%
    # Filter for genes within given range of component
    filter(start_comp <= start,
           start_comp < end,
           end_comp >= end,
           end_comp > start) %>%
    # Calculate gene coordinate shift caused by rescaffolding
    mutate(new_start =
             case_when(orientation == "+" ~ start - start_comp + start_scaf,
                       orientation == "-" ~ end_comp - start + start_scaf),
           new_end = new_start + length) %>%
    # Arrange dataframe by scaffold length, then start position
    arrange(desc(length_scaf), start_scaf)
  # Factor seqids by scaffold length and start position (for plotting)
  order_scaf <- unique(pull(gff_update, seqid_scaf))
  order_comp <- unique(pull(gff_update, seqid_comp))
  gff_update <- gff_update %>%
    mutate(seqid_scaf = factor(seqid_scaf, levels = order_scaf),
           seqid_comp = factor(seqid_comp, levels = order_comp))
  gff_list[[spc]] <- gff_update
  return(gff_list)
}
# Annotate orthologs with gene location information
idxOrtho <- function(spc, orthos, gff_list) {
  # Suffixes for merged dataframe
  suffs <- sapply(c(spc, spc_int), abbrevSpc, USE.NAMES = F)
  suffs <- paste0("_", suffs)
  # Given species
  gff1 <- gff_list[[spc]]
  # Species of interest
  if (any(grepl(spc_int, colnames(orthos)))) {
    spc_int1 <- grep(spc_int, colnames(orthos), value = T)
  } else {
    stop(paste0("Error: species of interest (", spc_int,
                ") not found in dataframe."))
  }
  gff2 <- gff_list[[spc_int1]]
  df1 <- orthos %>%
    filter(Species == spc) %>%
    # Filter for single-copy orthologs
    filter(!grepl(",", Orthologs)) %>%
    filter(all(!grepl(",", pick(spc_int1)))) %>%
    mutate_at(vars(spc_int1, Orthologs), parseJGI) %>%
    select(Orthogroup, matches(spc_int1), Orthologs)
  print(head(df1, n = 2))
  print(dim(df1))
  # Merge with GFF3 of given species
  df2 <- merge(df1, gff1, by.x = "Orthologs", by.y = "protein_ID", sort = F)
  print(head(df2, n = 2))
  print(dim(df2))
  colnames(df2)[colnames(df2) == "Orthologs"] <- spc
  # Merge with GFF3 of species of interest
  df3 <- merge(df2, gff2, by.x = spc_int1, by.y = "protein_ID",
               suffixes = suffs, sort = F)
  print(head(df3, n = 2))
  print(dim(df3))
  return(df3)
}

# Import data
# Species table
species_tab <- readSpecies(assembly_file)
# OrthoFinder results directory
# res_dir <- getOrthoDir(log_file)
res_dir <- gsub(".*ortho_prots", "ortho_prots", res_dir)
# Dictionary of OrthoFinder file prefixes and species names
ortho_dict <- dictOrtho(species_tab)
# Results for species of interest
orthos <- readOrthoTsv(spc_int, ortho_dict, res_dir)
# Import GFF3s as dataframes
gff_list <- sapply(unname(ortho_dict), readGff, species_tab, simplify = F)
# Dataframe describing scaffolding
agp <- readAgp(ragout_dir)
# FASTA indices of component and scaffolded FASTAs
idx_file_comp <- paste0(pull(filter(species_tab, grepl(spc_int, Species)),
                             Assembly), ".fai")
idx_file_scaf <- list.files(path = ragout_dir,
                            pattern = "_scaffolds\\.fasta\\.fai",
                            full.names = T)
idx_comp <- readFai(idx_file_comp, "_comp")
idx_scaf <- readFai(idx_file_scaf, "_scaf")

# Analysis
# Index AGP dataframe with component and scaffolded FASTA indices (lengths)
agp_idx <- idxAgp(agp, idx_comp, idx_scaf)
# Recalibrate gene positions with rescaffolded genome
gff_list_cal <- calGff(agp_idx, gff_list)
# Merge GFF3s with OrthoFinder results
# Original coordinates
o_idx_list <- sapply(grep(spc_int, names(gff_list), invert = T, value = T),
                     idxOrtho, orthos, gff_list, simplify = F)
# Calibrated rescaffolded coordinates
o_idx_list_cal <- sapply(grep(spc_int, names(gff_list_cal), invert = T, value = T),
                         idxOrtho, orthos, gff_list_cal, simplify = F)


rag_seqs <- unique(as.character(unlist(sapply(o_idx_list_cal, pull, "seqid_comp"))))
gene_density <- merge(gff_list$`Saccharina latissima SL-CT1-FG3 v1.0`, idx_comp,
                      by.x = "seqid", by.y = "seqid_comp") %>%
  select(seqid, length_comp, protein_ID) %>%
  group_by(seqid, length_comp) %>%
  summarize(n = n()) %>%
  mutate(ragout = seqid %in% rag_seqs)
ggplot(data = gene_density, mapping = aes(x = n)) +
  geom_histogram()
ggplot(data = gene_density, mapping = aes(x = length_comp, y = n)) +
  geom_point(aes(col = ragout))

# df <- o_idx_list_cal[["Undaria pinnatifida M23"]] %>%
#   select(!contains("length")) %>%
#   select(contains("seqid"), contains("_UP"), contains("new"))
# write.table(x = df, file = "for_jose_undaria.txt", quote = F, sep = "\t",
#             row.names = F)




orthos %>%
  filter(grepl("Und", Species)) %>%
  filter(all(!grepl(",", pick(colnames(orthos))))) %>%
  group_by(Species) %>%
  summarise(n = n())
ptn <- paste0("OrthologuesStats_one-to-one", "\\.tsv")
stats_path <- paste0(res_dir, "Comparative_Genomics_Statistics")
sc_stats_file <- list.files(path = stats_path, pattern = ptn, full.names = T)
test1 <- read.table(sc_stats_file, sep = "\t", header = T, row.names = 1)
heatmap(as.matrix(test1))
