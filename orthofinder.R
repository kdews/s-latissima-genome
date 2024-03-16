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
  # Ragout output directory
  ragout_dir <- "ragout-out/"
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
  # Ragout output directory
  ragout_dir <- line_args[5]
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


# Global variables
spc_int <- gsub("_"," ", spc_int)

assembly_file_cols <- c("Species", "Assembly", "Annotation", "Proteins")
agp_cols <- c("seqid_scaf", "start_scaf", "end_scaf",
              "component_number", "component_type",
              "seqid_comp", "start_comp", "end_comp", "orientation")

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
  checkPath(assembly_file)
  species_tab <- read.table(assembly_file, sep = "\t", fill = NA, header = F)
  colnames(species_tab) <- assembly_file_cols
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
    mutate(protein_ID = gsub(".*id=|;.*", "", attributes, ignore.case = T),
           length = abs(end-start)) %>%
    select(protein_ID, seqid, start, end, length)
  return(gff)
}
# Import NCBI AGP file (v2.0)
readAgp <- function(ragout_dir) {
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
  return(list(gff2, df2, df3))
}
# Plot
plotCal <- function(gff_update) {
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

# Import data
# Species table
species_tab <- readSpecies(assembly_file)
# OrthoFinder results directory
res_dir <- getOrthoDir(log_file)
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
idx_file_scaf <- list.files(path = ragout_dir, pattern = "scaffolds.*fai",
                            full.names = T)
idx_comp <- readFai(idx_file_comp, "_comp")
idx_scaf <- readFai(idx_file_scaf, "_scaf")

# Analysis
# Index AGP dataframe with component and scaffolded FASTA indices (lengths)
agp_idx <- idxAgp(agp, idx_comp, idx_scaf)
# Recalibrate gene positions with rescaffolded genome
gff_list_cal <- calGff(agp_idx, gff_list)
# Merge GFF3s with OrthoFinder results
o_idx_list <- sapply(grep(spc_int, names(gff_list), invert = T, value = T),
                     idxOrtho, orthos, gff_list, simplify = F)
o_idx_list_cal <- sapply(grep(spc_int, names(gff_list_cal), invert = T, value = T),
                     idxOrtho, orthos, gff_list_cal, simplify = F)


test_gff <- o_idx_list$`Undaria pinnatifida M23`[[1]]
test_df2 <- o_idx_list$`Undaria pinnatifida M23`[[2]]

length(which(test_df2$`Saccharina latissima SL-CT1-FG3 v1.0` %in% unname(ann_dict[test_gff$protein_ID])))

gff_list$`Ectocarpus siliculosus Ec 32 V2` %>%
  filter()

ann_idx_file <- "/scratch1/kdeweese/sterility/annotation_index.txt"
ann_idx <- read.table(ann_idx_file, header = T)
ann_dict <- c(as.character(ann_idx[,"proteinId"]))
names(ann_dict) <- as.character(ann_idx[,"transcriptId"])


dim(merge(test_df2, test_gff,
      by.x = "Saccharina latissima SL-CT1-FG3 v1.0", by.y = "protein_ID",
      suffixes = c("_UP", "_SL"),
      sort = F))

# df <- o_idx_list[["Undaria pinnatifida M23"]] %>%
#   select(!contains("length")) %>%
#   select(seqid, contains("_UP"), contains("_scaf"))
# write.table(x = df, file = "for_jose_undaria.txt", quote = F, sep = "\t",
#             row.names = F)

# Plots
((p <- plotCal(gff_list[["Saccharina latissima SL-CT1-FG3 v1.0"]])))
# test1 <- gff_list[["Saccharina latissima SL-CT1-FG3 v1.0"]] %>%
#   filter(seqid_scaf == "chr_chr1")
# test2 <- test1 %>%
#   mutate(temp = new_start, temp2 = start_scaf,
#          new_start = case_when(orientation == "-" ~ new_end,
#                                .default = new_start),
#          new_end = case_when(orientation == "-" ~ temp,
#                              .default = new_end),
#          start_scaf = case_when(orientation == "-" ~ end_scaf,
#                                 .default = start_scaf),
#          end_scaf = case_when(orientation == "-" ~ temp2,
#                               .default = end_scaf))
# ((p1 <- plotCal(test1)))
# ((p2 <- plotCal(test2)))


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


fname <- "/scratch1/kdeweese/latissima/genome_stats/ortho_prots/OrthoFinder/Results_Mar08/Orthologues/Orthologues_SlaSLCT1FG3_1_GeneCatalog_proteins_20210608_aa/SlaSLCT1FG3_1_GeneCatalog_proteins_20210608_aa__v__Undpi1_1_FilteredModels1_proteins_2023-12-12.tsv"
test <- read.table(fname, sep = "\t", header = T, comment.char = "")
dim(test)
test2 <- test %>%
  rowwise() %>%
  # filter(all(grepl(",", pick(colnames(test)))))
  filter(!grepl(",", SlaSLCT1FG3_1_GeneCatalog_proteins_20210608_aa))
  # filter(!grepl(",", Undpi1_1_FilteredModels1_proteins_2023.12.12))
dim(test2)


und_os <- orthos %>% filter(grepl("Und", Species)) %>% pull(Orthologs) %>% grep(",", ., value = T, invert = T)
und_os <- unname(unlist(sapply(und_os, parseJGI, simplify = F)))
und_gs <- unique(gff_list$`Undaria pinnatifida M23`$protein_ID)
length(which(und_os %in% und_gs))
