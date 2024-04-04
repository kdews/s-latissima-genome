# Clear environment
rm(list = ls())
# Required packages
suppressPackageStartupMessages(library(tidyverse, quietly = T, warn.conflicts = F))
library(ape, quietly = T, warn.conflicts = F)

# Input
if (interactive()) {
  wd <- "/scratch1/kdeweese/latissima/genome_stats"
  setwd(wd)
  # Cactus-formatted seqFile
  seqFile <- "s-latissima-genome/s_lat_alignment.txt"
  # Cactus output directory
  cactus_dir <- "cactus-steps-output"
  # Species of interest
  spc_int <- "Saccharina_latissima"
  # Output directory
  outdir <- "s-latissima-genome/"
} else {
  line_args <- commandArgs(trailingOnly = T)
  # Cactus-formatted seqFile
  seqFile <- line_args[1]
  # Cactus output directory
  cactus_dir <- line_args[2]
  # Species of interest
  spc_int <- line_args[3]
  # Output directory
  outdir <- line_args[4]
}
# Output filename
cactus_prefix <- basename(tools::file_path_sans_ext(seqFile))
rcp_file <- paste0(cactus_prefix, ".rcp")
# Prepend output directory to filename (if it exists)
if (dir.exists(outdir)) {
  rcp_file <- paste0(outdir, rcp_file)
}

# Checks existence of file or directory and errors if FALSE
checkPath <- function(test_path) {
  if (file.exists(test_path)) {
    return(TRUE)
  } else {
    stop(paste0("Cannot follow path (", test_path, ")."))
  }
}

# Parse Cactus seqFile
# Species and assembly paths
seq_cols <- c("Species", "Assembly")
if (!checkPath(seqFile)) stop(checkPath(seqFile))
seqs <- read.table(seqFile, sep = "\t", fill = NA, header = F,
                   # Omit phylogenetic tree line
                   comment.char = "(",
                   col.names = seq_cols)
# Species names (minus target species)
refs <- paste(grep(spc_int, seqs$Species, value = T, invert = T), collapse = ",")
# Reformat seqs for input
seqs <- seqs %>%
  mutate(Species=paste0(Species, ".fasta"))
# Phylogeny
t1 <- read.tree(seqFile)
tree <- write.tree(t1, tree.names = F)

# MAF file
maf_file <- paste0(cactus_prefix, ".maf")
if (!checkPath(seqFile)) stop(checkPath(seqFile))

# Construct Ragout-style recipe file
rcp <- data.frame(.target=spc_int,
                  .references=refs,
                  .tree=tree,
                  .maf=maf_file) %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column()
colnames(seqs) <- colnames(rcp)
rcp <- rbind(rcp, seqs)

# Write Ragout recipe file
write.table(rcp, rcp_file, quote = F, sep = " = ", col.names = F, row.names = F)

# Return recipe file to calling script
cat(rcp_file)
