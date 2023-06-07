# Clear environment
rm(list = ls())
# Set working directory
setwd("/scratch2/kdeweese/latissima/genome_stats/")
# Input files
busco_csvs <- list.files(pattern = "full_table\\.tsv$", recursive = T)
species_table <- "assemblies/pretty_names.txt"
lins <- unique(gsub(".*_", "", gsub("_odb.*", "", busco_csvs)))

# Functions
# Removes given pattern with gsub()
removePattern <- function(x, pattern) {
  gsub(pattern, "", x)
}

# Input data
# Extract species ids
species <- read.table(species_table, sep = "\t")
species[2] <- sapply(species[2], removePattern, "\\.fasta")
# Create data frame for each lineage
buscos <- list()
for (spc in species[[2]]) {
  fns <- grep(spc, busco_csvs, value = T)
  lin_busc <- list()
  for (lin in lins) {
    fn <- grep(lin, fns, value = T)
    tab <- read.delim(fn, skip = 2)
    tab <- tab[,1:3]
    lin_busc[[lin]] <- tab
  }
  buscos[[spc]] <- lin_busc
}

unique(unlist(lapply(lapply(buscos, "[[", lins[1]), "[[", 1)))
attributes(buscos)
names(buscos[grep("s_lat_genome|Sugarkelp", buscos)])

