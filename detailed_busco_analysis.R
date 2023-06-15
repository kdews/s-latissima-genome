# Clear environment
rm(list = ls())
# Required packages
library(tidyverse)
if (require(showtext, quietly = TRUE)) {
  showtext_auto()
  if (interactive()) showtext_opts(dpi=100) else showtext_opts(dpi=300)
}
# Set working directory
setwd("/scratch2/kdeweese/latissima/genome_stats/")
# BUSCO status codes and colors
codes_busco <- c("Missing", "Fragmented", "Complete", "Duplicated")
colors_busco <- c("#F04442", "#F0E442", "#56B4E9", "#3492C7")

# Functions
# Remove given pattern with gsub()
removePattern <- function(x, pattern) {
  gsub(pattern, "", x)
}
# Create data frame of BUSCO gene status per species for each lineage dataset
parseBuscos <- function(busco_csvs, species, lins) {
  buscos <- list()
  for (lin in lins) {
    fns <- grep(lin, busco_csvs, value = T)
    for (i in 1:dim(species)[1]) {
      spc <- species[,1][i]
      spc_fn <- species[,2][i]
      fn <- grep(spc_fn, fns, value = T)
      tab <- read.delim(fn, skip = 2)
      # Filter duplicate gene status entries
      tab <- unique(tab[,1:2])
      colnames(tab) <- c("BUSCO", spc)
      # Convert status variable to ordered factor
      tab[,spc] <- factor(tab[,spc], levels = codes_busco)
      if (spc == species[,1][1]) df <- tab else df <- merge(df, tab)
    }
    total <- dim(species)[1]
    buscos[[lin]]  <- df %>%
      rowwise() %>%
      mutate(Species_Present = sum(across(!BUSCO) == "Complete")) %>%
      pivot_longer(cols = !c(BUSCO, Species_Present),
                   names_to = "Species",
                   values_to = "Status")
  }
  return(buscos)
}
# Plot heat map comparing BUSCO content across genomes
heatmapBuscos <- function(lin, buscos) {
  df <- buscos[[lin]]
  p <- heatmap(df)
  # p <- ggplot(df, aes(x = BUSCO, y = Species)) +
  #   geom_tile(aes(fill = Status)) +
  #   theme(panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         panel.background = element_blank(),
  #         axis.line = element_line(color = "black")) +
  #   scale_fill_manual(values = colors_busco) +
  #   coord_flip() +
  #   ggtitle(tools::toTitleCase(lin))
  return(p)
}


# Input data
busco_csvs <- system("ls */*/full_table.tsv", intern = T)
species_table <- "assemblies/pretty_names.txt"
lins <- unique(gsub(".*_", "", gsub("_odb.*", "", busco_csvs)))
# Extract species ids
species <- read.table(species_table, sep = "\t")
species[,2] <- sapply(species[,2], removePattern, "\\.fasta")
spcs <- species[,2]
# Create data frame for each lineage
buscos <- parseBuscos(busco_csvs, species, lins)
# Generate heat map of BUSCO status for each genome in each lineage gene set
hmaps <- lapply(lins, heatmapBuscos, buscos)
print(hmaps)


