# Clear environment
rm(list = ls())
# Required packages
library(tidyverse)
# Set working directory
setwd("/scratch2/kdeweese/latissima/genome_stats/")
# Input
busco_csvs <- system("ls */*/full_table.tsv", intern = T)
species_table <- "assemblies/pretty_names.txt"
lins <- unique(gsub(".*_", "", gsub("_odb.*", "", busco_csvs)))
stat_codes <- c("Missing", "Fragmented", "Complete", "Duplicated")
  
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
    for (spc in spcs) {
      fn <- grep(spc, fns, value = T)
      tab <- read.delim(fn, skip = 2)
      # Filter duplicate gene status entries
      tab <- unique(tab[,1:2])
      # Convert status column to factors with 4 levels
      tab[,"Status"] <- factor(tab[,"Status"], levels = stat_codes)
      colnames(tab) <- c("BUSCO", spc)
      if (spc == spcs[1]) {
        df <- tab
      } else {
        df <- merge(df, tab)
      }
    }
    total <- length(species)
    # df <- df %>%
    #   rowwise() %>%
    #   mutate(Conservation = 
    #            ifelse(n_distinct(c_across(!BUSCO)) == 1,
    #                   "Conserved",
    #                   ifelse(n_distinct(c_across(!c(BUSCO, s_lat_genome))) == 1,
    #                          "Lost", "Unconserved")),
    #          Conservation =
    #            factor(Conservation,
    #                   levels = c("Unconserved", "Conserved", "Lost"))) %>%
    #   ungroup()
    df <- df %>%
      rowwise() %>%
      mutate(Species_Present = sum(across(!BUSCO) == "Complete")) %>%
      pivot_longer(cols = !c(BUSCO, Species_Present),
                   names_to = "Species",
                   values_to = "Status")
    buscos[[lin]] <- df
  }
  return(buscos)
}

# Input data
# Extract species ids
species <- read.table(species_table, sep = "\t")
species[,2] <- sapply(species[,2], removePattern, "\\.fasta")
spcs <- species[,2]
# Create data frame for each lineage
buscos <- parseBuscos(busco_csvs, spcs, lins)
head(buscos$eukaryota)


for (lin in lins) {
  df <- buscos[[lin]]
  df <- df[grep("sugar|s_lat", df[["Species"]], ignore.case = T),]
  p <- ggplot(df, aes(x = BUSCO, y = Species)) +
    geom_tile(aes(fill = Status)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black")) +
    coord_flip() +
    # scale_fill_brewer(palette = "Dark2") +
    # scale_color_brewer(palette = "Dark2") +
    # facet_grid(cols = vars(Species)) +
    ggtitle(tools::toTitleCase(lin))
  print(p)
}


