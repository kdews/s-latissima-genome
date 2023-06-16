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
  busc_mat <- list()
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
      # Simplify dataframe into numeric matrix
      nums <- as.matrix(as.numeric(tab[,spc]))
      dimnames(nums) <- list(tab[,"BUSCO"], spc)
      if (spc == species[,1][1]) {
        df <- tab
        mat <- nums
        } else {
          df <- merge(df, tab)
          mat <- cbind(mat, nums)
        }
    }
    total <- dim(species)[1]
    busc_mat[[lin]] <- mat
    buscos[[lin]]  <- df %>%
      rowwise() %>%
      mutate(All_Species_Present = 
               sum(across(!BUSCO) == "Complete") +
               sum(across(!BUSCO) == "Duplicated"),
             Species_Present =
               sum(across(!c(All_Species_Present, BUSCO,
                             `S. latissima`, `S. latissima2`)) ==
                     "Complete") + 
               sum(across(!c(All_Species_Present, BUSCO,
                             `S. latissima`, `S. latissima2`)) ==
                     "Duplicated")) %>%
      pivot_longer(cols = !c(BUSCO, Species_Present),
                   names_to = "Species",
                   values_to = "Status")
  }
  return(list(buscos, busc_mat))
}
# Filter for conserved BUSCO genes lost and gained between genome versions
filtAnalysis <- function(df) {
  df <- df %>%
    rowwise() %>%
    mutate(Lost = case_when((Species == "S. latissima" &
                               Status == "Missing" |
                               Status == "Fragmented" &
                               Species_Present == 4) ~ 1,
                            TRUE ~ 0),
           Regained = case_when((Species == "S. latissima2" &
                                   Lost == 1 &
                                   Status == "Complete" |
                                   Status == "Duplicated") ~ 1,
                                TRUE ~ 0))
  regained <- df[df$Regained == 1, "BUSCO"]$BUSCO
  df <- df %>%
    rowwise() %>%
    mutate(if_else(BUSCO %in% regained,
                   paste0("<span style='color: green'>", BUSCO, "</span>"),
                   BUSCO))
  
  return(df)
}
# Plot unclustered ggplot2 heat map
ggheatmapBuscos <- function(lin, buscos) {
  df <- buscos[[lin]]
  p <- ggplot(df, aes(x = Species, y = reorder(BUSCO, All_Species_Present))) +
    geom_tile(aes(fill = Status)) +
    theme(axis.text.y = ggtext::element_markdown(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black")) +
    scale_fill_manual(values = colors_busco) +
    ylab("BUSCO Gene") +
    ggtitle(tools::toTitleCase(lin))
  return(p)
}
# Plot clustered heat map comparing BUSCO content across genomes
heatmapBuscos <- function(lin, busc_mat) {
  mat <- busc_mat[[lin]]
  ttl <- tools::toTitleCase(lin)
  heatmap(mat, scale = "none", col = colors_busco,
          main = ttl, ylab = "BUSCO Gene")
  legend(x = "right", legend = codes_busco, fill = colors_busco)
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
busco_list <- parseBuscos(busco_csvs, species, lins)
buscos <- busco_list[[1]]
buscos <- lapply(buscos, filtAnalysis)
busc_mat <- busco_list[[2]]

# Generate heat map of BUSCO status for each genome in each lineage gene set
hmaps <- lapply(lins, ggheatmapBuscos, buscos)
print(hmaps)
# print(lapply(lins, heatmapBuscos, busc_mat))


