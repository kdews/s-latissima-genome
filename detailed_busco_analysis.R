# Clear environment
rm(list = ls())
# Required packages
library(tidyverse)
library(ggdendro)
library(gridExtra)
library(ggpubr)
if (require(showtext, quietly = TRUE)) {
  showtext_auto()
  if (interactive()) showtext_opts(dpi = 100) else showtext_opts(dpi = 300)
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
      pivot_longer(cols = !c(BUSCO, All_Species_Present, Species_Present),
                   names_to = "Species",
                   values_to = "Status")
  }
  return(list(buscos, busc_mat))
}
# Filter for conserved BUSCO genes lost and gained between genome versions
filtAnalysis <- function(df) {
  df <- df %>%
    rowwise() %>%
    mutate(Lost = case_when((Species == "S. latissima" & Species_Present == 4 & 
                               Status == "Missing") ~ 1,
                            (Species == "S. latissima" & Species_Present == 4 &
                               Status == "Fragmented") ~ 1,
                            TRUE ~ 0),
           Lost = as.factor(Lost))
  lost <- df[df$Lost == 1, "BUSCO"]$BUSCO
  df <- df %>%
    rowwise() %>%
    mutate(Regained = case_when((Species == "S. latissima2" & BUSCO %in% lost &
                                   Status == "Complete") ~ 1,
                                (Species == "S. latissima2" & BUSCO %in% lost &
                                   Status == "Duplicated") ~ 1,
                                TRUE ~ 0),
           Regained = as.factor(Regained))
  df <- df %>%
    rowwise() %>%
    mutate(Conservation = case_when(Lost == 1 ~ "Lost",
                                    Regained == 1 ~ "Regained",
                                    TRUE ~ "NA"),
           Conservation = factor(Conservation,
                                 levels = c("Regained", "Lost", "NA")))
  return(df)
}
# Extract hclust dendrogram across BUSCO genes
getDendro <- function(lin, busc_mat) {
  mat <- busc_mat[[lin]]
  h <- heatmap(mat, scale = "none", keep.dendro = T)
  hclust_mat <- as.hclust(h$Rowv)
  return(hclust_mat)
}
# Plot clustered ggplot2 heat map
heatmapBuscos <- function(lin, buscos, busc_mat, hclust_mats) {
  ttl <- tools::toTitleCase(lin)
  df <- buscos[[lin]]
  mat <- busc_mat[[lin]]
  hclust_mat <- hclust_mats[[lin]]
  hclust_order <- hclust_mat$order
  df$BUSCO <- factor(x = df$BUSCO,
                     levels = rownames(mat)[hclust_order],
                     ordered = T)
  p <- ggplot(df, aes(x = Species, y = BUSCO)) +
    geom_tile(aes(fill = Status)) +
    scale_fill_manual(values = colors_busco) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_blank(),
          axis.line = element_line(color = "black"),
          legend.position = "left") +
    ggtitle(ttl)
  return(p)
}
# Annotated regained genes with alpha aesthetic
regainHeatmap <- function(p) {
  p2 <- p + 
    aes(alpha = Conservation) +
    scale_alpha_manual(values = c(1, 1, 0.2)) + 
    guides(alpha = "none")
  return(p2)
}
# Create row dendrogram for BUSCO gene axis
rowDendro <- function(hclust_mat) {
  d <- ggdendrogram(hclust_mat, rotate = T) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())
  return(d)
}
# Arrange heatmap and dendrogram plots
heatArrange <- function(hmaps, dends) {
  p1 <- ggarrange(hmaps[[1]], hmaps[[2]], 
                  ncol = 1, nrow = 2,
                  common.legend = T, legend = "left")
  p2 <- ggarrange(dends[[1]], dends[[2]], 
                  ncol = 1, nrow = 2)
  p3 <- ggarrange(p1, p2, ncol = 2, nrow = 1, widths = c(1, 0.25))
  return(p3)
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
busc_mat <- busco_list[[2]]
buscos <- lapply(buscos, filtAnalysis)

# Generate heat map of BUSCO status for each genome in each lineage gene set
hclust_mats <- lapply(lins, getDendro, busc_mat)
names(hclust_mats) <- lins
hmaps <- lapply(lins, heatmapBuscos, buscos, busc_mat, hclust_mats)
hmaps_r <- lapply(hmaps, regainHeatmap)
dends <- lapply(hclust_mats, rowDendro)
fig <- heatArrange(hmaps, dends)
fig_r <- heatArrange(hmaps_r, dends)
# Save figures to files
showtext_opts(dpi = 300)
ggsave("busco_heat_across_genomes.png", fig, 
       width = 7, height = 10, units = "in", bg = "white")
ggsave("busco_heat_across_genomes_REGAIN.png", fig_r, 
       width = 7, height = 10, units = "in", bg = "white")
