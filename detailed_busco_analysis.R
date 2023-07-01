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
  spcs <- unique(df$Species)
  n_lost <- sum(df$Conservation == "Lost")
  n_reg <- sum(df$Conservation == "Regained")
  tot_df <- data.frame(BUSCO = rep("Genes", length(spcs)),
                       Species = spcs,
                       Status = "Complete",
                       Conservation = "NA")
  df <- rbind(df[,c("BUSCO", "Species", "Status", "Conservation")],
              tot_df)
  df <- df %>%
    rowwise() %>%
    mutate(Num =
             case_when((BUSCO == "Genes" & Species == "S. latissima") ~ n_lost,
                       (BUSCO == "Genes" & Species == "S. latissima2") ~ n_reg,
                       TRUE ~ NA),
           Num = as.factor(Num))
  return(df)
}
# Extract hclust dendrogram across BUSCO genes
getDendro <- function(lin, busc_mat) {
  mat <- busc_mat[[lin]]
  h <- heatmap(mat, scale = "none", keep.dendro = T)
  hclust_mat <- as.hclust(h$Rowv)
  return(hclust_mat)
}
# Arrange BUSCO column by hclust matrix
arrangeBuscos <- function(lin, buscos, busc_mat, hclust_mats) {
  df <- buscos[[lin]]
  mat <- busc_mat[[lin]]
  hclust_mat <- hclust_mats[[lin]]
  hclust_order <- hclust_mat$order
  df$BUSCO <- factor(x = df$BUSCO,
                     levels = c(rownames(mat)[hclust_order], "Genes"),
                     ordered = T)
  return(df)
}
# Plot clustered ggplot2 heat map
heatmapBuscos <- function(lin, buscos) {
  ttl <- tools::toTitleCase(lin)
  df <- buscos[[lin]]
  p <- ggplot(mapping = aes(x = Species, y = BUSCO)) +
    geom_tile(data = df,
              mapping = aes(fill = Status)) +
    scale_fill_manual(values = colors_busco) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          # axis.text.x = element_text(size = 8),
          axis.text.y = element_blank(),
          axis.line = element_line(color = "black"),
          legend.position = "left") +
    ggtitle(ttl)
  return(p)
}
# Annotated regained genes with alpha aesthetic
regainHeatmap <- function(lin, buscos, hmaps) {
  df <- buscos[[lin]]
  p <- hmaps[[lin]]
  p2 <- p +
    aes(alpha = Conservation) +
    scale_alpha_manual(values = c(1, 1, 0.2)) +
    geom_text(data = df, mapping = aes(label = Num, color = Num),
              size = 10, alpha = 1,
              vjust = "inward", hjust = "inward") +
    scale_color_manual(values = c("red", "green", "white")) +
    guides(alpha = "none", color = "none")
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
buscos <- lapply(lins, arrangeBuscos, buscos, busc_mat, hclust_mats)
names(buscos) <- lins
hmaps <- lapply(lins, heatmapBuscos, buscos)
names(hmaps) <- lins
hmaps_r <- lapply(lins, regainHeatmap, buscos, hmaps)
dends <- lapply(hclust_mats, rowDendro)
fig <- heatArrange(hmaps, dends)
fig_r <- heatArrange(hmaps_r, dends)

euk_reg <- buscos$eukaryota[buscos$eukaryota$Conservation == "Regained",]$BUSCO
stra_reg <- buscos$stramenopiles[buscos$stramenopiles$Conservation == "Regained",]$BUSCO
sug_csvs <- grep("s_lat|Sugar", busco_csvs, value = T)
for (fn in sug_csvs) {
  tab <- read.delim(fn, skip = 2)
  if (any(tab$X..Busco.id %in% euk_reg)) {
    euk_tab <- tab[tab$X..Busco.id %in% euk_reg, ]
  }
  if (any(tab$X..Busco.id %in% stra_reg)) {
    stra_tab <- tab[tab$X..Busco.id %in% stra_reg, ]
  }
}
euk_tab[, c("Sequence")]
stra_tab[, c("Sequence", "Description")]

# Save figures to files
showtext_opts(dpi = 300)
# width/height ratio should be 2/3: 8 in/12 in seems good
ggsave("busco_heat_across_genomes.png", fig,
       width = 8, height = 12, units = "in", bg = "white")
ggsave("busco_heat_across_genomes_REGAIN.png", fig_r,
       width = 8, height = 12, units = "in", bg = "white")
