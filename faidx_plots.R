# Clear environment
rm(list = ls())
# Required packages
library(tidyverse)
library(gridExtra)
library(ggpubr)
if (require(showtext, quietly = TRUE)) {
  showtext_auto()
  if (interactive()) showtext_opts(dpi = 100) else showtext_opts(dpi = 300)
}
# Input
wd <- "/scratch2/kdeweese/latissima/genome_stats"
setwd(wd)
species_file <- "species.txt"
species_tab <- read.table(species_file)
species <- species_tab$V1
fns <- paste0(species_tab$V2, ".fai")
for (i in 1:length(fns)) {
  print(species[i])
  idx_temp <- read.table(fns[i])
  idx_temp <- idx_temp[,1:2]
  idx_temp$Species <- species[i]
  if (i == 1) {
    idx <- idx_temp
  } else {
    idx <- rbind(idx, idx_temp)
  }
}

colnames(idx) <- c("Scaffold", "Length", "Species")
idx <- idx[idx$Length > 1000000,]

ggplot(data = idx) + 
  geom_density(mapping = aes(x = Length, fill = Species))

plot_list <- list()
for (i in 1:length(species)) {
  idx_plot <- idx[idx$Species == species[i],]
  plot_list[[i]] <- ggplot(data = idx_plot, mapping = aes(x = Length)) +
    geom_histogram(stat = "density", fill = "black") +
    geom_density(col = "lightblue", fill = NULL) +
    ggtitle(species[i]) +
    xlim(c(min(idx$Length), max(idx$Length)))
}

((
  p1 <- ggarrange(plotlist = plot_list, ncol = 2, nrow = 2)
))
