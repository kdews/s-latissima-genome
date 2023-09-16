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

# Functions
densOverlay <- function(idx) {
  p <- ggplot(data = idx, mapping = aes(x = `Scaffold Length (Mb)`)) + 
    geom_histogram(aes(y = after_stat(density), fill = Species)) +
    geom_density() +
    scale_fill_viridis_d() +
    scale_color_viridis_d() +
    facet_grid(~ Species) +
    theme(legend.position="none")
  return(p)
}

violinPlot <- function(idx, cutoff) {
  sum_idx <- idx %>%
    group_by(Species) %>%
    summarize(max=max(`Scaffold Length (Mb)`),
              quant=quantile(`Scaffold Length (Mb)`)[["75%"]],
              n = n())
  annot_labs <- c(paste("n =", sum_idx$n),
                     paste("max = ", round(sum_idx$max, digits = 1)),
                     paste("75% = ", round(sum_idx$quant, digits = 1)))
  annot_labs <- paste(annot_labs, collapse = "\n")
  print(sum_idx)
  print(annot_labs)
  y_pos <- max(idx$`Scaffold Length (Mb)`)*1.1
  p <- ggplot(data = idx,
              mapping = aes(x = Species, y = `Scaffold Length (Mb)`)) +
    geom_violin(mapping = aes(col = Species, fill = Species), alpha = 0.4) +
    # geom_boxplot() +
    geom_jitter(height = 0, width = 0.02, size = 0.8) +
    geom_hline(yintercept = cutoff) +
    annotate(geom = "text", x = sum_idx$Species, y = y_pos,
             label = paste(paste("n =", sum_idx$n),
                           paste("max = ", round(sum_idx$max, digits = 1)),
                           sep = "\n")) +
    scale_fill_viridis_d() +
    scale_color_viridis_d()
  return(p)
}

# Read in data
species_tab <- read.table(species_file)
species <- gsub("_"," ", species_tab$V1)
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
idx <- idx %>%
  mutate(`Scaffold Length (Mb)`= Length/1000000) %>%
  select(Scaffold, `Scaffold Length (Mb)`, Species) %>%
  arrange(`Scaffold Length (Mb)`)

# Plot distributions of scaffold lengths
p1 <- violinPlot(idx, 2)
idx <- idx[idx$Species != "Saccharina japonica",]
p2 <- violinPlot(idx, 2)
ggarrange(p1, p2, common.legend = T)

f <- quantile(idx$`Scaffold Length (Mb)`)
f[["100%"]]

