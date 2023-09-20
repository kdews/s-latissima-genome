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
violinPlot <- function(idx, n50 = NULL) {
  sum_idx <- idx %>%
    group_by(Species) %>%
    summarize(N50=Biostrings::N50(`Length (Mb)`),
              noverN50 = sum(`Length (Mb)` > N50),
              total = sum(`Length (Mb)`),
              n = n())
  # annot_labs <- paste(paste("n =", sum_idx$n),
  #                     paste("n>N50 =", sum_idx$noverN50),
  #                     sep = "\n")
  p <- ggplot(data = idx,
              mapping = aes(x = Species, y = `Length (Mb)`)) +
    geom_violin(mapping = aes(col = Species, fill = Species),
                alpha = 0.4, linewidth = 1) +
    geom_jitter(height = 0, width = 0.02, size = 0.8) +
    annotate(geom = "text", x = sum_idx$Species, y = max(idx$`Length (Mb)`)*1.1,
             label = paste("n =", sum_idx$n)) +
    annotate(geom = "text", x = sum_idx$Species, y = max(idx$`Length (Mb)`)*1.05,
             label = paste(round(sum_idx$total, digits = 1), "Mb")) +
    scale_fill_viridis_d() +
    scale_color_viridis_d()
  n50_lab <- paste(paste("N50 =", round(sum_idx$N50, digits = 1), "Mb"),
                   paste("n>N50 =", sum_idx$noverN50),
                   sep = "\n")
  p_n50 <- p +
    stat_summary(fun = Biostrings::N50,
                 geom = "point", col = "black", fill = "white",
                 size = 3, shape = 25) +
    stat_summary(fun = Biostrings::N50, label = n50_lab,
                 geom = "label", col = "black", size = 3,
                 position = position_nudge(x = 0.35, y = 2))
  if (missing(n50)) {
    return(p)
  } else if (n50) {
    return(p_n50)
  } else {
    print("Error: unrecognized argument to 'annot' variable.")
  }
  
}

# Read in data
species_tab <- read.table(species_file)
species <- gsub("_"," ", species_tab$V1)
fns <- paste0(species_tab$V2, ".fai")
for (i in 1:length(fns)) {
  idx_temp <- read.table(fns[i])
  idx_temp <- idx_temp[,1:2]
  idx_temp$Species <- species[i]
  if (i == 1) {
    idx <- idx_temp
  } else {
    idx <- rbind(idx, idx_temp)
  }
}
colnames(idx) <- c("ID", "Length", "Species")
idx <- idx %>%
  mutate(`Length (Mb)`= Length/1000000) %>%
  select(ID, `Length (Mb)`, Species) %>%
  group_by(Species) %>%
  arrange(desc(`Length (Mb)`), .by_group = T)

# Plot length distributions using violin plots
p1 <- violinPlot(idx, n50 = T)
idx <- idx %>%
  filter(Species != "Saccharina japonica")
p2 <- violinPlot(idx, n50 = T)
idx_sac <- idx %>%
  filter(Species == "Saccharina latissima") %>%
  mutate(Species = paste0(Species, " (", gsub("_.*", "", ID), ")"))
p3 <- violinPlot(idx_sac)
idx <- idx %>%
  filter(`Length (Mb)` >= 4)
p4 <- violinPlot(idx, n50 = T)
ggarrange(p1, p2, p3, p4, common.legend = T)
print(idx_sac, n = 21)
# print(idx %>% filter(`Length (Mb)` > 4), n = 100)

# print(idx[idx$Species == "Macrocystis pyrifera",], n=34)
# print(idx[idx$Species == "Undaria pinnatifida",], n=50)
# print(idx[idx$Species == "Saccharina latissima",], n=50)

# showtext_opts(dpi = 300)
# ggsave("scaffold_sizes_violin.png", p2, width = 8, height = 6)
