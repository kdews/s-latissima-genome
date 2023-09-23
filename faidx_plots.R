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
sumDf <- function(df) {
  sum_df <- df %>%
    group_by(Species) %>%
    summarize(N50=Biostrings::N50(`Length (Mb)`),
              noverN50 = sum(`Length (Mb)` > N50),
              total = sum(`Length (Mb)`),
              n = n())
  return(sum_df)
}

violinPlot <- function(idx, ttl, n50 = NULL) {
  sum_idx <- sumDf(idx)
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
    scale_color_viridis_d() +
    ggtitle(ttl)
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
# Create grouped data frame and convert "Length" column
colnames(idx) <- c("ID", "Length", "Species")
idx <- idx %>%
  mutate(`Length (Mb)`= Length/1000000) %>%
  select(ID, `Length (Mb)`, Species) %>%
  group_by(Species) %>%
  arrange(desc(`Length (Mb)`), .by_group = T) %>%
  filter(Species != "Saccharina japonica")
# Summarize data frame
sum_idx <- sumDf(idx)
# Filter contigs/scaffolds by size (>4 Mb)
idx_filt <- idx %>%
  filter(`Length (Mb)` >= 4)
# Summarize filtered data frame
sum_idx_filt <- sumDf(idx_filt)
target_len <- round(mean(sum_idx_filt %>%
                        filter(Species != "Saccharina latissima") %>%
                        pull(total)))
# Approximate genome length of related species by retrieving more S. lat contigs
n_contigs <- sum_idx %>% 
  filter(Species == "Saccharina latissima") %>%
  pull(total)
for (x in 1:n_contigs) {
  temp_len <- sum(idx %>%
                    filter(Species == "Saccharina latissima") %>%
                    slice_head(n = x) %>%
                    pull(`Length (Mb)`))
  if (x == 1) {
    sac_lens <- data.frame(n=c(x), total_len=c(temp_len))
  } else {
    sac_lens <- rbind(sac_lens, c(x, temp_len))
  }
  if (!exists("opt_len") & temp_len >= target_len) {
    opt_len <- temp_len
    opt_n <- x
  }
}
# Filter out smaller S. latissima contigs (< 1Mb)
retrieved_contigs <- idx %>% 
  filter(Species == "Saccharina latissima") %>%
  slice_head(n = opt_n) %>%
  filter(`Length (Mb)` >= 1)
# Combine retrieved contigs with filtered index
idx_filt <- unique(rbind(idx_filt, retrieved_contigs))

# Plot length distributions using violin plots
p1 <- violinPlot(idx, "Unfiltered", n50 = T)
p2 <- violinPlot(idx_filt, "Size filtered")
ps <- ggarrange(p1, p2, common.legend = T)

# Save plots
showtext_opts(dpi = 300)
ggsave("scaffold_sizes_violin.png", ps, width = 13, height = 6, bg = "white")
