# Clear environment
rm(list = ls())
# Required packages
library(tidyverse, quietly = T, warn.conflicts = F)
library(gridExtra, quietly = T, warn.conflicts = F)
library(ggpubr, quietly = T)
if (require(showtext, quietly = T)) {
  showtext_auto()
  if (interactive()) showtext_opts(dpi = 100) else showtext_opts(dpi = 300)
}

# Input
if (interactive()) {
  wd <- "/scratch2/kdeweese/latissima/genome_stats"
  setwd(wd)
  assembly_file <- "s-latissima-genome/species_table.txt"
  for_seqtk <- "s-latissima-genome/for_seqtk.txt"
  # Species of interest
  spc_int <- "Saccharina_latissima"
  # Output directory
  outdir <- "s-latissima-genome/"
} else {
  line_args <- commandArgs(trailingOnly = T)
  assembly_file <- line_args[1]
  for_seqtk <- line_args[2]
  # Species of interest
  spc_int <- line_args[3]
  # Output directory
  outdir <- line_args[4]
}
# Output plot filename
vio_plot <- "scaffold_sizes_violin.png"
# Append output directory to plot name (if it exists)
if (dir.exists(outdir)) vio_plot <- paste0(outdir, vio_plot)

# Functions
# Summarize index statistics for each assembly
sumDf <- function(df) {
  sum_df <- df %>%
    group_by(Species) %>%
    summarize(N50=Biostrings::N50(`Length (Mb)`),
              noverN50 = sum(`Length (Mb)` > N50),
              total = sum(`Length (Mb)`),
              n = n())
  return(sum_df)
}
# Create annotated violin plot of contig lengths by species
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
    ggtitle(ttl) +
    theme(axis.text.x = element_blank(),
          # axis.text.x = element_text(angle = 25),
          axis.title.x = element_blank())
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
species_tab <- read.table(assembly_file, sep = "\t", fill = NA, header = F)
colnames(species_tab) <- c("Species", "Assembly", "Annotation", "Info", "Info2")
spc_int <- gsub("_"," ", spc_int)
species <- species_tab$Species
fns <- paste0(species_tab$Assembly, ".fai")
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
  arrange(desc(`Length (Mb)`), .by_group = T)
# Summarize data frame
sum_idx <- sumDf(idx)
# Filter contigs/scaffolds by size (100 Mb > length > 4 Mb)
idx_filt <- idx %>%
  filter(`Length (Mb)` < 100) %>%
  filter(`Length (Mb)` >= 4)
# Summarize filtered data frame
sum_idx_filt <- sumDf(idx_filt)
# Calculate average filtered genome length of all other species
target_len <- round(mean(sum_idx_filt %>%
                        filter(!grepl(spc_int, Species)) %>%
                        pull(total)))
# Retrieve contigs from species of interest to approximate length of related species
n_contigs <- sum_idx %>% 
  filter(grepl(spc_int, Species)) %>%
  pull(total)
for (x in 1:n_contigs) {
  temp_len <- sum(idx %>%
                    filter(grepl(spc_int, Species)) %>%
                    slice_head(n = x) %>%
                    pull(`Length (Mb)`))
  if (x == 1) {
    spc_lens <- data.frame(n=c(x), total_len=c(temp_len))
  } else {
    spc_lens <- rbind(spc_lens, c(x, temp_len))
  }
  if (!exists("opt_len") & temp_len >= target_len) {
    opt_len <- temp_len
    opt_n <- x
  }
}
# Filter out smaller contigs (< 1Mb) from species of interest
retrieved_contigs <- idx %>% 
  filter(grepl(spc_int, Species)) %>%
  slice_head(n = opt_n) %>%
  filter(`Length (Mb)` >= 1)
# Combine retrieved contigs with filtered index
idx_filt <- unique(rbind(idx_filt, retrieved_contigs))

# Plot length distributions using violin plots
p1 <- violinPlot(idx, "Unfiltered", n50 = T)
p2 <- violinPlot(idx_filt, "Size filtered", n50 = T)
ps <- ggarrange(p1, p2, common.legend = T, legend = "bottom")
# Save plots
showtext_opts(dpi = 300)
ggsave(vio_plot, ps, width = 13, height = 6, bg = "white")
# Export filtered lists of contigs for each species
assembly_list <- c()
outlist <- c()
for (spc in unique(idx_filt$Species)) {
  contig_list <- idx_filt %>%
    filter(Species == spc) %>%
    pull(ID)
  fname <- species_tab %>%
    filter(Species == spc) %>%
    pull(Assembly)
  outfile <- paste0("chromosome_extract_",
                    basename(tools::file_path_sans_ext(fname)),
                    ".txt")
  write.table(contig_list, file = outfile,
              quote = F, col.names = F, row.names = F, sep = "\t")
  assembly_list <- append(assembly_list, fname)
  outlist <- append(outlist, outfile)
}
# Create data frame of assemblies and respective contig lists
df_for_subset <- data.frame(species=unique(idx_filt$Species),
                            assembly=assembly_list,
                            contig_list=outlist)
write.table(df_for_subset, file = for_seqtk,
            quote = F, col.names = F, row.names = F, sep = "\t")
print(paste("Table of results in:", for_seqtk))