## Initialization
# Load required packages
library(scales, quietly = T, warn.conflicts = F)
suppressPackageStartupMessages(library(tidyverse, quietly = T,
                                       warn.conflicts = F))
if (require(showtext, quietly = T, warn.conflicts = F)) {
  showtext_auto()
  if (interactive()) showtext_opts(dpi = 100) else showtext_opts(dpi = 300)
}

# Input
# Only take command line input if not running interactively
if (interactive()) {
  setwd("/project/noujdine_61/kdeweese/latissima/genome_stats")
  assembly_file <- "s-latissima-genome/species_table.txt"
  # Output directory
  outdir <- "s-latissima-genome/"
}
n_content_file <- "N_content_per_100_kb.png"
# Prepend output directory to file names (if it exists)
if (dir.exists(outdir)) {
  n_content_file <- paste0(outdir, n_content_file)
}
# Order species by decreasing relatedness
spc_order <- c("latissima", "japonica", "pyrifera", "pinnatifida", "Ectocarpus")
# Format species Latin name
formatSpc <- function(spc) {
  spc_f <- gsub("_", " ", spc)
  # Converts Ectocarpus siliculosus to Ectocarpus sp. Ec32
  spc_f <- gsub("Ectocarpus siliculosus", "Ectocarpus sp. Ec32", spc_f)
  return(spc_f)
}

# Data wrangling
# Import data frame of species and associated file names
species_tab <- read.table(assembly_file, sep = "\t", fill = NA, header = F)
species <- species_tab[,1]
fasta_base <- basename(tools::file_path_sans_ext(species_tab[,2]))
# Create list QUAST output directory names
all_dirs <- list.dirs(path = ".", recursive = F, full.names = F)
quast_dirs <- paste0("quast_", fasta_base)
quast_dirs <- sapply(quast_dirs, grep, all_dirs, value = T)
quast_files <- sapply(quast_dirs, list.files, pattern = "^report.tsv$",
                      full.names = T)
quast_list <- sapply(quast_files, read.table, fill = NA, header = F,
                     skip = 1, sep = "\t", comment.char = "", quote = "",
                     col.names = c("Stat", "Assembly", "Assembly_broken"),
                     simplify = F)
quast_list <- setNames(quast_list, species)
quast_df <- quast_list %>% bind_rows(.id = "Species")
lvls <- unlist(unname(sapply(spc_order, grep, unique(quast_df$Species),
                             value = T)))
quast_df <- quast_df %>% mutate(Species=factor(Species, levels = rev(lvls)))
n_content <- quast_df %>% filter(grepl("^# N", Stat))
p_n <- ggplot(data = n_content,
       mapping = aes(x = Species,
                     y = Assembly,
                     fill = Species, col = Species)) +
  geom_col(show.legend = F) +
  scale_x_discrete(labels = as_labeller(
    function(x) str_wrap(formatSpc(x), width = 15)
  )) +
  ylab(unique(n_content$Stat)) +
  theme_classic() +
  theme(axis.text.x = element_text(face = "italic"),
        legend.text = element_text(face = "italic"))
ggsave(p_n, filename = n_content_file)
