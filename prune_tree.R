## Initialization
# Load required packages
# library(ape)
suppressPackageStartupMessages(library(tidytree, quietly = T))
suppressPackageStartupMessages(library(treeio, quietly = T))
suppressPackageStartupMessages(library(tidyverse, quietly = T))
suppressPackageStartupMessages(library(ggtree, quietly = T))
suppressPackageStartupMessages(library(ggpubr, quietly = T))
if (require(showtext, quietly = T)) {
  showtext_auto()
}

# Input
if(interactive()) {
  # Set working dir
  setwd("/project/noujdine_61/kdeweese/latissima/genome_stats")
  assembly_file <- "s-latissima-genome/species_table.txt"
  tree_file <- "https://ars.els-cdn.com/content/image/1-s2.0-S1055790319300892-mmc1.txt"
  seq_file <- "s-latissima-genome/s_lat_alignment.txt"
  outdir <- "s-latissima-genome/"
} else if(sourced == T) {
  tree_file <- "https://ars.els-cdn.com/content/image/1-s2.0-S1055790319300892-mmc1.txt"
  seq_file <- "s_lat_alignment.txt"
  if (dir.exists(outdir)) seq_file <- paste0(outdir, seq_file)
} else if(length(commandArgs(trailingOnly = T)) > 0) {
  line_args <- commandArgs(trailingOnly = T)
  assembly_file <- line_args[1]
  tree_file <- line_args[2]
  seq_file <- line_args[3]
  outdir <- line_args[4]
} else {
  stop("4 positional arguments expected.")
}

# Set filename(s)
png_file <- "phylo_prune.png"
# Split filename for output naming
tree_name <- tools::file_path_sans_ext(basename(tree_file))
tree_ext <- tools::file_ext(tree_file)
out_tree <- paste0(tree_name, "_pruned.", tree_ext)
# Append output directory to output filenames (if it exists)
if (dir.exists(outdir)) {
  png_file <- paste0(outdir, png_file)
  out_tree <- paste0(outdir, out_tree)
}

# Format species to match tree names
unformatSpc <- function(spc) {
  spc <- gsub(" ", "_", spc)
  return(spc)
}
# Format species Latin name
formatSpc <- function(spc) {
  spc_f <- gsub("_", " ", spc)
  # Converts Ectocarpus siliculosus to Ectocarpus sp. Ec32
  spc_f <- gsub("Ectocarpus siliculosus", "Ectocarpus sp. Ec32", spc_f)
  # spc_f <- str_wrap(spc_f, width = 20)
  return(spc_f)
}

## Data
# Import species list
species_tab <- read.table(assembly_file, sep = "\t", fill = NA, header = F)
colnames(species_tab) <- c("Species", "Assembly")
species_tab <- species_tab %>% select(Species, Assembly)
# Import tree
t <- read.tree(tree_file)
t <- as.treedata(t)
# Reformat tip labels and save original tree labels to dictionary
f_labs <- formatSpc(tip.label(t))
# Use reformatted labels to match to species
spc_match <- unlist(sapply(f_labs, grep, species_tab$Species, value = T))
# Error if not all species in tree
if (!(any(species_tab$Species %in% spc_match))) {
  stop("Error: not all species found in tree.")
}
labs_df <- tibble(labs=tip.label(t), f_labs=f_labs) %>%
  mutate(label=unformatSpc(f_labs),
         spc=spc_match[f_labs],
         highlight=case_when(!is.na(spc) ~ 1))
spc_dict <- labs_df %>%
  filter(!if_any(everything(), .fns = is.na)) %>%
  pull(label, spc)
# Create output species table
# (original species replaced with formatted tree labels)
species_tab2 <- species_tab %>%
  mutate(Species=spc_dict[Species])
# Relabel tree
tip_dict <- labs_df %>% pull(label, labs)
tip.label(t) <- unname(tip_dict[tip.label(t)])
t <- full_join(t, labs_df, by = "label")
# Prune tree
tp <- keep.tip(t, unname(spc_dict))
write.tree(get.tree(tp), out_tree)
print(paste("Pruned tree written to:", out_tree))
# Generate formatted seqFile for Cactus, i.e.,
# NEWICK tree
# name1 path1
# name2 path2
# ...
# nameN pathN
write.tree(get.tree(tp), seq_file)
write.table(species_tab2, seq_file,
            append = T, sep = "\t",
            col.names = F, row.names = F, quote = F)
print(paste("Cactus-formated seqFile written to:", seq_file))

## Plot phylogenetic trees
t1 <- ggtree(t) +
  hexpand(0.2) +
  geom_tiplab(aes(label = f_labs, color = factor(highlight)),
              fontface = "italic",
              show.legend = F) +
  scale_color_manual(values = "red", na.value = "black") +
  theme_tree2()
t2 <- ggtree(tp) +
  hexpand(0.2) +
  geom_tiplab(aes(label = f_labs), fontface = "italic") +
  theme_tree2()
p <- ggarrange(t1, t2, labels = "AUTO")
ggsave(filename = png_file, plot = p)
print(paste("Full and pruned phylogenetic trees plotted in:", png_file))
