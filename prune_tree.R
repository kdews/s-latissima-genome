## Initialization
# Load required packages
library(ape)
if (require(showtext, quietly = TRUE)) {
  showtext_auto()
}
# Set working dir
setwd("/scratch2/kdeweese/latissima/genome_stats")
# Set filename(s)
if(length(commandArgs(trailingOnly = TRUE)) > 0) {
  line_args <- commandArgs(trailingOnly = TRUE)
} else {
  line_args <- c("brown_algae_newick_tree.txt", "species.txt")
}
tree_file <- line_args[1]
species_file <- line_args[2]
# Split filename for output naming
tree_ext <- tools::file_ext(tree_file)
out_tree <- gsub(paste0("\\.", tree_ext), paste0("_pruned", "\\.", tree_ext),
                 tree_file)

## Data
# Import tree and species list
species <- read.table(species_file)[,1]
t <- read.tree(tree_file)
# Prune tree
tp <- keep.tip(t, species)
write.tree(tp, out_tree)

# Plot phylogenetic trees
# Label only species of interest in full phylo
idx <- grep(paste(species, collapse = "|"), t$tip.label, invert = T)
t$tip.label[idx] <- ""
# Shorten species names for labels
t$tip.label <- gsub("[a-z]+_", "._", t$tip.label)
tp$tip.label <- gsub("[a-z]+_", "._", tp$tip.label)
# Save plots to png
png("phylo_prune.png", units = "in", width = 8, height = 5, res = 96)
par(mfrow = c(2, 1))
par(oma = rep(0, 4))
par(mar = c(2, 0, 0, 0))
# Full
plot.phylo(t, show.tip.label = T)
title(sub = "Full phylogeny")
axisPhylo()
# Pruned
par(mar = c(2, 17, 0, 6))
title(sub = "Pruned phylogeny")
plot.phylo(tp, show.tip.label = T)
axisPhylo()
dev.off()
