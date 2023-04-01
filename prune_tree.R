## Initialization
# Load required packages
library(ape)
# Set working dir
setwd("/scratch2/kdeweese/latissima/genome_stats")
# Set filename(s)
if(length(commandArgs(trailingOnly = TRUE)) > 0) {
  line_args <- commandArgs(trailingOnly = TRUE)
} else {
  line_args <- c("s-latissima-genome/s_lat_alignment.txt")
}
tree_file <- line_args[1]
# Import tree and species list
species <- read.table(tree_file, comment.char = "(")[,1]
t <- read.tree(tree_file)
# Prune tree
tp <- keep.tip(t, species)
trees <- list(t, tp)
t_list <- unlist(strsplit(write.tree(t), split = "\\(|,|;"))
t_list <- t_list[t_list != ""]
tp_list <- unlist(strsplit(write.tree(tp), split = "\\(|,|;"))
tp_list <- tp_list[tp_list != ""]
grep(paste(species, collapse = "|"), t_list, value = T)
tp_list



# par(cex = 1.2, cex.axis = 1.5)
par(mar = rep(3, 4))
par(oma = rep(0, 4))
for (tree in list(t, tp)) {
  plot(tree, align.tip.label = T)
  box(lty = 3)
  axisPhylo()
  edgelabels(round(tree$edge.length), bg = "white")
}

plot(trees, layout = 2, align.tip.label = T)
box(lty = 3)
axisPhylo()
# edgelabels(round(tree$edge.length), bg = "white")

