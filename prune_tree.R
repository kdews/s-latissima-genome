## Initialization
# Load required packages
library(ape)
if (require(showtext, quietly = T)) {
  showtext_auto()
}

# Input
if(interactive()) {
  # Set working dir
  setwd("/project/noujdine_61/kdeweese/latissima/genome_stats")
  line_args <- c(
    "s-latissima-genome/species_table.txt",
    "https://ars.els-cdn.com/content/image/1-s2.0-S1055790319300892-mmc1.txt",
    "s-latissima-genome/s_lat_alignment.txt",
    "s-latissima-genome/"
  )
} else if(length(commandArgs(trailingOnly = T)) > 0) {
  line_args <- commandArgs(trailingOnly = T)
} else {
  stop("4 positional arguments expected.")
}
species_file <- line_args[1]
tree_file <- line_args[2]
seq_file <- line_args[3]
outdir <- line_args[4]
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
spc_tab <- read.table(species_file, sep = "\t", fill = NA, header = F)
colnames(spc_tab) <- c("Species", "Assembly")
spc_tab <- spc_tab %>% select(Species, Assembly)
# Import tree
t <- read.tree(tree_file)
# Reformat tip labels and save original tree labels to dictionary
f_labs <- formatSpc(t$tip.label)
# Use reformatted labels to match to species
spc_match <- unlist(sapply(f_labs, grep, spc_tab$Species, value = T))
# Error if not all species in tree
if (!(any(spc_tab$Species %in% spc_match))) {
  stop("Error: not all species found in tree.")
}
labs_df <- tibble(labs=t$tip.label, f_labs=f_labs) %>%
  mutate(new_labs=unformatSpc(f_labs),
         spc=spc_match[f_labs])
# Create output species table (original species matched to tree labels)
spc_dict <- labs_df %>%
  filter(!if_any(everything(), .fns = is.na)) %>%
  pull(new_labs, spc)
spc_tab2 <- spc_tab %>% mutate(Species=spc_dict[Species])
# Relabel tree
tip_dict <- labs_df %>% pull(new_labs, labs)
t$tip.label <- tip_dict[t$tip.label]
# Prune tree
tp <- keep.tip(t, unname(spc_dict))
write.tree(tp, out_tree)
print(paste("Pruned tree written to:", out_tree))
# Generate formatted seqFile for Cactus, i.e.,
# NEWICK tree
# name1 path1
# name2 path2
# ...
# nameN pathN
write.tree(tp, seq_file)
write.table(spc_tab2, seq_file,
            append = T, sep = "\t",
            col.names = F, row.names = F, quote = F)
print(paste("Cactus-formated seqFile written to:", seq_file))

## Plot phylogenetic trees
trees <- list(t, tp)
class(trees) <- "multiPhylo"
png(png_file, units = "in", width = 8, height = 5, res = 96)
# Get x limits for each graph
xlims <- list()
ylims <- list()
for (i in 1:length(trees)) {
  tree <- trees[[i]]
  plt <- plot.phylo(tree, plot = F)
  xlims[[i]] <- plt$x.lim
  ylims[[i]] <- plt$y.lim
}
# Save plots to png
par(oma = rep(0, 4), mar = rep(1, 4))
layout(matrix(c(1:length(trees)), byrow = T, nrow = 2),
       heights = c(length(trees):1))
for (i in 1:length(trees)) {
  tree <- trees[[i]]
  xlim <- xlims[[1]]
  xlim[2] <- xlim[2]*0.8
  ylim <- ylims[[i]]
  # Label only species of interest in full phylo
  idx <- grep(spc_dict, tree$tip.label, invert = T)
  tree$tip.label[idx] <- ""
  # Shorten species names for labels
  tree$tip.label <- gsub("[a-z]+_", "._", tree$tip.label)
  if (i == 1) {
    ttl = "Phylogenetic tree pruning"
  } else {
    ttl = ""
  }
  plt <- plot.phylo(tree, show.tip.label = T,
                    x.lim = xlim, y.lim = c(0, ylim[2]+1),
                    main = ttl)
  # Place scale bar in lower right corner
  if (i == length(trees)) {
    xpos <- plt$x.lim[1] + diff(plt$x.lim)*0.5
    ypos <- plt$y.lim[1] + (ylim[1] - plt$y.lim[1])*0.5
    len <- round(mean(trees[[1]]$edge.length), 2)
    add.scale.bar(x = xpos, y = ypos, length = len)
  }
}
dev.off()
print(paste("Full and pruned phylogenetic trees plotted in:", png_file))
