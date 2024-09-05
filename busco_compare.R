## Initialization
# Load required packages
suppressPackageStartupMessages(library(tidytree, quietly = T))
suppressPackageStartupMessages(library(treeio, quietly = T))
suppressPackageStartupMessages(library(tidyverse, quietly = T))
suppressPackageStartupMessages(library(ggtree, quietly = T))
# suppressPackageStartupMessages(library(aplot, quietly = T))
if (require(showtext, quietly = T)) {
  showtext_auto()
  if (interactive()) showtext_opts(dpi=100) else showtext_opts(dpi=300)
}

## Input
if (interactive()) {
  setwd("/project/noujdine_61/kdeweese/latissima/genome_stats")
  line_args <- c("busco_summaries/eukaryota_odb10",
                 "s-latissima-genome/species_table.txt",
                 "s-latissima-genome/")
} else if (length(commandArgs(trailingOnly = T)) == 3) {
  line_args <- commandArgs(trailingOnly = T)
} else {
  stop("3 positional arguments expected.")
}
wd <- line_args[1]
assembly_file <- line_args[2]
outdir <- line_args[3]
busco_script <- "busco_figure.R"
tree_script <- "prune_tree.R"
quast_script <- "quast.R"

# Output
# Split lineage from working directory
lineage <- unlist(strsplit(wd, "/"))[2]
extens <- c("png", "tiff")
busc_plot_file <- paste0("F2_busco_", lineage, ".", extens)

# Prepend working directory to BUSCO script (if exists)
if (dir.exists(wd)) {
  busco_script <- paste(wd, busco_script, sep = "/")
}
# Prepend output directory to script and plot filenames (if it exists)
if (dir.exists(outdir)) {
  tree_script <- paste(outdir, tree_script, sep = "/")
  quast_script <- paste(outdir, quast_script, sep = "/")
  busc_plot_file <- paste0(outdir, busc_plot_file)
}

## Functions
# Labels for figure legend
lbs <- c("S"="Complete and single-copy",
         "D"="Complete and duplicated",
         "F"="Fragmented",
         "M"="Missing")
color_names <- unname(lbs)
# Format species Latin name
formatSpc <- function(spc) {
  spc_f <- gsub("_", " ", spc)
  # Converts Ectocarpus siliculosus to Ectocarpus sp. Ec32
  spc_f <- gsub("Ectocarpus siliculosus", "Ectocarpus sp. Ec32", spc_f)
  spc_f <- str_wrap(spc_f, width = 20)
  return(spc_f)
}

## Data wrangling
sourced <- T
# Source phylogenetic tree-pruning R script
cat(paste("Sourcing", tree_script))
source(tree_script)
# Source QUAST output parsing R script
cat(paste("Sourcing", quast_script))
source(quast_script)
# Source BUSCO-generated R script
cat(paste("Sourcing", busco_script))
suppressWarnings(source(busco_script))
my_colors <- setNames(my_colors, color_names)
# Read table of species names and FASTA file names
species_tab <- read.table(assembly_file, sep = "\t", fill = NA, header = F)
colnames(species_tab)[1:2] <- c("Species", "Assembly")
# Match species table names to BUSCO IDs and tree tip labels
species_tab <- species_tab %>%
  mutate(BUSCO_id=paste(basename(tools::file_path_sans_ext(Assembly)),
                        lineage, sep = "_"),
         # Add underscores to match tree species names
         Species=gsub(" ", "_", Species))
tree_spc <- sapply(tip.label(tp), grep, species_tab$Species, value=T)
tree_spc <- setNames(names(tree_spc), unname(tree_spc))
species_dict <- species_tab %>%
  mutate(Species=tree_spc[Species]) %>%
  pull(Species, BUSCO_id)
df <- df %>%
  # Match BUSCO ID with species name
  mutate(Species=unname(species_dict[my_species]),
         # Add more detail to labels
         label=factor(lbs[category], levels = unname(lbs)))
# Label summary of values in each category for each species
my_labs <- df %>%
  pivot_wider(id_cols = Species,
              names_from = category,
              values_from = my_values) %>%
  mutate(my_labs = paste(paste0("C:", S+D),
                         paste0("[S:", S, ","),
                         paste0("D:", D, "],"),
                         paste0("F:", F, ","),
                         paste0("M:", M))) %>%
  pull(my_labs)

# Convert Species to factor ordered by tree tip labels
df <- df %>% mutate(Species=factor(Species, levels = tip.label(tp)))

df2 <- quast_df %>%
  mutate(Species=gsub(" ", "_", Species),
         Species=tree_spc[Species],
         Species=factor(Species, levels = tip.label(tp))) %>%
  pivot_wider(id_cols = Species,
              names_from = Stat,
              values_from = Assembly)

# Format plot titles and axis labels
lin_split <- unlist(strsplit(lineage, "_"))
lin <- tools::toTitleCase(lin_split[1])
odb_v <- toupper(lin_split[2])
x_lab1 <- paste(lin, odb_v, "BUSCOs")
x_lab2 <- paste0("(n = ", total_buscos, ")")
x_lab <- paste(x_lab1, x_lab2, sep = "\n")

# Plotting
# Bar plot of BUSCO score breakdown
busc_plot <- ggplot(data = df, mapping = aes(x = my_percentage, y = Species)) +
  geom_col(mapping = aes(fill = label),
           position = position_stack(reverse = T)) +
  scale_fill_manual(name = NULL, values = my_colors,
                    labels = scales::label_wrap(width = 15)) +
  # guides(fill = guide_legend(ncol = 2)) +
  # annotate("text", x = 1:length(my_labs), y = 3, label = my_labs,
  #          hjust = 0, size = 5) +
  scale_x_continuous(name = x_lab,
                     labels = scales::label_percent(scale = 1),
                     expand = expansion(mult = 0, add = c(0, 2))) +
  scale_y_discrete(labels = as_labeller(formatSpc)) +
  coord_cartesian(clip = "off") +
  theme_classic2() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
tree_plot <- ggtree(tp) +
  geom_tiplab(aes(label = f_labs), align = T, fontface = "italic") +
  geom_treescale() +
  hexpand(0.5) +
  scale_x_continuous(expand = expansion(mult = 0, add = 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15), add = 0)) +
  coord_cartesian(clip = "off")
n_plot <- ggplot(data = df2, mapping = aes(x = `# N's per 100 kbp`,
                                           y = Species)) +
  geom_col() +
  scale_x_continuous(label = scales::label_comma(), 
                     expand = expansion(mult = 0, add = c(0, 2))) +
  coord_cartesian(clip = "off") +
  theme_classic2() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
pp <- ggarrange(n_plot, busc_plot, nrow = 1, legend = "right",
                legend.grob = get_legend(busc_plot), align = "h")
p <- ggarrange(tree_plot, pp, widths = c(1, 2), nrow = 1, align = "h")
# Save plot with message to user
cat(paste("Saving pretty BUSCO plot of", lineage, "to", busc_plot_file))
h <- 5
w <- h*2.5
sapply(busc_plot_file, ggsave,
       plot = p, bg = "white", width = w, height = h,
       simplify = F)
