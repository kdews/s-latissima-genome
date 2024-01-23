## Initialization
# Load required packages
library(tidyverse, quietly = TRUE)
if (require(showtext, quietly = TRUE)) {
  showtext_auto()
  if (interactive()) showtext_opts(dpi=100) else showtext_opts(dpi=300)
}
# Input
if (interactive()) {
  setwd("/scratch1/kdeweese/latissima/genome_stats")
  line_args <- c("busco_summaries/eukaryota_odb10",
                 "s-latissima-genome/species_table.txt",
                 "s-latissima-genome/")
} else if (length(commandArgs(trailingOnly = TRUE)) == 3) {
  line_args <- commandArgs(trailingOnly = TRUE)
} else {
  stop("3 positional arguments expected.")
}
wd <- line_args[1]
spec_file <- line_args[2]
outdir <- line_args[3]
# Split lineage from working directory
lineage <- unlist(strsplit(wd, "/"))[2]
busc_plot_file <- paste0("busco_", lineage, ".png")
# Append output directory to plot name (if it exists)
if (dir.exists(outdir)) busc_plot_file <- paste0(outdir, busc_plot_file)

## Data wrangling
# Source BUSCO-generated R script
suppressWarnings(source(paste(wd, "busco_figure.R", sep = "/")))
# Read table of species names and FASTA file names
spec_names <- read.table(spec_file, sep = "\t", fill = NA, header = F)[,1:2]
names(spec_names) <- c("Species", "my_species")
# Wrap species names for plot
spec_names$Species <- str_wrap(spec_names$Species, width = 20)
# Strip lineage names from species
lvls <- gsub(paste0("_", lineage), "", levels(df$my_species))
df$my_species <- factor(gsub(paste0("_", lineage), "", df$my_species),
                        levels = lvls)
# Match FASTA file names to BUSCO species names for merge
spec_names$my_species <- tools::file_path_sans_ext(basename(spec_names$my_species))
spec_names$my_species <- factor(spec_names$my_species, levels = lvls)
df <- merge(df, spec_names, sort = FALSE)
# Labels for figure legend
lbs <- c(" Complete (C) and single-copy (S)  ",
           " Complete (C) and duplicated (D)",
           " Fragmented (F)  ",
           " Missing (M)")
cat_labs <- data.frame(label=factor(lbs, levels = rev(lbs)),
                       category=levels(df$category))
df <- merge(df, cat_labs, sort = FALSE)
# Order species factor by highest complete (S) percentage
ord_spc <- df[df$category == "S",]
ord_spc <- ord_spc[order(ord_spc$my_percentage), "Species"]
df$Species <- factor(df$Species, levels = ord_spc)
df <- df[order(df$Species),]
# Format plot titles and axis labels
cnt <- paste0("n = ", total_buscos)
yl <- paste0("\n", "% BUSCOs", " (", cnt, ")")
ttl <- unlist(strsplit(lineage, "_"))
ttl[1] <- tools::toTitleCase(ttl[1])
ttl[2] <- paste0("(", ttl[2], ")")
ttl <- paste(ttl, collapse = " ")
# Label summary of values in each categeory for each species
vals <- as.data.frame(pivot_wider(df[,c("Species", "my_values", "category")],
                    names_from = category,
                    values_from = my_values))
vals <- vals %>% mutate(summ = paste(paste0("C:", S+D),
                                     paste0("[S:", S, ","),
                                     paste0("D:", D, "],"),
                                     paste0("F:", F, ","),
                                     paste0("M:", M)))
vals <- vals$summ

## Plot
# Bar plot of BUSCO score breakdown
busc_plot <- ggplot(data = df,
                    aes(x = Species, y = my_percentage,
                        fill = label)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(title = my_title, subtitle = ttl, x = NULL, y = yl) +
  scale_fill_manual(values = rev(my_colors)) +
  annotate("text", x = 1:length(vals), y = 3, label = vals, hjust = 0, size = 5) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        text = element_text(size = 15),
        axis.text = element_text(color = "black"),
        axis.text.y = element_text(face = "italic")) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, reverse = TRUE))
# Save plot with message to user
print(paste("Saving pretty BUSCO plot of", lineage, "to", busc_plot_file))
ggsave(file = busc_plot_file, plot = busc_plot,
       width = my_width, height = my_height, units = my_unit)
