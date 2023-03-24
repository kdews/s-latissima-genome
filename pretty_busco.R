## Initialization
# Load required packages
library(ggplot2)
wd <- "/scratch2/kdeweese/latissima/genome_stats/"
setwd(wd)
# line_args <- commandArgs(TRUE)
# wd <- line_args[1]
wd <- "busco_summaries/eukaryota_odb10"
lineage <- gsub(".*/", "", wd)
# Source BUSCO-generated R script
source(paste(wd, "busco_figure.R", sep = "/"))
# Read table of species names and FASTA file names
spec_names <- read.table("assemblies/pretty_names.txt", sep = "\t")
names(spec_names) <- c("species", "my_species")

## Data
# Strip lineage names from species
lvls <- gsub(paste0("_", lineage), "", levels(df$my_species))
df$my_species <- factor(gsub(paste0("_", lineage), "", df$my_species),
                        levels = lvls)
# Match FASTA file names to BUSCO species names for merge
spec_names$my_species <- gsub("\\..*", "", spec_names$my_species)
spec_names$my_species <- factor(spec_names$my_species, levels = lvls)
df <- merge(df, spec_names, sort = FALSE)
# Labels for figure legend
lbs <- c("Complete and single-copy", "Complete and duplicated",
         "Fragmented", "Missing")
cat_labs <- data.frame(label=factor(lbs, levels = rev(lbs)),
                       category=unique(df$category))
df <- merge(df, cat_labs, sort = FALSE)
# Order species factor by highest complete (S) percentage
ord_spc <- df[df$category == "S",]
ord_spc <- ord_spc[order(ord_spc$my_percentage), "species"]
df$species <- factor(df$species, levels = ord_spc)
# Format plot axis labels
cnt <- paste0("n=", total_buscos)
yl <- paste0("% BUSCOs", "\n", lineage, " (", cnt, ")")
# Filter out species
drop_spec <- c("U. pinnatifida 20", "F. vesiculosus",
               "C. okamuranus", "E. siliculosus")
df <- df[grep(paste(drop_spec, collapse = "|"), df$species, invert = TRUE),]

## Plot
busc_plot <- ggplot(data = df,
                    aes(x = species, y = my_percentage, fill = label)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(x = NULL, y = yl) +
  scale_fill_manual(values = rev(my_colors)) +
  theme(legend.position = "top",
        legend.title = element_blank()) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, reverse = TRUE)) +
  geom_text(data = df[df$my_values > 0,],
              aes(label = my_values),
            position = position_stack(vjust = 0.5))

ggsave(file = paste(lineage, "busco_plot.png", sep = "_"),
       plot = busc_plot, width = 7, height = 4)
