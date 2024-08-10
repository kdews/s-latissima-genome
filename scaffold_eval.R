# Clear environment
rm(list = ls())
# Required packages
suppressPackageStartupMessages(library(tidyverse, quietly = T,
                                       warn.conflicts = F))
library(gridExtra, quietly = T, warn.conflicts = F)
library(ggpubr, quietly = T, warn.conflicts = F)
library(scales, quietly = T, warn.conflicts = F)
if (require(showtext, quietly = T, warn.conflicts = F)) {
  showtext_auto()
  if (interactive()) showtext_opts(dpi = 100) else showtext_opts(dpi = 300)
}

# Inputs
if (interactive()) {
  wd <- "/project/noujdine_61/kdeweese/latissima/genome_stats"
  setwd(wd)
  assembly_file <- "s-latissima-genome/species_table.txt"
  # Species of interest
  spc_int <- "Saccharina_latissima"
  # Outgroup species
  out_spc <- "Ectocarpus_sp."
  # Output directory
  outdir <- "s-latissima-genome/"
} else {
  line_args <- commandArgs(trailingOnly = T)
  assembly_file <- line_args[1]
  # Species of interest
  spc_int <- line_args[2]
  # Outgroup species
  out_spc <- line_args[3]
  # Output directory
  outdir <- line_args[4]
}
# Order for species factor
spc_order <- c("Ectocarpus", "pinnatifida", "pyrifera", "japonica", "latissima")
# Output plot filenames
vio_plot <- "scaffold_sizes_violin.png"
bar_plot <- "scaffold_sizes_bar.png"
# Prepend output directory to plot name (if it exists)
if (dir.exists(outdir)) {
  vio_plot <- paste0(outdir, vio_plot)
  bar_plot <- paste0(outdir, bar_plot)
}
# Functions
# Abbreviate species genus name
abbrevSpc <- function(spc) {
  spc <- unlist(strsplit(spc, "_| "))
  let1 <- substr(spc[1], 1, 1)
  spc[1] <- paste0(let1, ".")
  spc_a <- paste(spc, collapse = " ")
  return(spc_a)
}
# Format species Latin name
formatSpc <- function(spc) {
  spc_f <- gsub("_", " ", spc)
  # Converts Ectocarpus siliculosus to Ectocarpus sp. Ec32
  spc_f <- gsub("Ectocarpus siliculosus", "Ectocarpus sp. Ec32", spc_f)
  spc_f <- str_wrap(spc_f, width = 20)
  return(spc_f)
}
# Extract numbers from contig IDs for filtering
fixChrom <- function(contigs) {
  contigs <-
    as.character(as.numeric(str_remove_all(str_remove_all(contigs, ".*_"),
                                           "[^0-9]")))
  return(contigs)
}
# Find N_len-length gaps (N's) from FASTA file using Biostrings
findGaps <- function(fasta_file, N_len) {
  fasta <- readDNAStringSet(fasta_file)
  gap_ptn <- paste(rep("N", N_len), collapse = "")
  gap_matches <- vmatchPattern(gap_ptn, fasta)
  start_comp <- startIndex(gap_matches)
  fasta_gaps <- tibble(start_comp) %>%
    mutate(seqid_comp = names(fasta)) %>%
    unnest_longer(start_comp) %>%
    mutate(start_comp = as.numeric(start_comp),
           end_comp = start_comp + N_len,
           gap_length = N_len)
  # %>%
  #   # Convert genomic position columns from bp to Mb
  #   mutate_at(.vars = vars(grep("start|end|length", colnames(.), value = T)),
  #             .funs = ~ .x*1e-6)
  return(fasta_gaps)
}
# Summarize index statistics for each assembly
sumDf <- function(df) {
  sum_df <- df %>%
    # Convert "Length" column from bp to Mb
    mutate(`Length (Mb)` = Length*1e-6) %>%
    group_by(Species) %>%
    summarize(N50=Biostrings::N50(`Length (Mb)`),
              L50=sum(`Length (Mb)` > N50),
              total=sum(`Length (Mb)`),
              n=n())
  return(sum_df)
}
# Create annotated curve
filtCurve <- function(spc_lens, spc_int, opt_n, real_n) {
  opt_len <- spc_lens[opt_n, "total_len"]
  real_len <- spc_lens[real_n, "total_len"]
  ttl <- paste(spc_int, "contigs filtered to approximate target genome length")
  p0 <- ggplot(data = spc_lens, aes(x = n, y = total_len)) + 
    geom_line() +
    annotate(geom = "point", shape = 23, size = 3, fill = "magenta",
             x = opt_n, y = opt_len) +
    annotate(geom = "point", shape = 23, size = 3, fill = "blue",
             x = real_n, y = real_len) +
    annotate(geom = "label", x = opt_n*2, y = opt_len,
             label = paste("Target length =", round(opt_len), "Mb"),
             color = "magenta") +
    annotate(geom = "label", x = real_n*2.5, y = real_len,
             label = paste("Real length =", round(real_len), "Mb"),
             color = "blue") +
    labs(title = ttl, x = "Number of scaffolds", y = "Genome length (Mb)") +
    theme_light()
  return(p0)
}
# Create annotated violin plot of contig lengths by species
violinPlot <- function(idx, sum_idx, n50 = NULL) {
  # Positioning functions for annotations
  # Calculate maximum y in density distributions of each species
  find_max_y <- function(spc) idx %>% filter(Species == spc) %>% pull(Length) %>%
    log10() %>% density() %>% .$x %>% max()
  max_y <- 10^(max(sapply(unique(idx$Species), find_max_y)))
  # size_fun <- function(x) max(density(x)$x)*1.1
  # Calculate N50 from original size distribution, then convert to log10
  n50_line_fun <- function(x) log10(Biostrings::N50(10^x))
  n50_fun <- function(x) n50_line_fun(x)*1.08
  # Plot annotations
  # Total size and number of scaffolds+contigs
  size_lab <- paste(paste(round(sum_idx$total, digits = 1), "Mb"),
                    paste("n =", prettyNum(sum_idx$n, big.mark = ",")),
                    sep = "\n")
  # N50/L50 of each assembly (in Mb)
  n50_lab <- paste(paste("N50 =", round(sum_idx$N50, digits = 1), "Mb"),
                   paste("L50 =", sum_idx$L50),
                   sep = "\n")
  # Scale to use for violin and sina plots
  vio_scale <- "width"
  p <- ggplot(data = idx, mapping = aes(x = Species, y = Length)) +
    # Convert y-axis from bp to log10 bp scale
    scale_y_log10(labels = label_log()) +
    labs(y = "Length (bp)") +
    # Violin: equalize violin width between groups, lower alpha for points
    geom_violin(mapping = aes(col = Species, fill = Species), scale = vio_scale,
                trim = F, linewidth = 1, alpha = 0.4, show.legend = F) +
    # Points: shrink size and jitter width (not height) for visibility
    # geom_jitter(height = 0, width = 0.1, alpha = 0.05) + # size = 0.2, ) +
    # ggforce::geom_sina(mapping = aes(col = Species),
    #                    scale = vio_scale, alpha = 0.8, show.legend = F) +
                       # , alpha = Length
    geom_boxplot(width = 0.1) +
    # Label assembly total size and n()
    geom_text(data = sum_idx, mapping = aes(x = levels(Species), y = max_y*5),
             label = size_lab) +
    theme_linedraw() +
    # Italicize Latin species names
    theme(axis.text.x = element_text(face = "italic"))
  # Label N50/L50 on top of violins and mark N50 with dashed line on graph
  p_n50 <- p +
    stat_summary(mapping = aes(x = as.numeric(Species) - 0.45,
                               xend = as.numeric(Species) + 0.45),
                 geom = "segment", linetype = "dashed", alpha = 0.8,
                 fun = n50_line_fun) + 
    stat_summary(mapping = aes(x = as.numeric(Species) - 0.33), geom = "label",
                 label = n50_lab, fill = "white", size = 2.5, fun = n50_fun)
  # max_len <- as.integer(max(idx$Length))
  # if (max_len > 10) {
  #   low_rng <- p_n50 + coord_cartesian(ylim = c(0, 1)) +
  #     scale_y_continuous(breaks = scales::breaks_width(1)) +
  #     theme(axis.title.y.left = element_blank(),
  #           title = element_blank())
  #   hi_rng <- p_n50 + coord_cartesian(ylim = c(1, max_len*1.3)) +
  #     scale_y_continuous(breaks = scales::breaks_width(10)) +
  #     theme(axis.ticks.x.bottom = element_blank(),
  #           axis.text.x.bottom = element_blank(),
  #           axis.title = element_blank())
  #   p_n50 <- ggarrange(hi_rng, low_rng, nrow = 2, heights = c(1, 2),
  #                      align = "v", legend = "none")
  #   p_n50 <- annotate_figure(p_n50, left = "Length (Mb)")
  #   return(p_n50)
  # }
  # if (max_len > 40) {
  #   low_rng <- p_n50 + coord_cartesian(ylim = c(0, sec_max_len)) +
  #     theme(axis.title.y.left = element_blank(),
  #           title = element_blank())
  #   hi_rng <- p_n50 + coord_cartesian(ylim = c(max_len*0.9, max_len*1.3)) +
  #     scale_y_continuous(breaks = scales::breaks_width(10)) +
  #     theme(axis.ticks.x.bottom = element_blank(),
  #           axis.text.x.bottom = element_blank(),
  #           axis.title = element_blank())
  #   p_n50 <- ggarrange(hi_rng, low_rng, nrow = 2, heights = c(1, 4),
  #                      align = "v", legend = "none")
  #   p_n50 <- annotate_figure(p_n50, left = "Length (Mb)")
  # } else {
  #   p_n50 <- p_n50 + coord_cartesian(ylim = c(0, max_len*1.3))
  # }
  if (missing(n50)) {
    return(p)
  } else if (n50) {
    return(p_n50)
  } else {
    print("Error: unrecognized argument to 'annot' variable.")
  }
}
# Create annotated graphs of contig length distribution by species
distPlot <- function(idx, limit = 1e6) {
  idx <- idx %>%
    arrange(Species, desc(Length)) %>%
    group_by(Species) %>%
    mutate(n=n(),
           ID=1:unique(n),
           n50=Biostrings::N50(Length),
           high=case_when(Length >= n50 ~ "over_N50"),
           cutoff=case_when(Length >= limit ~ paste0(">", limit),
                            .default = paste0("<", limit))) %>%
    ungroup() %>%
    mutate(Species=factor(Species, levels = rev(levels(Species))))
  p <- ggplot(data = idx, mapping = aes(x = ID, y = Length)) +
    # geom_hline(aes(yintercept=n50), col = "red", lty = 2, ) +
    # geom_col(aes(col = high, fill = high), show.legend = F) +
    geom_hline(yintercept = limit, col = "red", lty = 2) +
    geom_col(aes(col = cutoff, fill = cutoff), show.legend = F) +
    scale_color_manual(name = "", aesthetics = c("color", "fill"),
                       values = c("lightgrey", "darkblue")) +
    # Convert length to log10 scale
    scale_y_log10(labels = label_log()) +
    # scale_y_log10(labels = label_log(), expand = c(0, 0)) +
    # scale_x_continuous(expand = c(0, 0)) +
    labs(x = "Scaffolds + contigs", y = "Length (bp)") +
    facet_wrap(~Species, scales = "free_x") +
    theme_bw() +
    # Italicize Latin species names
    theme(strip.text = element_text(face = "italic"),
          strip.background = element_blank())
  return(p)
}

# Read in data
species_tab <- read.table(assembly_file, sep = "\t", fill = NA, header = F)
colnames(species_tab) <- c("Species", "Assembly", "Annotation")
species_tab <- species_tab[,na.omit(colnames(species_tab))]
spc_int <- formatSpc(spc_int)
out_spc <- formatSpc(out_spc)
species <- species_tab$Species
species <- sapply(species, formatSpc, USE.NAMES = F)
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
# Clean up data frame of combined length indices
colnames(idx) <- c("ID", "Length", "Species")
spc_lvls <- sapply(spc_order, grep, unique(idx$Species), value = T)
idx <- idx %>%
  # Sort "Species" column with factor for plotting
  mutate(Species = factor(Species, levels = spc_lvls)) %>%
  select(ID, Length, Species) %>%
  group_by(Species) %>%
  arrange(desc(Length), .by_group = T) %>%
  # Convert ID characters into numeric
  mutate(ID_num = as.numeric(str_remove_all(str_remove_all(ID, ".*_"),
                                            "[^0-9]"))) %>%
  # Fix SDR scaffold number
  replace_na(list(ID_num = 29))
# Summarize data frame
sum_idx <- sumDf(idx)

# Plots
# Plot length distributions using violin plots
p_l <- violinPlot(idx, sum_idx, n50 = T)
p_d <- distPlot(idx, 1e6)
# Save plots
print(paste("Saving violin plot to:", vio_plot))
ggsave(vio_plot, p_l, bg = "white", width = 8, height = 7)
print(paste("Saving barplot plot to:", bar_plot))
ggsave(bar_plot, p_d, bg = "white", width = 10, height = 6)


compVersions <- function(version_list) {
  df <- sapply(version_list, read.table, simplify = F) %>%
    bind_rows(.id = "Species") %>%
    mutate(Species=factor(Species, levels = names(version_list))) %>%
    rename(ID=V1, Length=V2) %>%
    select(Species, ID, Length)
  sum_df <- sumDf(df)
  p <- violinPlot(df, sum_df, n50 = T)
  return(p)
}
ec <- list(
  less_2kb="assemblies/old_genomes/Ectsi_genome_InfEg2kb.fasta.fai",
  JGI="assemblies/old_genomes/Ectsil1_AssemblyScaffolds_Repeatmasked.fasta.fai",
  ORCAE_V1="assemblies/old_genomes/Ectsi_genome_V2_cleaned.tfa.fai",
  ORCAE_V2="assemblies/EctsiV2_genome.fasta.fai",
  split_V2="assemblies/EctsiV2_genome_split_artificial.fasta.fai"
)
sj <- list(
  JGI="assemblies/old_genomes/Sacja1_AssemblyScaffolds_Repeatmasked.fasta.fai",
  ORCAE="assemblies/SJ_v6_2_chromosome.fa.fai",
  ORCAE_split="assemblies/SJ_v6_2_chromosome_split_artificial.fa.fai"
)
p_ec <- compVersions(ec)
p_sj <- compVersions(sj)
