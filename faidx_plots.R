# Clear environment
rm(list = ls())
# Required packages
suppressPackageStartupMessages(library(tidyverse, quietly = T, warn.conflicts = F))
library(gridExtra, quietly = T, warn.conflicts = F)
library(ggpubr, quietly = T)
suppressPackageStartupMessages(library(Biostrings, quietly = T, warn.conflicts = F))
if (require(showtext, quietly = T)) {
  showtext_auto()
  if (interactive()) showtext_opts(dpi = 100) else showtext_opts(dpi = 300)
}

# Input
if (interactive()) {
  wd <- "/project/noujdine_61/kdeweese/latissima/genome_stats"
  setwd(wd)
  assembly_file <- "s-latissima-genome/species_table.txt"
  for_seqtk <- "s-latissima-genome/for_seqtk.txt"
  # Species of interest
  spc_int <- "Saccharina_latissima"
  # Outgroup species
  out_spc <- "Ectocarpus_siliculosus"
  # Output directory
  outdir <- "s-latissima-genome/"
} else {
  line_args <- commandArgs(trailingOnly = T)
  assembly_file <- line_args[1]
  for_seqtk <- line_args[2]
  # Species of interest
  spc_int <- line_args[3]
  # Outgroup species
  out_spc <- line_args[4]
  # Output directory
  outdir <- line_args[5]
}
# Output plot filenames
filt_plot <- paste0(spc_int, "_filtering.png")
vio_plot <- "scaffold_sizes_violin.png"
# Prepend output directory to plot name (if it exists)
if (dir.exists(outdir)) {
  filt_plot <- paste0(outdir, filt_plot)
  vio_plot <- paste0(outdir, vio_plot)
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
    group_by(Species) %>%
    summarize(N50=Biostrings::N50(`Length (Mb)`),
              noverN50 = sum(`Length (Mb)` > N50),
              total = sum(`Length (Mb)`),
              n = n())
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
violinPlot <- function(idx, ttl, n50 = NULL) {
  sum_idx <- sumDf(idx)
  max_len <- as.integer(max(idx$`Length (Mb)`))
  sec_max_len <- as.integer(sort(idx$`Length (Mb)`, decreasing = T)[2])
  p <- ggplot(data = idx,
              mapping = aes(x = Species, y = `Length (Mb)`)) +
    geom_violin(mapping = aes(col = Species, fill = Species),
                alpha = 0.4, linewidth = 1) +
    geom_jitter(height = 0, width = 0.02, size = 0.8) +
    annotate(geom = "text", x = sum_idx$Species, y = max_len+3,
             label = paste(round(sum_idx$total, digits = 1), "Mb\n",
                           "n =", sum_idx$n)) +
    scale_fill_viridis_d() +
    scale_color_viridis_d() +
    ggtitle(ttl) +
    # Remove legend
    theme(legend.position = "none")
  n50_lab <- paste(
    # paste(round(sum_idx$total, digits = 1), "Mb"),
    # paste("n =", sum_idx$n),
    paste("N50 =", round(sum_idx$N50, digits = 1), "Mb"),
    paste("n>N50 =", sum_idx$noverN50),
    sep = "\n")
  p_n50 <- p +
    # Mark N50 with white triangle
    stat_summary(fun = Biostrings::N50,
                 geom = "point", col = "black", alpha = 0.5,
                 # col = "black", fill = "white", size = 3, shape = 25,
                 size = 10, shape = "—") +
    # Label N50
    stat_summary(fun = Biostrings::N50, label = n50_lab,
                 geom = "label", col = "black", size = 3,
                 mapping = aes(group = Species, y = max_len-3))
                 # position = position_nudge(x = 0.35, y = 3),
  if (max_len > 40) {
    low_rng <- p_n50 +
      coord_cartesian(ylim = c(0, sec_max_len)) +
      theme(axis.title.y.left = element_blank(),
            title = element_blank())
    hi_rng <- p_n50 +
      coord_cartesian(ylim = c(max_len-5, max_len+5)) +
      scale_y_continuous(breaks = scales::breaks_width(10)) +
      theme(axis.ticks.x.bottom = element_blank(),
            axis.text.x.bottom = element_blank(),
            axis.title = element_blank())
    p_n50 <- ggarrange(hi_rng, low_rng, nrow = 2, heights = c(1, 4),
                       align = "v", legend = "none")
    p_n50 <- annotate_figure(p_n50, left = "Length (Mb)")
  }
  if (missing(n50)) {
    return(p)
  } else if (n50) {
    p_n50 <- list(p_n50, get_legend(p))
    return(p_n50)
  } else {
    print("Error: unrecognized argument to 'annot' variable.")
  }
}
# Create annotated graphs of contig length distribution by species
distPlot <- function(spc, idx) {
  idx <- idx %>%
    mutate(`Length (log10 bp)` = log10(`Length (Mb)`*1e6))
  y_min <- min(idx$`Length (Mb)`)
  y_max <- max(idx$`Length (Mb)`)
  idx_filt <- idx %>%
    filter(Species == spc) %>%
    arrange(desc(`Length (log10 bp)`)) %>%
    mutate(ID_num = factor(ID_num, levels = ID_num))
  idx_high <- idx_filt %>%
    filter(`Length (log10 bp)` >= Biostrings::N50(`Length (log10 bp)`))
  p <- ggplot(data = idx_filt,
              mapping = aes(x = ID_num, y = `Length (log10 bp)`)) +
    geom_col(col = "lightgrey", fill = "lightgrey") +
    geom_hline(yintercept = Biostrings::N50(idx_filt$`Length (log10 bp)`),
               col = "red", lty = 2) +
    geom_col(data = idx_high,
             mapping = aes(x = ID_num, y = `Length (log10 bp)`),
             col = "darkblue", fill = "darkblue") +
    annotate(geom = "text", col = "white",
             label = paste0("n = ", dim(idx_high)[1]),
             x = dim(idx_high)[1]/2,
             y = mean(idx$`Length (log10 bp)`)) +
    coord_cartesian(ylim = c(y_min, log10(35*1e6))) +
    labs(title = spc) +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  return(p)
}

# Read in data
species_tab <- read.table(assembly_file, sep = "\t", fill = NA, header = F)
colnames(species_tab) <- c("Species", "Assembly", "Annotation")
species_tab <- species_tab[,na.omit(colnames(species_tab))]
# spc_int <- gsub("_"," ", spc_int)
# out_spc <- gsub("_"," ", out_spc)
spc_int <- abbrevSpc(spc_int)
out_spc <- abbrevSpc(out_spc)
species <- species_tab$Species
# Keep only Latin species names and abbreviate
species <- word(species, start = 1, end = 2, sep = " ")
species <- sapply(species, abbrevSpc, USE.NAMES = F)
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
# Create grouped data frame and convert "Length" column from bp to Mb
colnames(idx) <- c("ID", "Length", "Species")
idx <- idx %>%
  mutate(`Length (Mb)`= Length/1000000) %>%
  select(ID, `Length (Mb)`, Species) %>%
  group_by(Species) %>%
  arrange(desc(`Length (Mb)`), .by_group = T) %>%
  mutate(ID_num = as.numeric(str_remove_all(str_remove_all(ID, ".*_"),
                                            "[^0-9]"))) %>%
  replace_na(list(ID_num = 29))
# Summarize data frame
sum_idx <- sumDf(idx)
# Filter contigs/scaffolds by size (100 Mb > length > 4 Mb)
idx_filt <- idx %>%
  # Removes ORCAE artificial chromosomes
  filter(ID_num != 0) %>%
  filter(`Length (Mb)` < 100) %>%
  filter(`Length (Mb)` >= 4)
# Summarize filtered data frame
sum_idx_filt <- sumDf(idx_filt)
# Retrieve more scaffolds from species of interest
# Calculate average filtered genome length of all other species, minus outgroup
target_len <- round(mean(sum_idx_filt %>%
                        filter(!grepl(spc_int, Species)) %>%
                        filter(!grepl(out_spc, Species)) %>%
                        pull(total)))
# Retrieve contigs from species of interest to approximate length of related species
n_contigs <- sum_idx %>% 
  filter(grepl(spc_int, Species)) %>% 
  pull(n)
  # pull(total)
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
# Save actual n_contigs of species of interest for plotting
real_n <- dim(retrieved_contigs)[1]
# Combine retrieved contigs with filtered index
idx_filt <- unique(rbind(idx_filt, retrieved_contigs))
sum_idx_filt_2 <- sumDf(idx_filt)

faidx <- read.table("assemblies/split_scaff_test/SJ_v6_2_chromosome_chr0_split.fa.fai")
test <- findGaps("assemblies/split_scaff_test/SJ_v6_2_chromosome_chr0.fa", 200)
test <- test %>% mutate(contig_len = c(start_comp[[1]]-1, diff(start_comp)-200)) %>%
  filter(contig_len > 1) %>%
  rowid_to_column(var = "index") %>%
  # filter(contig_len %in% faidx$V2) %>%
  select(index, contig_len)
faidx <- faidx %>% rowid_to_column(var = "index") %>%
  # filter(V2 %in% test$contig_len) %>%
  select(index, V2)

fa_mer <- merge(faidx, test, by.x = "V2", by.y = "contig_len", sort = F,
                # all.x = T,
                suffixes = c(".faidx", ".200bp"))
fa_mer <- fa_mer %>% mutate(diff_index = index.faidx - index.200bp) %>%
  arrange(index.faidx)
fa_mer_filt <- fa_mer %>%
  filter(diff_index >= 0) %>%
  filter(diff_index < 1000)
ggplot(data = fa_mer_filt, mapping = aes(x = index.faidx, y = index.200bp)) +
  geom_point() +
  theme_classic()

old_fasta_file <- "assemblies/old_genomes/GCA_000978595.1/GCA_000978595.1_SJ6.1_genomic.fna"
old_fasta <- readDNAStringSet(old_fasta_file)
old_gaps <- letterFrequency(old_fasta, letters = "N")
length(which(old_gaps==0)) + length(which(old_gaps>0))
which(old_gaps==200)
gap_ptn <- DNAString(paste(rep("N", 200), collapse = ""))
old_gap_matches <- names(vmatchPattern(gap_ptn, old_fasta))
length(names(old_fasta))
old_fasta$`JXRI01000608.1 Saccharina japonica cultivar Ja scaffold609, whole genome shotgun sequence`
old_gaps_scaf <- old_gaps[old_gaps>0]
median(old_gaps_scaf)
median(old_gaps)
fasta_file <- "assemblies/split_scaff_test/SJ_v6_2_chromosome_chr0.fa"
fasta <- readDNAStringSet(fasta_file)
letterFrequency(fasta, letters = "N")
gap_ptn <- DNAString(paste(rep("N", 200), collapse = ""))
gap_matches <- matchPattern(gap_ptn, fasta$chr0)
strsplit(fasta, gap_ptn)
ggplot(data = test, mapping = aes(x = start_comp, y = contig_len)) +
  geom_col()
# Plots
# Plot filtering of species of interest on curve of length vs. n_contigs
p0 <- filtCurve(spc_lens, spc_int, opt_n, real_n)
# Plot length distributions using violin plots
((p_l <- violinPlot(idx, "", n50 = T)))
idx_no0 <- idx %>% filter(fixChrom(ID) != "0")
violinPlot(idx_no0, "", n50 = T)
p1_l <- violinPlot(idx, "Unfiltered", n50 = T)
p2_l <- violinPlot(idx_filt, "Size filtered", n50 = T)
ps <- ggarrange(p1_l[[1]], p2_l[[1]], legend.grob = p2_l[[2]], legend = "right")
p_list <- lapply(species, distPlot, idx)
p3 <- ggarrange(plotlist = p_list, align = "v")
# Save plots
showtext_opts(dpi = 300)
ggsave(filt_plot, p0, width = 10, height = 6)
ggsave(vio_plot, ps, width = 20, height = 6, bg = "white")
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