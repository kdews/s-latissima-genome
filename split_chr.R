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

# Split chr0s on gaps
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

sja <- read.table("assemblies/old_genomes/GCA_000978595.1/GCA_000978595.1_SJ6.1_genomic.fna.fai")
ecto <- read.table("assemblies/old_genomes/GCA_000310025.1/ncbi_dataset/data/GCA_000310025.1/GCA_000310025.1_ASM31002v1_genomic.fna.fai")
min(sja$V2)
min(ecto$V2)