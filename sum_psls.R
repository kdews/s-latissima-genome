# Clear environment
rm(list = ls())
# Required packages
library(Hmisc, quietly = T, warn.conflicts = F)
suppressPackageStartupMessages(library(tidyverse, quietly = T, warn.conflicts = F))
library(ggpubr, quietly = T, warn.conflicts = F)
suppressPackageStartupMessages(library(ggpmisc, quietly = T, warn.conflicts = F))
library(RColorBrewer, quietly = T)
library(BiocManager, quietly = T)

# Input
# Only take command line input if not running interactively
if (interactive()) {
  wd <- "/project/noujdine_61/kdeweese/latissima/genome_stats"
  setwd(wd)
  # Cactus seqFile
  seq_file <- "s-latissima-genome/s_lat_alignment.txt"
  # Species of interest
  spc_int <- "Saccharina_latissima"
  # Output directory
  outdir <- "s-latissima-genome/"
} else {
  line_args <- commandArgs(trailingOnly = T)
  seq_file <- line_args[1]
  spc_int <- line_args[2]
  outdir <- line_args[3]
}
# Output files
match_sums_file <- "align_sums.tsv"
max_match_file <- "max_matches.tsv"
align_report_file <- "alignment_report.tsv"
# Prepend output directory to file names (if it exists)
if (dir.exists(outdir)) {
  match_sums_file <- paste0(outdir, match_sums_file)
  max_match_file <- paste0(outdir, max_match_file)
  align_report_file <- paste0(outdir, align_report_file)
  
}

# Functions
# Set column names for PSL data frames
psl_col <- c("matches", "misMatches", "repMatches", "nCount", "qNumInsert",
             "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName",
             "qSize", "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd",
             "blockCount", "blockSizes", "qStarts", "tStarts")
# Total genome length of species of interest
spc_int_len <- 615545555
# Get PSL file list from seqFile
listFiles <- function(seq_file) {
  # Data wrangling
  # Import data frame of species and associated file names
  species_tab <- read.table(seq_file, sep = "\t", skip = 1)
  species <- species_tab$V1
  # Designate species of interest as query genome in all combinations
  comp_spc <- grep(spc_int, species, invert = T, value = T)
  spc_int_vs <- paste(spc_int, comp_spc, sep = "_vs_")
  # Keep all combinations of remaining species
  comp_spc_vs <- gtools::permutations(v = comp_spc, n = length(comp_spc), r = 2)
  comp_spc_vs <- paste(comp_spc_vs[,1], comp_spc_vs[,2], sep = "_vs_")
  species_vs <- c(spc_int_vs, comp_spc_vs)
  # Create list of PSL file names
  psl_files <- sapply(species_vs, grep, list.files(pattern = "\\.psl",
                                                   path = ".",
                                                   full.names = T), value = T)
  # Select out PSLs with new filtering
  psl_files <- grep("new", psl_files, value = T, invert = T)
  return(psl_files)
}
# Gets genome name from "genome_vs_genome" character named: query_vs_target
getGen <- function(df_name, gen_type) {
  if (gen_type == "query") {
    n <- 1
  } else if (gen_type == "target") {
    n <- 2
  } else {
    stop("Error: <gen_type> must be either 'query' or 'target'.")
  }
  gen_name <- unlist(strsplit(df_name, "_vs_"))[[n]]
  return(gen_name)
}
# Abbreviate species genus name
abbrevSpc <- function(spc) {
  spc <- unlist(strsplit(spc, "_| "))
  let1 <- substr(spc[1], 1, 1)
  spc[1] <- paste0(let1, ".")
  spc_a <- paste(spc, collapse = " ")
  return(spc_a)
}
# Fixes chromosome labels for plotting
fixChrom <- function(contigs) {
  contigs <-
    as.character(as.numeric(str_remove_all(str_remove_all(contigs, ".*_"),
                                           "[^0-9]")))
  return(contigs)
}
# Order given ID column by size column in data frame
sizeSort <- function(df, id_col, size_col) {
  x <- df %>% select(all_of(c(id_col, size_col))) %>%
    arrange(desc(pick(all_of(size_col)))) %>%
    pull(all_of(id_col)) %>%
    unique()
  return(x)
}
# Order PSL data frame by contig size by converting contig names to factors
orderPsl <- function(df_name, df_list) {
  df <- df_list[[df_name]]
  query <- getGen(df_name, "query")
  target <- getGen(df_name, "target")
  # Extract numeric information from IDs
  df <- df %>% mutate(qNum=fixChrom(qName), tNum=fixChrom(tName))
  # Order ID character vectors by size
  query_contigs <- sizeSort(df, "qName", "qSize")
  target_contigs <- sizeSort(df, "tName", "tSize")
  query_nums <- sizeSort(df, "qNum", "qSize")
  target_nums <- sizeSort(df, "tNum", "tSize")
  # Capture size rank of each target ID with index
  target_idx <- df %>% select(tNum, tSize) %>%
    arrange(desc(tSize)) %>%
    unique() %>%
    # Filter out artificial chromosomes for ranking
    filter(tNum != "0") %>%
    # Use row ids for index
    rowid_to_column(var="index") %>%
    select(tNum, index) %>%
    # Force any artificial chromosomes to 0 index
    rbind(data.frame(index = 0, tNum = "0")) %>%
    # Convert to named vector
    deframe()
  df <- df %>% mutate(
    # Add integer index column using named vector
    index=as.integer(target_idx[tNum]),
    # Convert all IDs to factors ordered by size
    qName=factor(qName, levels = query_contigs),
    tName=factor(tName, levels = target_contigs),
    qNum=factor(qNum, levels = query_nums),
    tNum=factor(tNum, levels = target_nums)
  ) %>%
    arrange(qName)
  return(df)
}
# Summarize PSL table by summing matches for each contig pair
sumPsl <- function(df_name, df_list) {
  df <- df_list[[df_name]]
  query <- getGen(df_name, "query")
  target <- getGen(df_name, "target")
  # Group by columns to keep (ID pairs)
  df_sum <- df %>% group_by(index, qName, tName, qNum, tNum, qSize, tSize) %>%
    # Sum matches per ID pair
    summarize(total_matches=as.numeric(sum(matches)), .groups = "drop") %>%
    # Add columns for query and target
    mutate(query=query, target=target, .before = 1) %>%
    # Calculate and sort by percentage of query covered by matches
    mutate(qPercent=total_matches/qSize) %>% arrange(qNum, desc(qPercent))
  return(df_sum)
}
# Collapse list of sumPsl data frames into one data frame
collapseSums <- function(psl_sums) {
  # Collapse list of data frames
  match_sums <- bind_rows(psl_sums, .id = "Alignment") %>%
    arrange(Alignment) %>% select(!Alignment)
  return(match_sums)
}
# Pull out maximal matching contigs, where most of query contig maps onto target
maxMatches <- function(match_sums) {
  max_matches <- match_sums %>% group_by(query, target, qNum) %>%
    filter(qPercent == max(qPercent)) %>%
    ungroup() %>% arrange(query, target, tNum)
  return(max_matches)
}
# Summarize statistics per species alignment
perSpecies <- function(df) {
  sum_df <- df %>% group_by(query, target) %>%
    summarize(n_and_sum=paste(n(), sum(qSize)*1e-6, sep = ", "),
              .groups = "drop") %>%
    pivot_wider(names_from = "target", values_from = "n_and_sum")
  return(sum_df)
}
# Summarize statistics per homolog with species of interest
perHomolog <- function(df) {
  # Subset species of interest
  match_lens <- df %>% filter(query == spc_int) %>%
    # Abbreviate species names
    rowwise() %>% mutate(Species = abbrevSpc(target))
    # # Filter out artificial chromosomes
    # mutate(qNum=as.character(qNum), tNum=as.character(tNum)) %>%
    # filter(tNum != "0")
  # Convert species column to ordered factor sorted by species relatedness
  match_lens_sum <- match_lens %>% group_by(index, tNum, tSize, Species) %>%
    # Sums per ID (i.e., per contig)
    summarize(sum_homolog=sum(qSize),
              sum_match=sum(total_matches),
              `n homologs`=n(), .groups = "drop") %>% 
    arrange(Species, desc(tSize)) 
  return(match_lens_sum)
}
# Report alignment statistics in table
alnReport <- function(match_lens_sum) {
  align_report <- match_lens_sum %>% 
    # Convert bp to Mb
    mutate_at(grep("sum|size", colnames(.), ignore.case = T, value = T),
              ~.*1e-6) %>%
    # Per species statistics
    group_by(Species) %>%
    summarize(`Average homologs mapped per chromosome`=
                paste(round(smean.sd(`n homologs`), 2), collapse = " ± "),
              `Maximum homologs mapped per chromosome`=max(`n homologs`),
              `Minimum homologs mapped per chromosome`=min(`n homologs`),
              `Total homologs mapped`=sum(`n homologs`),
              `Total length of homologs`=sum(sum_homolog),
              `Average exact matches (Mb) per chromosome`=
                paste(round(smean.sd(sum_match), 2), collapse = " ± "),
              `Average exact matches (%) per chromosome`=
                paste(round(smean.sd(sum_match/tSize*100), 2), collapse = " ± "),
              `Maximum exact matches (Mb) per chromosome`=max(sum_match),
              `Minimum exact matches (Mb) per chromosome`=min(sum_match),
              `Total exact matches (Mb)`=sum(sum_match),
              .groups = "keep") %>%
    # Coerce all values to characters
    ungroup() %>% mutate(across(everything(), as.character)) %>%
    # Transpose table
    pivot_longer(cols = !Species,
                 names_to = " ", values_to = "Value") %>%
    pivot_wider(names_from = "Species", values_from = "Value")
  # write.table(align_report, file = align_report_file,
  #             quote = F, row.names = F, sep = "\t")
  # print(paste("Table of alignment statistics in:", align_report_file))
  return(align_report)
}
# Plot number and length of scaffold matches vs. chromosome length
plotLens <- function(match_lens_sum) {
  # Clean up data frame for plotting
  # Order "Species" column by decreasing relatedness
  spc_order <- c("japonica", "pyrifera", "pinnatifida", "siliculosus")
  lvls <- unname(sapply(spc_order, grep, unique(match_lens$Species), value = T))
  match_lens_sum <- match_lens_sum %>%
    mutate(Species=factor(Species, levels = lvls),
           # Convert bp to Mb
           `Reference chromosome length (Mb)`=tSize*1e-6,
           `summed homolog lengths (Mb)`=sum_homolog*1e-6,
           `summed match lengths (Mb)`=sum_match*1e-6,
           # Change Species variable name for legend
           Reference=Species)
  # Function for individual plots
  lenPlot <- function(my_var, lab_pos) {
    p <- ggplot(data = match_lens_sum,
                mapping = aes(x = `Reference chromosome length (Mb)`,
                              y = .data[[my_var]], group = Reference,
                              col = Reference, fill = Reference)) +
      geom_point(alpha = 0.5) +
      stat_poly_line(formula = y~x+0, alpha = 0.1) +
      stat_poly_eq(formula = y~x+0, mapping = use_label(labels = c("R2", "eq")),
                   label.x = lab_pos[["x"]], label.y = lab_pos[["y"]]) +
      scale_x_continuous(breaks = pretty_breaks()) +
      scale_y_continuous(breaks = pretty_breaks()) +
      theme_bw() +
      theme(legend.text = element_text(face = "italic"))
    return(p)
  }
  # Function to annotate left of plot with italicized name of species of interest
  annotSpcInt <- function(p) {
    p <- annotate_figure(
      p,
      left = text_grob(abbrevSpc(spc_int), face = "italic", rot = 90)
    )
    return(p)
  }
  p1 <- lenPlot("summed homolog lengths (Mb)", c(x = "right", y = "bottom"))
  leg1 <- get_legend(p1)
  p1 <- annotSpcInt(p1 + theme(legend.position = "none"))
  p2 <- lenPlot("summed match lengths (Mb)", c(x = "right", y = "top"))
  p2 <- annotSpcInt(p2 + theme(legend.position = "none"))
  p <- ggarrange(p1, p2, nrow = 2, align = "hv",
                 legend.grob = leg1, legend = "right", labels = "AUTO")
  p <- setNames(list(p), spc_int)
  return(p)
}

# Use PSL summary file if available
if (file.exists(match_sums_file)) {
  match_sums <- read.table(match_sums_file, header = T, sep = "\t")
} else {
  # Wrangle raw data
  psl_files <- listFiles(seq_file)
  # Read PSL files into list of data frames
  psl_list_raw <- sapply(psl_files, read.table, col.names = psl_col, simplify = F)
  # Remove file extensions from list names
  names(psl_list_raw) <- gsub(".*ment_|.psl", "", names(psl_list_raw))
  # Convert scaffold names to size-ordered factors
  psl_list <- sapply(names(psl_list_raw), orderPsl, psl_list_raw, simplify = F)
  # Summarize by contig vs. contig of each syntenic comparison
  psl_sums <- sapply(names(psl_list), sumPsl, psl_list, simplify = F)
  match_sums <- collapseSums(psl_sums)
  write.table(match_sums, match_sums_file, quote = F, row.names = F, sep = "\t")
  print(paste("Table of match sums in:", match_sums_file))
}

# Select maximal matching contigs
match_sums_no_0 <- match_sums %>% filter(tNum != "0", qNum != "0")
max_matches <- maxMatches(match_sums)
max_matches_no_0 <- maxMatches(match_sums_no_0)

613375999/spc_int_len
# Summary for all alignments
(all_sum <- perSpecies(match_sums))
(all_sum_no_0 <- perSpecies(match_sums_no_0))
(max_sum <- perSpecies(max_matches))
(max_sum_no_0 <- perSpecies(max_matches_no_0))

# Summarize alignment statistics between species of interest and references
unique_n <- match_sums %>% filter(query == spc_int) %>%
  select(qName, qSize) %>% unique() %>% pull(qName) %>% length()
unique_len <- match_sums %>% filter(query == spc_int) %>%
  select(qName, qSize) %>% unique() %>% pull(qSize) %>% sum()
unique_len_perc <- unique_len/spc_int_len*100

# _sums or max_ doesn't matter here because of unique()
match_sums %>% filter(query == spc_int) %>%
  select(qName, qPercent) %>% unique() %>% mutate(qPercent=qPercent*100) %>%
  # group_by(qName) %>%
  summarize(`Average exact match per query (%)`=
              paste(round(smean.sd(qPercent), 2), collapse = " ± "),
            `Median exact match per query (%)`=median(qPercent),
            `Maximum exact match per query (%)`=max(qPercent),
            `Minimum exact match per query (%)`=min(qPercent))


match_lens_sum <- perHomolog(match_sums)
max_match_lens_sum <- perHomolog(max_matches)

# Save report of alignment statistics by reference species used
(align_report <- alnReport(match_lens_sum))
(max_align_report <- alnReport(max_match_lens_sum))

