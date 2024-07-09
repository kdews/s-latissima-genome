# Clear environment
rm(list = ls())
# Required packages
library(Hmisc, quietly = T, warn.conflicts = F)
suppressPackageStartupMessages(library(tidyverse, quietly = T,
                                       warn.conflicts = F))

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
core_scafs_file <- "core_scafs.tsv"
max_match_file <- "max_matches.tsv"
lens_file <- "lens_by_chrom.tsv"
align_report_file <- "alignment_report.tsv"
# Prepend output directory to file names (if it exists)
if (dir.exists(outdir)) {
  match_sums_file <- paste0(outdir, match_sums_file)
  core_scafs_file <- paste0(outdir, core_scafs_file)
  max_match_file <- paste0(outdir, max_match_file)
  lens_file <- paste0(outdir, lens_file)
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
    mutate(qPercent=total_matches/qSize, tPercent=total_matches/tSize) %>% 
    arrange(qNum, desc(qPercent))
  return(df_sum)
}
# Collapse list of sumPsl data frames into one data frame
collapseSums <- function(psl_sums) {
  # Collapse list of data frames
  match_sums <- bind_rows(psl_sums, .id = "Alignment") %>%
    arrange(Alignment) %>% select(!Alignment)
  return(match_sums)
}
# Intersection of homologous scaffolds between species
getIntersect <- function(match_sums) {
  core_scafs <- match_sums %>% filter(query == spc_int) %>%
    select(target, qName, qSize) %>% unique %>%
    group_by(qName, qSize) %>% summarize(n=n(), .groups = "drop") %>%
    filter(n == 4) %>% select(qName, qSize)
  return(core_scafs)
}
# Pull out maximal matching contigs, where most of query contig maps onto target
maxMatches <- function(match_sums) {
  max_matches <- match_sums %>% group_by(query, target, qNum) %>%
    filter(qPercent == max(qPercent)) %>% ungroup %>%
    arrange(query, target, tNum)
  return(max_matches)
}
# Summarize statistics per species alignment
perSpecies <- function(max_matches) {
  max_matches_sum <- max_matches %>% group_by(query, target) %>%
    summarize(n_and_sum=paste(n(), sum(qSize)*1e-6, sep = ", "),
              .groups = "drop") %>%
    pivot_wider(names_from = "target", values_from = "n_and_sum")
  return(max_matches_sum)
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
# Function to calculate mean and standard deviation and paste with sep "±" 
calcMeanSd <- function(x) {
  x_ <- paste(round(smean.sd(x), 2), collapse = " ± ")
  return(x_)
}
# Function to dynamically summarize variable
makeSummary <- function(df, col_nm, new_col) {
  av_col <- paste("Average", new_col, "per chromosome")
  av_p_col <- paste("Average", gsub(" \\(.*\\)", "", new_col),
                    "(%) per chromosome")
  max_col <- paste("Maximum", new_col, "per chromosome")
  min_col <- paste("Minimum", new_col, "per chromosome")
  tot_col <- paste("Total", new_col)
  sum_df <- df %>% summarize({{av_col}}:=calcMeanSd({{col_nm}}),
                             {{av_p_col}}:=calcMeanSd({{col_nm}}/tSize*100),
                             {{max_col}}:=max({{col_nm}}),
                             {{min_col}}:=min({{col_nm}}),
                             {{tot_col}}:=sum({{col_nm}}),
                             .groups = "drop")
  if (!grepl("exact", new_col)) {
    sum_df <- sum_df[,grep("%", colnames(sum_df), invert = T, value = T)]
  }
  return(sum_df)
}
# Make report
makeReport <- function(df) {
  rpt <- cbind(makeSummary(df, `n homologs`, "n homologs mapped"),
               makeSummary(df, sum_homolog, "length of homologs (Mb)"),
               makeSummary(df, sum_match, "exact matches (Mb)"))
  # Remove any duplicated grouping columns
  rpt <- rpt[!duplicated(colnames(rpt))]
  return(rpt)
}
# Report alignment statistics in table
alnReport <- function(df) {
  df <- df %>% 
    # Convert bp to Mb
    mutate_at(grep("sum|size", colnames(.), ignore.case = T, value = T),
              ~.*1e-6)
  # Per-species statistics
  spc_report <- df %>% group_by(Species) %>% makeReport
  # Statistics across all species
  tot_report <- df %>% makeReport %>% mutate(Species="All", .before = 1)
  # Combine reports
  align_report <- rbind(spc_report, tot_report) %>%
    # Coerce all values to characters
    mutate(across(everything(), as.character)) %>%
    # Pivot longer for plotting
    pivot_longer(cols = !Species,
                 names_to = "Statistic", values_to = "Value")
  return(align_report)
}
# Write transposed alignment report to file
writeReport <- function(df) {
  # Transpose table
  df <- df %>% pivot_wider(names_from = "Species", values_from = "Value")
    # column_to_rownames(var = "Statistic")
  write.table(df, file = align_report_file, quote = F, sep = "\t",
              row.names = F)
  print(paste("Table of alignment statistics in:", align_report_file))
}

# Summarize input data
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
# Combine list of data frames into one labeled by species
match_sums <- collapseSums(psl_sums)
write.table(match_sums, match_sums_file, quote = F, row.names = F, sep = "\t")
print(paste("Table of match sums in:", match_sums_file))
# Save core intersecting scaffolds
core_scafs <- getIntersect(match_sums)
write.table(core_scafs, core_scafs_file, quote = F, row.names = F, sep = "\t")
print(paste("Table of core scaffolds:", core_scafs_file,
      "(n =", dim(core_scafs)[1], sum(core_scafs$qSize)*1e-6, "Mb)"))
# Select maximal matching contigs
max_matches <- maxMatches(match_sums)
print(paste("Table of maximal matches in:", max_match_file))
write.table(max_matches, max_match_file, quote = F, row.names = F, sep = "\t")
# Per chromosome view
max_match_lens_sum <- perHomolog(max_matches)
write.table(max_match_lens_sum, lens_file, quote = F, row.names = F, sep = "\t")
print(paste("Table of n and length aligned per chromosome in:", lens_file))

# Summary table
# Save report of alignment statistics by reference species used
max_align_report <- alnReport(max_match_lens_sum)
writeReport(max_align_report)

# Summary for all alignments
max_sum <- perSpecies(max_matches)
