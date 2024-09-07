# Clear environment
rm(list = ls())
# Required packages
suppressPackageStartupMessages(library(tidyverse, quietly = T, warn.conflicts = F))
# library(ggrepel)
suppressPackageStartupMessages(library(ggpmisc, quietly = T, warn.conflicts = F))
library(ape, quietly = T, warn.conflicts = F)
library(viridisLite, quietly = T, warn.conflicts = F)
library(RColorBrewer, quietly = T, warn.conflicts = F)
library(gridExtra, quietly = T, warn.conflicts = F)
library(ggpubr, quietly = T, warn.conflicts = F)
library(ggpmisc, quietly = T, warn.conflicts = F)
suppressPackageStartupMessages(library(Biostrings, quietly = T, warn.conflicts = F))
if (require(showtext, quietly = T, warn.conflicts = F)) {
  showtext_auto()
  if (interactive())
    showtext_opts(dpi = 100)
  else
    showtext_opts(dpi = 300)
}

# Input
if (interactive()) {
  setwd("/project/noujdine_61/kdeweese/latissima/genome_stats")
  # Table of species in analysis
  seqFile <- "s-latissima-genome/s_lat_alignment.txt"
  # Ragout output directory
  ragout_dir <- "ragout-out-all-solid-refine-update"
  # Species of interest
  spc_int <- "Saccharina_latissima"
  # Output directory
  outdir <- "s-latissima-genome/"
} else {
  line_args <- commandArgs(trailingOnly = T)
  # Table of species in analysis
  seqFile <- line_args[1]
  # Ragout output directory
  ragout_dir <- line_args[2]
  # Species of interest
  spc_int <- line_args[3]
  # Output directory
  outdir <- line_args[4]
}

# Output plot file names
extens <- c("png", "tiff")
outfiles <- list(
  len_plot = "F4A_ragout_length_dist",
  map_plot = "FS3_ragout_pseudo_mapping",
  comp_len_plot = "FS4_ragout_comp_length_dist",
  comp_map_plot = "ragout_comp_pseudo_mapping"
)
outfiles <- sapply(outfiles, paste, extens, sep = ".", simplify = F)

# If exists, prepend output directory to output file names
redirect <- function(filename, outdir) {
  if (dir.exists(outdir))
    filename <- paste0(outdir, filename)
  return(filename)
}
outfiles <- sapply(outfiles, redirect, outdir, simplify = F)

# Global variables
types <-
  c("input",
    "pseudochromosomes",
    "chimera",
    "remainder",
    "excluded")
label_names <- c("input", "gap")

# Functions
# Checks existence of file or directory and errors if FALSE
checkPath <- function(test_path) {
  if (file.exists(test_path)) {
    return(TRUE)
  } else {
    stop(paste0("Cannot follow path (", test_path, ")."))
  }
}
# Fixes chromosome labels for plotting
fixChrom <- function(contigs) {
  contigs <-
    as.character(as.numeric(str_remove_all(
      str_remove_all(contigs, ".*_"),
      "[^0-9]"
    )))
  return(contigs)
}
# Function to extract title from directory name
extTtl <- function(ragout_dir) {
  ttl_dict <- list(
    all = "No filtering,",
    `no-sac` = "Short contigs filtered out (except S. latissima),",
    filt = "Short contigs filtered out,",
    solid = "unbroken scaffolds",
    chimera = "broken scaffolds"
  )
  ttl <- gsub("ragout-out-", "", ragout_dir)
  vals <- list()
  for (val in names(ttl_dict))
    if (grepl(val, ttl))
      vals <- append(vals, val)
  for (val in vals)
    ttl <- gsub(val, ttl_dict[[val]], ttl)
  ttl <- gsub("-refine", "", ttl)
  ttl <- gsub("-", " ", ttl)
  return(ttl)
}
# Read Cactus-formatted seqFile
readSpecies <- function(seqFile) {
  seq_cols <- c("Species", "Assembly")
  checkPath(seqFile)
  seqs <- read.table(
    seqFile,
    sep = "\t",
    fill = NA,
    header = F,
    # Omit phylogenetic tree line
    comment.char = "(",
    col.names = seq_cols
  )
  colnames(seqs) <- seq_cols[1:length(seqs)]
  seqs <- seqs[, na.omit(colnames(seqs))]
  return(seqs)
}
# Find 10-kb N gaps from FASTA file using Biostrings
findGaps <- function(fasta_file, def_len) {
  # gap_len <- 1e4
  fasta <- readDNAStringSet(fasta_file)
  Nfreq <- as_tibble(alphabetFrequency(fasta)) %>%
    mutate(seqid_comp = names(fasta)) %>%
    filter(!grepl("^chr0\\.|^chr_00\\.", seqid_comp)) %>%
    pull(N, seqid_comp)
  gap_len <- min(Nfreq[Nfreq > 1])
  if (gap_len > 2e4) {
    gap_len <- def_len
  }
  print(paste("Gap size:", gap_len, "bp"))
  gap_ptn <- paste(rep("N", gap_len), collapse = "")
  gap_matches <- vmatchPattern(gap_ptn, fasta)
  start_comp <- startIndex(gap_matches)
  fasta_gaps <- tibble(start_comp) %>%
    mutate(seqid_comp = names(fasta),
           Length = width(fasta),) %>%
    unnest_longer(start_comp) %>%
    mutate(
      start_comp = as.numeric(start_comp),
      end_comp = start_comp + gap_len,
      gap_length = gap_len,
    ) %>%
    filter(!grepl("^chr0\\.|^chr_00\\.", seqid_comp)) %>%
    group_by(seqid_comp, Length) %>%
    group_modify( ~ add_row(
      .x,
      start_comp = 0,
      end_comp = 0,
      gap_length = 0,
      .before = 1
    )) %>%
    arrange(desc(Length), start_comp) %>%
    mutate(Contig_Size = case_when(
      is.na(lead(start_comp)) ~ Length - end_comp,
      .default = lead(start_comp) - end_comp
    )) %>%
    filter(Contig_Size > 1) %>%
    mutate(
      # Account for start = 0 rows added in each group
      Number_of_Gaps = n() - 1,
      Gap_Pattern = Number_of_Gaps * gap_len,
      N_Frequency = Nfreq[seqid_comp]
    ) %>%
    ungroup() %>%
    mutate(seqid_comp = factor(seqid_comp, levels = unique(seqid_comp)))
  # %>%
  #   mutate(start_comp = as.numeric(start_comp),
  #          end_comp = start_comp + gap_len,
  #          gap_length = gap_len,
  #          label = factor("gap", levels = label_names, ordered = T)) %>%
  #   # Convert genomic position columns from bp to Mb
  #   mutate_at(.vars = vars(grep("start|end|length", colnames(.), value = T)),
  #             .funs = ~ .x*1e-6)
  return(fasta_gaps)
}
# Find gap size of each reference
findGapSize <- function(def_len, named_fas) {
  all_gaps <- sapply(named_fas, findGaps, def_len, simplify = F)
  df_gaps <- all_gaps %>%
    bind_rows(.id = "Species")
  p <- ggplot(
    data = df_gaps,
    mapping = aes(
      x = N_Frequency,
      y = Gap_Pattern,
      color = Species,
      group = Species
    )
  ) +
    geom_point() +
    stat_poly_eq(mapping = use_label("eq"),
                 formula = y ~ x + 0) +
    stat_poly_line(fullrange = T,
                   se = T,
                   formula = y ~ x + 0) +
    labs(title = def_len) +
    theme_bw()
  return(p)
}
# Import NCBI AGP file (v2.0)
readAgp <- function(ragout_dir, pre_gaps) {
  agp_cols <- c(
    "seqid_scaf",
    "start_scaf",
    "end_scaf",
    "component_number",
    "component_type",
    "seqid_comp",
    "start_comp",
    "end_comp",
    "orientation"
  )
  # Named vector for renaming to gap colnames
  gap_cols <-
    c(agp_cols[1:5], "gap_length", "gap_type", "linkage", "evidence")
  names(gap_cols) <- agp_cols
  # AGP file describing scaffolding
  checkPath(ragout_dir)
  agp_file <-
    list.files(path = ragout_dir,
               pattern = ".*agp$",
               full.names = T)
  agp_raw <- read.table(agp_file, col.names = agp_cols)
  agp <- agp_raw %>%
    filter(component_type == "W") %>%
    select(!contains("component")) %>%
    mutate(label = factor("input", levels = label_names, ordered = T))
  gaps_df <- agp_raw %>%
    rename_with(
      .fn = function(x)
        gap_cols[x]
    ) %>%
    filter(component_type == "N") %>%
    select(!contains(c("component", "linkage", "type"))) %>%
    # Quantify number of species used as evidence for gap
    rowwise() %>%
    mutate(n_evidence = length(unlist(strsplit(evidence, ",")))) %>%
    ungroup() %>%
    mutate(label = factor("gap", levels = label_names, ordered = T))
  agp_combo <- merge(agp, gaps_df, all = T)
  agp_combo <- agp_combo %>%
    # Convert genomic position columns to numeric type and bp to Mb
    mutate_at(.vars = vars(grep(
      "start|end|length", colnames(.), value = T
    )),
    .funs = ~ as.numeric(.x) * 1e-6)
  seqid_df <- agp_combo %>%
    select(seqid_scaf,
           seqid_comp,
           start_scaf,
           end_scaf,
           start_comp,
           end_comp,
           orientation) %>%
    filter_all(all_vars(!is.na(.))) %>%
    unique()
  pre_gaps <- merge(
    pre_gaps,
    seqid_df,
    sort = F,
    by = "seqid_comp",
    suffixes = c("", "_OG")
  )
  pre_gaps <- pre_gaps %>%
    # Keep only v1.0 gaps if they are included in rescaffolded assembly
    filter(start_comp >= start_comp_OG &
             end_comp <= end_comp_OG) %>%
    mutate(
      start_scaf = case_when(
        orientation == "+" ~
          start_scaf + (start_comp - start_comp_OG),
        orientation == "-" ~
          start_scaf + (end_comp_OG - end_comp),
        .default = NA
      ),
      end_scaf = start_scaf + gap_length
    ) %>%
    select(!contains("OG"))
  agp_combo <- merge(agp_combo, pre_gaps, all = T, sort = F)
  agp_combo <- agp_combo %>% arrange(seqid_scaf, start_scaf)
  agp_gaps <-
    agp_combo %>% filter(label == "gap" & !is.na(seqid_comp))
  for (i in 1:dim(agp_gaps)[1]) {
    agp_gaps_row <- agp_gaps[i,]
    row_id <- which(
      agp_combo$label == "input" &
        agp_combo$seqid_comp == agp_gaps_row$seqid_comp &
        agp_combo$start_scaf <= agp_gaps_row$start_scaf &
        agp_combo$end_scaf >= agp_gaps_row$end_scaf
    )
    agp_row1 <- agp_combo[row_id,] %>%
      mutate(end_scaf = agp_gaps_row$start_scaf)
    agp_row2 <- agp_combo[row_id,] %>%
      mutate(start_scaf = agp_gaps_row$end_scaf)
    agp_combo <- agp_combo %>%
      add_row(agp_row2)
    agp_combo[row_id,] <- agp_row1
  }
  agp_combo <- agp_combo %>%
    group_by(seqid_scaf) %>%
    # Sum N content of each scaffold
    mutate(N_scaf = sum(gap_length, na.rm = T)) %>%
    ungroup()
  return(agp_combo)
}
# Import FASTA index file (.fasta.fai)
readFai <- function(idx_file, suff = NULL) {
  checkPath(idx_file)
  idx <- read.table(idx_file)
  idx <- idx[, 1:2]
  colnames(idx) <- c("seqid", "length")
  idx <- idx %>%
    arrange(desc(length)) %>%
    mutate(seqid = factor(seqid, levels = seqid))
  if (!missing(suff)) {
    colnames(idx) <- paste0(colnames(idx), suff)
  }
  return(idx)
}
# Generate annotated index of pre- and post-Ragout scaffolds
genIdx <- function(ragout_dir, agp_list, seqs) {
  agp <- agp_list[[ragout_dir]] %>% filter(label == "input")
  # Pre-Ragout scaffold FASTA
  pre_genome_file <-
    pull(filter(seqs, grepl(spc_int, Species)), Assembly)
  # Pre-Ragout scaffold index (later filtered for plotting)
  pre_idx_file <- paste0(pre_genome_file, ".fai")
  pre_idx <- readFai(pre_idx_file)
  # Ragout rescaffolded scaffold index
  scaf_file <- list.files(path = ragout_dir,
                          pattern = "_scaffolds\\.fasta\\.fai$",
                          full.names = T)
  scaf <- readFai(scaf_file)
  # Unplaced scaffold index
  unplc_file <- list.files(path = ragout_dir,
                           pattern = "_unplaced\\.fasta\\.fai$",
                           full.names = T)
  unplc <- readFai(unplc_file)
  # Create vector of chimeric scaffolds
  chims <-
    unique(gsub("\\[.*\\]", "", grep("\\[", unplc$seqid, value = T)))
  # Create vector of rescaffolded and unplaced seqids
  rag_seqids <- c(agp$seqid_comp, gsub("\\[.*\\]", "", unplc$seqid))
  # Label scaffolds with "type" column before combining
  # Type: rescaffolded
  scaf <- scaf %>% mutate(type = "pseudochromosomes")
  # Type: remainder
  unplc <- unplc %>% mutate(type = "remainder") %>%
    # Filter out chimeric scaffolds from remainder
    rowwise() %>%
    filter(!grepl("\\[", seqid)) %>%
    ungroup()
  pre_idx <- pre_idx %>%
    # Type: chimera
    mutate(type = case_when(
      seqid %in% chims ~ "chimera",
      # Type: input
      (seqid %in% rag_seqids &
         !seqid %in% chims) ~ "input",
      # Type: excluded
      .default = "excluded"
    ))
  # Combine all scaffold types
  idx <- rbind(scaf, unplc, pre_idx) %>%
    mutate(
      type = factor(type, levels = types, ordered = T),
      seqid = factor(seqid, levels = unique(as.character(seqid))),
      # Convert scaffold lengths to log10(bp) scale
      length = log10(length)
    ) %>%
    arrange(desc(length))
  return(idx)
}
# Annotate AGP dataframe with component synteny to chr0
annotChr0 <- function(idx_agp) {
  comp_chr0 <- read.table("comp_chr0.txt", header = T, sep = "\t")
  comp_chr0 <- comp_chr0 %>%
    mutate(qName = paste0("scaffold_", qNum))
  dict_chr0 <- comp_chr0$qPercent
  names(dict_chr0) <- comp_chr0$qName
  idx_agp <- idx_agp %>%
    mutate(chr0 = case_when(seqid_comp %in% names(dict_chr0) ~
                              dict_chr0[as.character(seqid_comp)],
                            .default = 0))
  return(idx_agp)
}
# Index AGP dataframe with scaffold/contig lengths
idxAgp <- function(ragout_dir, agp_list, idx_list) {
  # AGP dataframe
  agp <- agp_list[[ragout_dir]]
  # Rescaffolded pseudochromosome index
  idx <-
    idx_list[[ragout_dir]] %>% filter(type == "pseudochromosomes")
  # Genome v1.0 index
  pre_idx <- idx_list[[ragout_dir]] %>%
    filter(type %in% c("input", "chimera")) %>%
    arrange(desc(length)) %>%
    mutate(seqid = factor(seqid, levels = unique(as.character(seqid))))
  # Create dictionary of lengths and seqids from each index
  idx_dict <- idx$length
  names(idx_dict) <- idx$seqid
  pre_idx_dict <- pre_idx$length
  names(pre_idx_dict) <- pre_idx$seqid
  idx_agp <- agp %>%
    mutate(
      # Order component seqid factors by length (for plotting)
      seqid_comp = factor(seqid_comp, levels = levels(pre_idx$seqid)),
      # Create length columns
      length_scaf = as.numeric(idx_dict[seqid_scaf]),
      length_comp = as.numeric(pre_idx_dict[seqid_comp]),
      # Convert length out of log10(bp) scale to Mb scale
      length_scaf = (10 ^ length_scaf) * 1e-6,
      length_comp = (10 ^ length_comp) * 1e-6
    )
  idx_agp <- annotChr0(idx_agp)
  idx_agp <- idx_agp %>%
    arrange(desc(length_scaf), start_scaf)
  # Factor seqids by scaffold length and start position (for plotting)
  # order_scaf <- unique(pull(idx_agp, seqid_scaf))
  # order_comp <- fixChrom(levels(pull(idx_agp, seqid_comp)))
  # Arrange data by scaffold length, then start position
  idx_agp <- idx_agp %>%
    mutate(
      seqid_comp_num = fixChrom(seqid_comp),
      seqid_comp_num = factor(seqid_comp_num,
                              levels = unique(seqid_comp_num)),
      seqid_scaf_num = as.factor(as.numeric(factor(
        seqid_scaf,
        levels =
          unique(seqid_scaf)
      ))),
      perc_N = (N_scaf / length_scaf) * 100
    )
  return(idx_agp)
}
# Takes factor vector and returns factor-named vector of colors (i.e., dict)
getColors <- function(vec, palFun, pal = NULL) {
  n_colors <- 3
  if (length(levels(vec)) > n_colors)
    n_colors <- length(levels(vec))
  if (!is.factor(vec))
    stop("Error: vector class is not factor.")
  if (identical(palFun, brewer.pal)) {
    if (missing(pal))
      pal = "Paired"
    vec_colors <- palFun(n_colors, pal)
    names(vec_colors) <- levels(vec)
    # Color chimera red using "Paired" palette
    if (any(names(vec_colors) %in% c("chimera"))) {
      new_ord <-
        c("input",
          "pseudochromosomes",
          "remainder",
          "excluded",
          "chimera")
      if (any(names(vec_colors) %in% c("total")))
        new_ord <- c(new_ord, "total")
      names(vec_colors) <- new_ord
    }
  } else if (identical(palFun, viridis)) {
    if (missing(pal))
      pal = "turbo"
    vec_colors <- palFun(n = n_colors, option = pal)
    names(vec_colors) <- levels(vec)
  }
  return(vec_colors)
}
# Get common legend for variable in list of dataframes using a plotting function
getCommon <- function(plotFun, var_name, df_list, ...) {
  var_count <- sapply(df_list,
                      function(x)
                        length(unique(pull(
                          x, matches(var_name)
                        ))),
                      simplify = F)
  legend_name <- names(which.max(var_count))
  common_legend <- get_legend(plotFun(legend_name, df_list, ...))
  return(common_legend)
}
# Barplot of scaffolds indexed by length
idxPlot <-
  function(ragout_dir,
           idx_list,
           leg_pos = "right",
           pal = "Paired",
           exclude = NULL) {
    x_max <- max(unlist(sapply(idx_list,
                               function(x)
                                 length(unique(pull(x, seqid))),
                               simplify = F)), na.rm = T)
    y_max <-
      max(unlist(sapply(idx_list, pull, length, simplify = F)), na.rm = T)
    common_lims <- coord_cartesian(
      xlim = c(0, x_max + 1),
      ylim = c(0, y_max + 1),
      expand = F
    )
    idx <- idx_list[[ragout_dir]]
    if (!missing(exclude))
      idx <- idx %>% filter(!type %in% exclude)
    # Filter out input scaffolds that are chimera, remainder, or excluded
    if (any(idx$type %in% c("input"))) {
      filt_seqids <- idx %>% filter(type != "input") %>% pull(seqid) %>%
        gsub("\\[.*\\]", "", .)
      idx <-
        idx %>% filter(!(type == "input" & seqid %in% filt_seqids))
    }
    # Colors for plot
    vec_colors <- getColors(idx$type, brewer.pal, pal)
    idx_sum <- idx %>%
      group_by(type) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      # Add total n row
      mutate(type = factor(type, levels = c(types, "total"), ordered = T)) %>%
      rbind(list(factor(
        "total",
        levels = c(types, "total"),
        ordered = T
      ),
      sum(.$n)))
    # Colors for table
    tbl_colors <- getColors(idx$type, brewer.pal, pal)
    tbl_colors[["total"]] <- "white"
    fnt_colors <- rep("black", length(tbl_colors))
    names(fnt_colors) <- names(tbl_colors)
    fnt_colors[names(fnt_colors) %in% c("pseudochromosomes", "excluded")] <-
      "white"
    tbl_colors <- tbl_colors[as.character(unique(idx_sum$type))]
    fnt_colors <- fnt_colors[as.character(unique(idx_sum$type))]
    annot <- annotate(
      geom = "table",
      x = x_max,
      y = y_max + 0.75,
      label = idx_sum,
      table.colnames = F,
      table.theme =
        ttheme_default(
          padding = unit(c(2, 2), "mm"),
          base_size = 12,
          core = list(
            bg_params =
              list(fill = tbl_colors, col = "black"),
            fg_params =
              list(col = fnt_colors)
          )
        )
    )
    idx <- idx %>% arrange(type)
    p <- ggplot(data = idx, mapping = aes(fill = type)) +
      geom_col(mapping = aes(x = as.numeric(row.names(idx)), y = length),
               width = 1) +
      common_lims +
      scale_fill_manual(values = vec_colors) +
      labs(x = "Scaffold index",
           y = "Scaffold length (log10 bp)",) +
      theme_bw() +
      # theme_minimal() +
      annot +
      theme(
        legend.position = leg_pos,
        legend.title = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.background = element_blank()
      )
    return(p)
  }
# Combine pre- and post-Ragout index plots
combPlot <-
  function(ragout_dir,
           idx_list,
           leg_pos = "right",
           ttls = NULL,
           pal = "Paired",
           legend = F) {
    p1 <- idxPlot(
      ragout_dir,
      idx_list = idx_list,
      leg_pos = leg_pos,
      pal = pal,
      exclude = c("pseudochromosomes")
    ) +
      labs(subtitle = "Before",
           x = NULL,
           y = NULL)
    p2 <- idxPlot(
      ragout_dir,
      idx_list = idx_list,
      leg_pos = leg_pos,
      pal = pal,
      exclude = c("input")
    ) +
      labs(subtitle = "After",
           x = NULL,
           y = NULL)
    p_list <- list(p1, p2)
    if (legend) {
      common_legend <-
        getCommon(idxPlot, "type", idx_list, leg_pos = leg_pos)
      p <- ggarrange(
        plotlist = p_list,
        nrow = length(p_list),
        align = "hv",
        legend.grob = common_legend,
        legend = leg_pos
      )
      # p <- annotate_figure(p,
      #                      bottom = "Scaffold index",
      #                      left = "Scaffold length (log10 bp)")
    } else {
      p <- ggarrange(
        plotlist = p_list,
        nrow = length(p_list),
        align = "hv",
        legend = "none"
      )
    }
    if (!is.null(ttls)) {
      ttl <- ttls[[ragout_dir]]
      p <-
        annotate_figure(p, top = text_grob(str_wrap(ttl, width = 35)))
    }
    mar_cm <- 0.5
    p <- annotate_figure(
      p,
      bottom = text_grob("Scaffold index", size = 12),
      left = text_grob("Scaffold length (log10 bp)", size = 12, rot = 90)
    ) +
      theme(plot.margin = margin(mar_cm, mar_cm, mar_cm, mar_cm, "cm"))
    return(p)
  }
# Add ggpubr border to plot
addBorder <- function(p) {
  p <- p + border(size = 0.2)
  return(p)
}
# Function to summarize indexed AGP dataframe
sumAgp <- function(idx_agp) {
  sum_df <- idx_agp %>%
    select(seqid_scaf_num, N_scaf, length_scaf, perc_N) %>% unique()
  tot_len <- sum_df %>% pull("length_scaf") %>% sum()
  tot_N <- sum_df %>% pull("N_scaf") %>% sum()
  tot_perc_N <- (tot_N / tot_len) * 100
  sum_df <- sum_df %>%
    mutate(perc_N = paste0(round(perc_N, digits = 2), "%"))
  return(list(
    sum_df = sum_df,
    tot_len = tot_len,
    tot_N = tot_N,
    tot_perc_N = tot_perc_N
  ))
}
# Plot of ragout rescaffolding
ragoutPlot <-
  function(ragout_dir,
           idx_agp_list,
           leg_pos = "right",
           ttls = NULL,
           pal = "Paired",
           labels = T) {
    ttl <- ttls[[ragout_dir]]
    idx_agp <- idx_agp_list[[ragout_dir]]
    # Common x-axis limits for aligning all plots
    x_max <-
      max(unlist(sapply(idx_agp_list, pull, length_scaf)), na.rm = T) * 1.1
    common_x_lims <- xlim(c(0, x_max + 1))
    # Positions
    xpos <- x_max * 0.8
    ypos <- length(levels(idx_agp$seqid_scaf_num))
    # Summary of AGP dataframe
    sum_list <- sumAgp(idx_agp)
    sum_df <- sum_list$sum_df
    # Color by labels with consistent colors in all graphs
    myColors <- getColors(idx_agp$label, brewer.pal, pal)
    mar_cm <- 0.5
    p <- ggplot(data = idx_agp) +
      geom_segment(
        mapping = aes(
          x = start_scaf,
          xend = end_scaf,
          y = seqid_scaf_num,
          yend = seqid_scaf_num,
          # col = perc_N,
          col = label,
          linewidth = label
        )
      ) +
      annotate(
        geom = "text",
        x = sum_df$length_scaf,
        y = sum_df$seqid_scaf_num,
        label = sum_df$perc_N,
        hjust = -0.3
      ) +
      annotate(
        geom = "text",
        x = xpos,
        y = levels(idx_agp$seqid_scaf_num)[ypos - 2],
        label = paste(paste(
          round(sum_list$tot_len, digits = 2), "Mb"
        ),
        paste0(
          round(sum_list$tot_perc_N, digits = 2), "% N's"
        ),
        sep = "\n")
      ) +
      # annotate(
      #   geom = "text",
      #   x = xpos,
      #   y = levels(idx_agp$seqid_scaf_num)[ypos - 1],
      #   label = paste0(round(sum_list$tot_perc_N, digits = 2), "% N's")
      # ) +
      labs(# title = str_wrap(ttl, width = 35),
        # subtitle = "Scaffold mapping onto pseudochromosomes",
        x = "Scaffold length (Mb)",
        y = "Pseudochromosome index") +
      theme_classic() +
      scale_color_manual(values = myColors) +
      scale_linewidth_manual(values = c(3, 1)) +
      theme(
        legend.position = leg_pos,
        legend.text = element_text(size = rel(1.3)),
        legend.title = element_blank(),
        plot.margin = margin(mar_cm, mar_cm, mar_cm, mar_cm, "cm"),
        axis.title = element_text(size = rel(1.2))
      ) +
      common_x_lims
    # Remove labels from graph if not desired (labels = F)
    if (!labels)
      p <- p + labs(x = "", y = "")
    return(p)
  }
# Run all functions on data
runAnalysis <-
  function(in_dirs,
           seqs,
           pre_gaps,
           plot1,
           plot2,
           leg_pos = "right") {
    # Extract titles for plots
    ttls <- sapply(in_dirs, extTtl)
    lab_lets <- LETTERS[1:length(ttls)]
    print(paste(lab_lets, paste0(ttls, " (", in_dirs, ")"), sep = ": "))
    # Wrangle data
    # Import AGP
    agp_list <- sapply(in_dirs, readAgp, pre_gaps, simplify = F)
    # Import FASTA indices
    idx_list <-
      sapply(in_dirs, genIdx, agp_list, seqs, simplify = F)
    # Index AGP
    idx_agp_list <-
      sapply(in_dirs, idxAgp, agp_list, idx_list, simplify = F)
    # Plots
    # Bar plot length distribution before and after rescaffolding
    p_list <- sapply(
      in_dirs,
      combPlot,
      idx_list = idx_list,
      leg_pos = leg_pos,
      simplify = F
    )
    # # Add borders if more than one plot
    # if (length(names(p_list)) > 1) {
    #   p_list <- sapply(p_list, addBorder, simplify = F)
    # }
    common_legend <-
      getCommon(idxPlot, "type", idx_list, leg_pos = leg_pos)
    all_bar <-
      ggarrange(plotlist = p_list,
                align = "hv",
                labels = "AUTO")
    # Save plot
    ht <- 7
    wd <- 7
    if (length(p_list) > 1) {
      wd <- wd * length(p_list) * 0.3
      ht <- ht * length(p_list) * 0.25
    }
    sapply(
      plot1,
      ggsave,
      plot = all_bar,
      height = ht,
      width = wd,
      bg = "white",
      simplify = F
    )
    print(paste("Saved plot:", plot1))
    
    # Line graph mapping of original scaffolds onto pseudochromosomes
    dot_list <-
      sapply(
        in_dirs,
        ragoutPlot,
        idx_agp_list = idx_agp_list,
        leg_pos = leg_pos,
        # labels = F,
        ttls = ttls,
        simplify = F
      )
    common_legend <-
      getCommon(ragoutPlot, "label", idx_agp_list, leg_pos = leg_pos)
    if (length(dot_list) > 1) {
      all_dot <- ggarrange(
        plotlist = dot_list,
        legend.grob = common_legend,
        legend = leg_pos,
        align = "hv",
        labels = "AUTO"
      )
    } else {
      all_dot = dot_list[[1]]
    }
    # Save plot2
    ht <- 6
    wd <- 7
    if (length(dot_list) > 1) {
      wd <- wd * length(dot_list) * 0.3
      ht <- ht * length(dot_list) * 0.3
    }
    sapply(
      plot2,
      ggsave,
      plot = all_dot,
      height = ht,
      width = wd,
      bg = "white",
      simplify = F
    )
    print(paste("Saved plot:", plot2))
    # # Combine all plots into figure
    # comp_rag <- ggarrange(all_bar, all_dot, align = "hv", nrow = 2)
    # comp_rag <- annotate_figure(comp_rag,
    #                             top =
    #                               text_grob("Reference-based scaffold ordering",
    #                                         face = "bold", size = 14))
    return(list(
      idx_list = idx_list,
      agp_list = agp_list,
      idx_agp_list = idx_agp_list
    ))
  }

# Analysis
# Import data
seqs <- readSpecies(seqFile)
# Tabulate original genome scaffold gaps
pre_genome_file <-
  pull(filter(seqs, grepl(spc_int, Species)), Assembly)
pre_gaps <- findGaps(pre_genome_file, 1e4)
ragout_dirs <- list.files(pattern = "ragout-out-")
# ragout_dirs <- grep("refine|filt", ragout_dirs, value = T)
# ragout_dirs <- grep("solid.*refine", ragout_dirs, value = T)
ragout_dirs <- grep("refine", ragout_dirs, value = T)
result_list1 <-
  runAnalysis(ragout_dir,
              seqs,
              pre_gaps,
              outfiles$len_plot,
              outfiles$map_plot,
              leg_pos = "right")
result_list2 <-
  runAnalysis(
    ragout_dirs,
    seqs,
    pre_gaps,
    outfiles$comp_len_plot,
    outfiles$comp_map_plot,
    leg_pos = "right"
  )
# # All species
# named_fas <- pull(seqs, Assembly, Species)
# ps <-
#   sapply(
#     c(1000, 900, 800, 700, 600, 500, 400, 300, 200, 100, 50),
#     findGapSize,
#     named_fas[c("Saccharina_japonica", "Ectocarpus_sp._Ec32")],
#     simplify = F
#   )
# ps <- ggarrange(plotlist = ps, common.legend = T)
