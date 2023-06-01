# Set working directory
setwd("/scratch2/kdeweese/latissima/genome_stats/")
# Input files
# busco_tables <- list.files(pattern = "full_table\\.tsv$", recursive = T)
species_table <- "assemblies/pretty_names.txt"



species <- read.table(species_table, sep = "\t")
gsubfunc <- function(x) {
  gsub("\\.fasta", "", x)
}
species[2] <- sapply(species[2], gsubfunc)
lat_busco_tables <- grep("s_lat_genome|Sugarkelp", busco_tables, value = T)
buscos <- list()
for (spc in species[[2]]) {
  fns <- grep(spc, busco_tables, value = T)
  for (fn in fns) {
    print(fn)
    tab <- read.delim(fn, skip = 2)
    tab <- tab[1:3]
    colnames(tab)[2:3] <- colnames(tab)[2:3]
    # buscos <- append(buscos, tab)
  }
  # lin <- gsub(paste(species_ids, collapse = "|"), "", fn)
  # species <- gsub(".*busco_|_euk.*|_stra.*", "", fn)
  # species <- sub("_euk.*|_stra.*", "", species)
  # print(lin)
}

