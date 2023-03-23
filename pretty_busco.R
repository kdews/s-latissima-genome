wd <- "/scratch2/kdeweese/latissima/genome_stats/"
setwd(wd)
rsource <- "busco_summaries/eukaryota_odb10/busco_figure.R"
source(rsource)

df
spec_names <- read.table("assemblies/pretty_names.txt")
spec_names$V2 <- gsub("\\..*", "", spec_names$V2)
merge(df, spec_names, by.x = "my_species", by.y = "V1")
