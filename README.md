# Genome assembly of the brown macroalga *Saccharina latissima* (North American sugar kelp)
Scripts to score and compare the genome assembly of *S. latissima* to related species.

## 1. Input files
Fetch assemblies and annotations from JGI website and ORCAE given a list of JGI portal names and ORCAE links.
```
# Usage: give JGI username and password, or give pre-generated curl login file for JGI
sbatch fetch_assemblies.sbatch <portal_list> [username] [password]
```
```
sbatch fetch_assemblies.sbatch <portal_list> [curl_login_file]
# Example
sbatch fetch_assemblies.sbatch s-latissima-genome/portal_names.txt jgi_login
```
Upon successful download, output directory (`assemblies/`) and `<assembly_file>` (`species_table.txt`) will be created.
Assembly file format:
```
Species_name	Genome_PATH	Annotation_PATH	Gene_Info_PATH
```

## 2. Evaluation of scaffoldedness and contig size filtering
Generates violin plots of contig size for each genome, then filters out extremely large or small (>1 Mb) contigs.
```
# Usage
sbatch chromosome_extract.sbatch <assembly_file> <species_of_interest>
```
```
# Example
sbatch chromosome_extract.sbatch s-latissima-genome/species_table.txt Saccharina_latissima
```
Note: Ensure `<species_of_interest>` is in the format "Genus_specificname" and does not contain spaces.

Resulting filtered genomes will be tabulated next to species names in `filt_species_table.txt`. Violin plots of contig sizes before and after filtering will be saved in `scaffold_sizes_violin.png`.

![alt text](https://github.com/kdews/s-latissima-genome/blob/main/scaffold_sizes_violin.png)

## 3. Genome scoring with BUSCO and QUAST
Run BUSCO and QUAST on each assembly listed in `<assembly_file>`:
```
# Usage
bash genome_stats.sh <assembly_file> aug_busco.sbatch quast.sbatch
```
```
# Examples
bash genome_stats.sh s-latissima-genome/species_table.txt aug_busco.sbatch quast.sbatch
bash genome_stats.sh s-latissima-genome/filt_species_table.txt aug_busco.sbatch quast.sbatch
```

Visualize BUSCO results.
```
# Usage
bash genome_stats.sh <assembly_file> busco_compare.sbatch
```
```
# Example
bash genome_stats.sh s-latissima-genome/species_table.txt busco_compare.sbatch
```

## 4. Multi-species whole genome alignment with Cactus Progressive Aligner
### Prune brown macroalgae phylogeny to include only species in analysis
Default [tree](https://ars.els-cdn.com/content/image/1-s2.0-S1055790319300892-mmc1.txt) sourced from [Starko, S. et al. 2019](https://doi.org/10.1016/j.ympev.2019.04.012).
```
# Usage
sbatch prune_tree.sbatch <assembly_file> [tree]
```
```
# Example
sbatch prune_tree.sbatch s-latissima-genome/species_table.txt 1-s2.0-S1055790319300892-mmc1.txt
```
Output `<seqFile>` formatted for Cactus: `s_lat_alignment.txt`. Phylogeny before and after pruning will be plotted in `phylo_prune.png`.
![alt text](https://github.com/kdews/s-latissima-genome/blob/main/phylo_prune.png)

### Run Cactus aligner
#### Prepare scripts for stepwise pipeline
```
# Usage
sbatch cactus_prepare.sbatch <seqFile>
```
```
# Example
sbatch cactus_prepare.sbatch s-latissima-genome/s_lat_alignment.txt
```
#### Run scripts sequentially
```
sbatch cactus_run_prepared.sbatch
```

### Visualize alignment

