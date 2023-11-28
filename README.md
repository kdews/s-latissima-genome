# Genome assembly of the brown macroalga *Saccharina latissima* (North American sugar kelp)
Scripts to score and compare the genome assembly of *S. latissima* to related species.

## 1. Input files
Fetch assemblies and annotations from JGI website given a list of JGI portal names.
```
# Give JGI username and password
sbatch fetch_assemblies.sbatch <portal_list> [username] [password]
# Give pre-generated curl login file for JGI
sbatch fetch_assemblies.sbatch <portal_list> [curl_login_file]
```
Upon successful download, output directory ```assemblies/``` and ```assembly_file``` (default: ```assemblies/species_table.txt```) will be created.
Assembly file format:
```
Species_name	Genome_PATH	Annotation_PATH
```

## 2. Genome scoring with BUSCO and QUAST
Run BUSCO and QUAST on each assembly listed in ```assembly_file```:
```
bash genome_stats.sh <assembly_file> aug_busco.sbatch quast.sbatch
```

Visualize BUSCO results.
```
bash genome_stats.sh <assembly_file> busco_compare.sbatch
```

## 3. Multi-species whole genome alignment with Cactus Progressive Aligner
### Prune brown macroalgae phylogeny to include only species in analysis
Default [tree](https://ars.els-cdn.com/content/image/1-s2.0-S1055790319300892-mmc1.txt) sourced from [Starko, S. et al. 2019](https://doi.org/10.1016/j.ympev.2019.04.012).
```
sbatch prune_tree.sbatch <assembly_file> [tree]
```
Output will be in ```<tree>_pruned``` and ```s_lat_alignment.txt```.

### Run Cactus aligner

