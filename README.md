# Genome assembly of the brown macroalga *Saccharina latissima* (North American sugar kelp)
Scripts to score and compare the genome assembly of *S. latissima* to related species.

## 1. Input files
### Fetch input files from portals
Fetch assemblies and annotations from [JGI PhycoCosm](https://phycocosm.jgi.doe.gov/phycocosm/home) and [ORCAE](https://bioinformatics.psb.ugent.be/orcae) given a list of JGI portal names and ORCAE links.
##### Usage
###### Option 1: Give JGI username and password
> sbatch fetch_assemblies.sbatch [\<portal_list\>](portal_names.txt) [username] [password]
###### Option 2: Give pre-generated curl login file for JGI
> sbatch fetch_assemblies.sbatch [\<portal_list\>](portal_names.txt) [curl_login_file]
##### Example
```
sbatch s-latissima-genome/fetch_assemblies.sbatch s-latissima-genome/portal_names.txt jgi_login
```
Upon successful download, the output directory `assemblies/` and [\<assembly_file\>](species_table.txt) will be created.

Assembly file format:

Species name  | Genome FASTA | Annotation file | Protein FASTA | Gene info file | ...
------------- | ------------ | --------------- | ------------- | -------------- | ---
Species 1     | PATH1        | PATH1           | PATH1         | PATH1          | ...
Species 2     | PATH2        | PATH2           | PATH2         | PATH2          | ...

### Split artificial chromosomes in assemblies
Run [split.scaffolds.to.contigs](https://github.com/MadsAlbertsen/miscperlscripts/blob/master/split.scaffolds.to.contigs.pl) on artificial chromosomes in *S. japonica* and *Ectocarpus* sp. assemblies.
##### Usage
> sbatch split_artificial.sbatch [\<assembly_file\>](species_table.txt)

##### Example
```
sbatch s-latissima-genome/split_artificial.sbatch s-latissima-genome/species_table.txt
```
The script will also update [\<assembly_file\>](species_table.txt) "Assembly" column with the path to the split assembly FASTA for each species.

## 2. Gene content and assembly statistics
#### Run BUSCO and QUAST on each assembly listed in `<assembly_file>`
##### Usage
> bash genome_stats.sh [\<assembly_file\>](species_table.txt) [[path/to/aug_busco.sbatch]](aug_busco.sbatch) [[path/to/quast.sbatch]](quast.sbatch)
##### Example
```
bash s-latissima-genome/genome_stats.sh s-latissima-genome/species_table.txt s-latissima-genome/aug_busco.sbatch s-latissima-genome/quast.sbatch
```
#### Visualize BUSCO results
##### Usage
> bash genome_stats.sh [\<assembly_file\>](species_table.txt) [[path/to/busco_compare.sbatch]](busco_compare.sbatch)
##### Example
```
bash s-latissima-genome/genome_stats.sh s-latissima-genome/species_table.txt s-latissima-genome/busco_compare.sbatch
```
For each lineage, a plot of BUSCO results across all genomes will be saved to `busco_<lineage>.png`.
![alt text](busco_eukaryota_odb10.png)
![alt text](busco_stramenopiles_odb10.png)

## 3. Evaluation of genome contiguity (scaffoldedness)
Generates violin plots of contig size for each genome.
##### Usage
> sbatch scaffold_eval.sbatch [\<assembly_file\>](species_table.txt) \<species_of_interest\> \<outgroup_species\>

>> **Note**: Ensure `<species_of_interest>` and `<outgroup_species>` are in the format "Genus_specificname" and do not contain spaces.
##### Example
```
sbatch s-latissima-genome/scaffold_eval.sbatch s-latissima-genome/species_table.txt Saccharina_latissima Ectocarpus_sp.
```
Violin plot of scaffold and contig sizes per genome will be saved in `scaffold_sizes_violin.png`.
![alt text](scaffold_sizes_violin.png)

## 4. Multi-species whole genome alignment with [Progressive Cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md)
### Prune brown macroalgae phylogeny to include only species in analysis
Default [tree](https://ars.els-cdn.com/content/image/1-s2.0-S1055790319300892-mmc1.txt) sourced from [Starko, S. et al. 2019](https://doi.org/10.1016/j.ympev.2019.04.012).
##### Usage
> sbatch [prune_tree.sbatch](prune_tree.sbatch) [\<assembly_file\>](species_table.txt) [[tree]](https://ars.els-cdn.com/content/image/1-s2.0-S1055790319300892-mmc1.txt)
##### Example
```
sbatch s-latissima-genome/prune_tree.sbatch s-latissima-genome/species_table.txt 1-s2.0-S1055790319300892-mmc1.txt
```
Output `<seqFile>` formatted for Cactus: `s_lat_alignment.txt`.
Phylogeny before and after pruning will be plotted in `phylo_prune.png`.
![Phylogenetic graphs showing pruning](phylo_prune.png)

### Run Progressive Cactus alignment
#### Prepare scripts for stepwise pipeline
##### Usage
> sbatch [cactus_prepare.sbatch](cactus_prepare.sbatch) [\<seqFile\>](s_lat_alignment.txt)
##### Example
```
sbatch s-latissima-genome/cactus_prepare.sbatch s-latissima-genome/s_lat_alignment.txt
```
#### Run scripts sequentially
> sbatch [cactus_run_prepared.sbatch](cactus_run_prepared.sbatch) [sbatch_list_file](cactus_sbatch_list.txt)
##### Example
```
sbatch s-latissima-genome/cactus_run_prepared.sbatch s-latissima-genome/cactus_sbatch_list.txt
```
### Visualize alignments and calculate statistics
#### Run [halSynteny](https://github.com/ComparativeGenomicsToolkit/hal) to extract syntenic blocks between each genome pair
##### Usage
> sbatch [halSynteny.sbatch](halSynteny.sbatch) \<inHal\>
##### Example
```
sbatch s-latissima-genome/halSynteny.sbatch cactus-steps-output/s_lat_alignment.hal
```
#### Generate report and alignment graphs of pairwise halSynteny PSL files
Description of PSL file format taken from [Ensembl](https://useast.ensembl.org/info/website/upload/psl.html).
##### Usage
> sbatch [align_analysis.sbatch](align_analysis.sbatch) [\<seqFile\>](s_lat_alignment.txt)
##### Example
```
sbatch s-latissima-genome/align_analysis.sbatch s-latissima-genome/s_lat_alignment.txt
```
#### Summary statistics
Alignment report: alignment_report.tsv

<!-- align_table_start -->

|Statistic                                                 |Ectocarpus_siliculosus |Undaria_pinnatifida |Macrocystis_pyrifera |Saccharina_japonica |All    |
|:---------------------------------------------------------|:----------------------|:-------------------|:--------------------|:-------------------|:------|
|Mean n homologous scaffolds mapped per chromosome         |31.57                  |31.03               |18.18                |35.09               |26.94  |
|SD n homologous scaffolds mapped per chromosome           |11.38                  |16.71               |15.5                 |60.59               |31.23  |
|Median n homologous scaffolds mapped per chromosome       |31.5                   |36                  |16                   |25.5                |24     |
|Total n homologous scaffolds mapped                       |884                    |1024                |1091                 |1123                |4122   |
|Mean length of homologous scaffolds (Mb) per chromosome   |19                     |17.77               |10.11                |19.17               |15.28  |
|SD length of homologous scaffolds (Mb) per chromosome     |6.95                   |8.62                |8.34                 |20.77               |12.52  |
|Median length of homologous scaffolds (Mb) per chromosome |17.46                  |17.23               |9.37                 |15.23               |15.23  |
|Total length of homologous scaffolds (Mb)                 |531.98                 |586.3               |606.34               |613.38              |2338   |
|Mean exact matches (Mb) per chromosome                    |0.62                   |1.07                |1.52                 |5.68                |2.13   |
|SD exact matches (Mb) per chromosome                      |0.38                   |0.74                |1.39                 |5.7                 |3.31   |
|Mean exact matches (%) per chromosome                     |9.32                   |6.76                |18.92                |41.14               |19.19  |
|Median exact matches (Mb) per chromosome                  |0.52                   |0.95                |1.12                 |4.68                |1.19   |
|Total exact matches (Mb)                                  |17.24                  |35.36               |90.93                |181.89              |325.43 |


<!-- align_table_end -->

![alt text](align_length_and_n.png)

#### Cluster query scaffolds by repeated alignments to same reference chromosome
![alt text](elbow_k.png)
![alt text](clust_complexheatmap.png)

#### Heatmap colored by alignment coverage
![alt text](align_clustered.png)

## 5. Homology-based rescaffolding
### Extract MAF from Cactus alignment using [HAL tools](https://github.com/ComparativeGenomicsToolkit/hal)
##### Usage
> sbatch [hal2maf.sbatch](hal2maf.sbatch) \<halFile\>
##### Example
```
sbatch s-latissima-genome/hal2maf.sbatch cactus-steps-output/s_lat_alignment.hal
```
### Rescaffold assembly v1.0 with [Ragout](https://github.com/fenderglass/Ragout)
##### Usage
> sbatch [ragout.sbatch](ragout.sbatch) [\<seqFile\>](s_lat_alignment.txt) \<cactus_dir\> [ragout_prefix]
##### Example
```
sbatch s-latissima-genome/ragout.sbatch s-latissima-genome/s_lat_alignment.txt cactus-steps-output ragout-out
```
![alt text](ragout_length_dist.png)
![alt text](ragout_pseudo_mapping.png)
![alt text](ragout_comp_length_dist.png)
![alt text](ragout_comp_pseudo_mapping.png)

## 6. Identify orthologous genes with OrthoFinder
### Run OrthoFinder on all brown algal species
##### Usage
> sbatch [orthofinder.sbatch](orthofinder.sbatch) [\<assembly_file\>](species_table.txt)
##### Example
```
sbatch s-latissima-genome/orthofinder.sbatch s-latissima-genome/species_table.txt
```





















