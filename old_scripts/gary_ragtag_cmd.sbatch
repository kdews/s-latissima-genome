minimap2 -x map-pb -a -Q /scratch/gmolano/chloroplast_assembly_gary/MZ156032.1_Macrocystis_pyrifera_plastid_complete_genome.fasta /project/noujdine_61/jdiesel/CANU_results_9.1/A01_clean_Kelp_canu1.9/asm.correctedReads.fasta.gz > CI_03_corrected_reads_on_MZ156032.1_Macrocystis_pyrifera_plastid_complete_genome.sam

canu -assemble -p chloroplast -d CI_03_chloro genomeSize=130000 minMemory=48 merylMemory=48 useGridutgovl=True gridOptionsutgovl="--time=24:00:00 --partition=cegs --ntasks=12" gridOptions="--time=24:00:00 --partition=cegs" gridOptionsJobName="CI_03_cc"  -pacbio-corrected CI_03_corrected_reads_mapped_to_putative_CI_03_chloro.fa >canu_CI_03_chloro.log 2>&1

#>tig00000002 len=178698 reads=483 class=contig suggestRepeat=yes suggestBubble=no suggestCircular=yes trim=17088-147301
#>tig00000003 len=61292 reads=5 class=contig suggestRepeat=yes suggestBubble=no suggestCircular=no trim=0-61292
#>tig00000004 len=136743 reads=137 class=contig suggestRepeat=yes suggestBubble=no suggestCircular=yes trim=4654-134859