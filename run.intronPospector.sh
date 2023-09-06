names="H0TRep2_m4_n10_toGenome H0TRep1_m4_n10_toGenome H0NRep2_m4_n10_toGenome H0NRep1_m4_n10_toGenome                  
  H0CRep2_m4_n10_toGenome H0CRep1_m4_n10_toGenome"                                                                        
  
  for name in $names ; do
    intronProspector --genome-fasta=../refs/GRCh38.primary_assembly.genome.fa --skip-missing-targets
    --min-confidence-score=1.0 --intron-calls=$name.introns.tsv $name.bam &
  done
  wait
  intronProspectorMerge --intron-bed=H0.introns.bed9 H0*.introns.tsv
  cut -f 1-6 H0.introns.bed9  > H0.introns.bed
