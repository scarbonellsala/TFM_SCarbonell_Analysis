#!/bin/bash

module load Python/3.8.2-GCCcore-9.3.0
module load BEDTools/2.29.2-GCC-9.3.0
module load SAMtools/1.11-GCC-9.3.0

#flair 12346 -r reads.fa -g genome.fa -f annotation.gtf -o flair.output --temp_dir temp_flair [optional arguments]
#(module numbers: align=1, correct=2, collapse=3, collapse-range=3.5, quantify=4, diffExp=5, diffSplice=6)

minimap_path="[path_to_minimap2]/"
reads_manifest="[path_to_reads_manifest.tsv]"
isoforms="[path_to_flair.collapse.isoforms.fa]"

main_dir="[path_to_outputdir]"
data="$main_dir/data"

mkdir -p $data
cd $data/

flair quantify -r $reads_manifest -i $isoforms --sample_id_only --threads 8 -m $minimap_path
