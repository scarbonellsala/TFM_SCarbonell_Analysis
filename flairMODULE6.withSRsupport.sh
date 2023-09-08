#!/bin/bash

set -xe
module load Python/3.8.2-GCCcore-9.3.0
module load BEDTools/2.29.2-GCC-9.3.0
module load SAMtools/1.11-GCC-9.3.0

#flair 12346 -r reads.fa -g genome.fa -f annotation.gtf -o flair.output --temp_dir temp_flair [optional arguments]
#(module numbers: align=1, correct=2, collapse=3, collapse-range=3.5, quantify=4, diffExp=5, diffSplice=6)


isoforms_bed="[path_to_flair.collapse.isoforms.bed]"
counts_matrix="[path_to_flair.quantify.counts.tsv]"

main_dir="[path_to_outputdir]"
data_test="$main_dir/data_test"

mkdir -p $data_test
cd $data_test/

flair diffSplice -i $isoforms_bed -q $counts_matrix --threads 4

