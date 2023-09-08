#!/bin/bash

module load Python/3.8.2-GCCcore-9.3.0
module load BEDTools/2.29.2-GCC-9.3.0
module load SAMtools/1.11-GCC-9.3.0

#flair align -g genome.fa -r <reads.fq.gz>|<reads.fq>|<reads.fa> [options]

genom="[path_to_GRCh38.p13.primary_assembly.genome.fa.gz]"
minimap_path="[path_to_minimap2]/"
jobarray_path="[path_to_TFM.jobARRAY.tsv]"

IFS=$'\t' read -r sample fastq_path <<< "$(sed -n ${SGE_TASK_ID}p ${jobarray_path})"

main_dir="[path_to_outputdir]"
data="$main_dir/data"

mkdir -p $data/$sample
cd $data/$sample

flair align -g $genom --threads ${NSLOTS} -r $fastq_path -m $minimap_path
