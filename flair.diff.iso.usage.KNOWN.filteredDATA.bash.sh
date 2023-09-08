#!/bin/bash

set -xe
module load Python/3.8.2-GCCcore-9.3.0
module load BEDTools/2.29.2-GCC-9.3.0
module load SAMtools/1.11-GCC-9.3.0


# usage: diff_iso_usage counts_matrix colname1 colname2 diff_isos.txt
# set the main_dir to be your output directory


# LRonly

main_dir="[path_to_output_dir]"
known_diffIsoformUsage="$main_dir/known_diffIsoformUsage"

mkdir -p $known_diffIsoformUsage
cd $known_diffIsoformUsage/

diff_iso_usage known_sample_counts.tsv H0C H0T diff_isos.cytosolVSwhole.filtered.txt
diff_iso_usage known_sample_counts.tsv H0N H0T diff_isos.nucleusVSwhole.filtered.txt
diff_iso_usage known_sample_counts.tsv H0C H0N diff_isos.cytosolVSnucleus.filtered.txt

# SRsupport

main_dir="/users/project/gencode_006070_no_backup/scarbonell/TFM/long_reads/FLAIR/all//filtering/data_SRsupport/"
known_diffIsoformUsage="$main_dir/known_diffIsoformUsage"

mkdir -p $known_diffIsoformUsage
cd $known_diffIsoformUsage/

diff_iso_usage known_sample_counts.tsv H0C H0T diff_isos.cytosolVSwhole.filtered.txt
diff_iso_usage known_sample_counts.tsv H0N H0T diff_isos.nucleusVSwhole.filtered.txt
diff_iso_usage known_sample_counts.tsv H0C H0N diff_isos.cytosolVSnucleus.filtered.txt
