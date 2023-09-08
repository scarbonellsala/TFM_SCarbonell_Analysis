#!/bin/bash

set -xe
conda activate

# usage: diffsplice_fishers_exact events.quant.tsv colname1 colname2 out.fishers.tsv

# From the path to outputdir

# Define the list of event files
event_files=("flair.diffsplice.alt3.events.quant.tsv" "flair.diffsplice.alt5.events.quant.tsv" "flair.diffsplice.es.events.quant.tsv" "flair.diffsplice.ir.events.quant.tsv")

# Define the list of replicate pairs
replicates=("H0CRep1 H0CRep2" "H0NRep1 H0NRep2" "H0TRep1 H0TRep2")

# Loop through each replicate pair
for replicate_pair in "${replicates[@]}"; do
    replicate_array=($replicate_pair)
    replicate1="${replicate_array[0]}"
    replicate2="${replicate_array[1]}"
    
    # Loop through each event file
    for event_file in "${event_files[@]}"; do
        output_file="out.fishers.${replicate1}.${replicate2}.${event_file%.tsv}.cytosol.tsv"
        
        # Execute diffsplice_fishers_exact command
        ~bin/flair/bin/diffsplice_fishers_exact "$event_file" "$replicate1" "$replicate2" "$output_file"
    done
done



