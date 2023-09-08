#!/bin/bash

set -x
conda activate SQANTI3.env
export PYTHONPATH=/users/rg/scarbonell/bin/SQANTI3-5.1.2/cDNA_Cupcake/sequence/:/users/rg/scarbonell/bin/SQANTI3-5.1.2/cDNA_Cupcake/


genome="[path_to_GRCh38.p13.primary_assembly.genome.fa.gz]"
annotation="[path_to__gencode.v43.primary_assembly.annotation.gtf]"
soforms_gtf="[path_to_flair.collapse.isoforms.gtf]"

main_dir="[path_to_outputfiles]"
data_SRsupport="$main_dir/data_SRsupport"

mkdir -p $data_SRsupport
cd $data_SRsupport/

python /users/rg/scarbonell/bin/SQANTI3-5.1.2/sqanti3_qc.py $isoforms_gtf $annotation $genome --force_id_ignore
