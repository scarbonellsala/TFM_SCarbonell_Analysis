#!/bin/bash

set -x
conda activate SQANTI3.env
export PYTHONPATH=/users/rg/scarbonell/bin/SQANTI3-5.1.2/cDNA_Cupcake/sequence/:/users/rg/scarbonell/bin/SQANTI3-5.1.2/cDNA_Cupcake/

genome="/nfs/users/rg/projects/references/Genome/H.sapiens/GRCh38/GRCh38.p13.primary_assembly.genome.fa"
annotation="/nfs/users/rg/projects/references/Annotation/H.sapiens/gencode43/gencode.v43.primary_assembly.annotation.gtf"
isoforms_gtf="/nfs/users/project/gencode_006070_no_backup/scarbonell/TFM/long_reads/FLAIR/all/data_SRsupport/flair.collapse.isoforms.gtf"

main_dir="[path_to_outputfiles]"
data_SRsupport="$main_dir/data_SRsupport"

mkdir -p $data_SRsupport
cd $data_SRsupport/

python /users/rg/scarbonell/bin/SQANTI3-5.1.2/sqanti3_qc.py $isoforms_gtf $annotation $genome --force_id_ignore
