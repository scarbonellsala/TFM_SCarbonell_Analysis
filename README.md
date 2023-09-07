# TFM_SCarbonell_Analysis

## Workflow:

![Workflow](workflow.png)

## To process sequencing long-read data:
- **ONT Basecall and QC**:
Guppy v6 SUP was used to basecall the data using default parameters.
Default parameters had been used to run Nanoplot and MultiQC.

- **Long Read Splicing analysis**
FLAIR v2.0 modules had been used to process long-read sequencing data in 2 different modes: LR-only and SR-supported. These program modules were run as a job in the cluster.

TODO: EXAMPLE script here

- **Short-read splicing junctions (SJs)**:
intronPorspector was used to generate the SJs file from short-read data, running the following bash script:

```
run.intronProspector.sh
```
  
- **Data QC and filtering**:

SQANTI3 was used to QC the long-read data models obtained by FLAIR. 

TODO: EXAMPLE script here

SQANTI categories were used to filter for Known (Gencode v43 annotated Genes), running the following bash script:

```
run.sqanti.filter.known.sh
```

## To process sequencing short-read data:

GRAPE-NF 

```
module load Java/11

nextflow -bg run grape-nf -r dev --rsemSkipCi --index /nfs/users/rg/scarbonell/TFM/short_reads/metadata.tsv --genome /nfs/users/rg/projects/references/Genome/H.sapiens/GRCh38/GRCh38.p13.primary_assembly.genome.fa.gz --annotation /nfs/users/rg/projects/references/Annotation/H.sapiens/gencode43/gencode.v43.primary_assembly.annotation.gtf.gz --rg-platform ILLUMINA --rg-center-name CRG -resume -c /software/rg/grape/config/rg.singularity.dsl2.config > pipeline.log
![image](https://github.com/scarbonellsala/TFM_SCarbonell_Analysis/assets/75372182/7b347584-4745-4702-a555-d11b9258ffcf)

```

## R-scripts:

- **PCA_scaled_centered.R**: Build scaled and centred PCAs for LRonly and SRsupported data.
- **Rscript_different_isoform_usage.R**: Perform different isoform usage analysis.
- **Splicing_events.R**: Analyze splicing events (IR, ES, Alt3, and Alt5).
- **Upset_plots_addingGENCODE.R**: Prepare upset plots of isoform counts and GENCODE for LRonly and SRsupported data.
- **PLOT_regression_plots_replicates.R**: Generate regression plots for LRonly and SRsupported data.

## Bash-scripts:

- **run.ggsashimi.sh**: Generate ggsashimi plots.
- **run.sqanti.filter.known.sh**: Filter known genes using SQANTI categories.
- **run.intronProspector.sh**: from BAM/SAM generates a list of splicing junctions (SJs)
