# TFM_SCarbonell_Analysis

## Workflow:

![Workflow](workflow.png)

## To process sequencing long-read data:
- **ONT Basecall and QC**:
Guppy v6 SUP was used to basecall the data using default parameters.
Default parameters had been used to run Nanoplot and MultiQC.

- **Long Read Splicing analysis**
FLAIR v2.0 modules had been used to process long-read sequencing data. These program modules were run as a job in the cluster.

TODO: EXAMPLE script here

- **Short-read splcing junctions (SJs)**:
intronPorspector was used to generate the SJs file from short-read data, running the following bash script:

```
run.intronProspector.sh
```
  
- **Data filtering**:

SQANTI

SQANTI filter for Known

```
run.sqanti.filter.known.sh
```

## To process sequencing long-read data:

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
