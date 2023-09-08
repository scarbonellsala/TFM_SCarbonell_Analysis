# TFM_SCarbonell_Analysis

## Workflow:

![Workflow](workflow.png)

## To process sequencing long-read data:
- **ONT Basecall and QC**:
Guppy v6 SUP was used to basecall the data using default parameters.
Default parameters had been used to run Nanoplot and MultiQC.

- **Long Read Splicing analysis**
FLAIR v2.0 modules had been used to process long-read sequencing data in 2 different modes: LR-only and SR-supported. These program modules were run as a job in the cluster.

Flair align, correct and collapse:

```
flairMODULE123.sh
```

Flair quantify:

```
flairMODULE4.sh
```

- **Short-read splicing junctions (SJs)**:
intronProspector was used to generate the SJs file from short-read data, running the following bash script:

```
run.intronProspector.sh
```
  
- **Data QC and filtering**:

SQANTI3 was used to QC the long-read data models obtained by FLAIR. 

```
sqanti3.withSRsupport.sh
```

SQANTI categories were used to filter for Known (Gencode v43 annotated Genes), running the following bash script:

```
run.sqanti.filter.known.sh
```

## To process sequencing short-read data:

GRAPE-NF has been used to process short-read sequencing data

```
module load Java/11

nextflow -bg run grape-nf -r dev --rsemSkipCi --index /nfs/users/rg/scarbonell/TFM/short_reads/metadata.tsv --genome /nfs/users/rg/projects/references/Genome/H.sapiens/GRCh38/GRCh38.p13.primary_assembly.genome.fa.gz --annotation /nfs/users/rg/projects/references/Annotation/H.sapiens/gencode43/gencode.v43.primary_assembly.annotation.gtf.gz --rg-platform ILLUMINA --rg-center-name CRG -resume -c /software/rg/grape/config/rg.singularity.dsl2.config > pipeline.log
![image](https://github.com/scarbonellsala/TFM_SCarbonell_Analysis/assets/75372182/7b347584-4745-4702-a555-d11b9258ffcf)

```

## Data analysis:

special modules from FLAIR were used to perform splicing events and isoform usage analysis.

```
flairMODULE6.withSRsupport.sh
```

```
flair.diff.iso.usage.KNOWN.filteredDATA.bash.sh
```

## Data visualization:

For long-read sequencing data FLAIR align module was run for each sample independently to obtain separate BAM files for data visualization.

```
flair.align.singlebams.sh
```

While separate BAM files from short-read sequencing data were obtained directly from Grape-nf output. 

Samtools was used to merge sample replicate BAM files from different compartments prior to data visualization.

USCS browser and IGV were used to check splicing events and isoform usage. The exact coordinates to plot were extracted and used to run sushi function on ggsashimi

```
run.ggsashimi.sh
```

## R-scripts:

- **PCA_scaled_centered.R**: Build scaled and centred PCAs for LRonly and SRsupported data.
- **Rscript_different_isoform_usage.R**: Perform different isoform usage analysis.
- **Splicing_events.R**: Analyze splicing events (IR, ES, Alt3, and Alt5).
- **Upset_plots_addingGENCODE.R**: Prepare upset plots of isoform counts and GENCODE for LRonly and SRsupported data.
- **PLOT_regression_plots_replicates.R**: Generate regression plots for LRonly and SRsupported data.

## Bash-scripts:

- **flairMODULE123.sh**: Bash script to run FLAIR along, correct and quantify modules from independent long-read fastq by simultaneously processing all samples and building a unique set of transcript models.
- **flairMODULE4.sh**: Bash script to run FLAIR quantify module which performs isoform quantifications matrix from transcript models file.
- **run.intronProspector.sh**: From BAM/SAM generates a list of splicing junctions (SJs)
- **sqanti3.withSRsupport.sh**: Bash script to run SQANTI3 (example for short-read supported data)
- **run.sqanti.filter.known.sh**: Filter known genes using SQANTI categories.
- **flairMODULE6.withSRsupport.sh**: Bash script to run FLAIR diffSplice - module 6 (example for short-read supported data)
- **flair.diff.iso.usage.KNOWN.filteredDATA.bash.sh**: Bash script to run FLAIR different isoform usage tool (example for known genes filtered data)
- **flair.align.singlebams.sh**: from long-read fastq generates single BAMs for each sample using FLAIR align module 1
- **run.ggsashimi.sh**: Uses sushi function to generate ggsashimi plots.
