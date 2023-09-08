# TFM_SCarbonell_Analysis

## Workflow:

![Workflow](workflow.png)

## Processing Long-Read Data:

- **ONT Basecall and QC:**
  - Guppy v6 SUP was used to basecall the data with default parameters.
  - Nanoplot and MultiQC were run using default parameters.

- **Long-Read Splicing Analysis:**
  - FLAIR v2.0 modules were used to process long-read sequencing data in two modes: LR-only and SR-supported.
  - Run Flair align, correct, and collapse:
    ```
    flairMODULE123.sh
    ```
  - Run Flair quantify:
    ```
    flairMODULE4.sh
    ```

- **Short-Read Splicing Junctions (SJs):**
  - intronProspector was used to generate SJs from short-read data using the following bash script:
    ```
    run.intronProspector.sh
    ```

- **Data QC and Filtering:**
  - SQANTI3 was used to QC the long-read data models obtained by FLAIR.
  - Filter for Known (Gencode v43 annotated Genes) using SQANTI categories with the following script:
    ```
    run.sqanti.filter.known.sh
    ```

## Processing short-read data:

Short-read sequencing data was processed using the GRAPE-NF pipeline.

To process short-read sequencing data with GRAPE-NF, use the following command:

```
module load Java/11

nextflow -bg run grape-nf -r dev --rsemSkipCi --index /nfs/users/rg/scarbonell/TFM/short_reads/metadata.tsv --genome /nfs/users/rg/projects/references/Genome/H.sapiens/GRCh38/GRCh38.p13.primary_assembly.genome.fa.gz --annotation /nfs/users/rg/projects/references/Annotation/H.sapiens/gencode43/gencode.v43.primary_assembly.annotation.gtf.gz --rg-platform ILLUMINA --rg-center-name CRG -resume -c /software/rg/grape/config/rg.singularity.dsl2.config > pipeline.log
![image](https://github.com/scarbonellsala/TFM_SCarbonell_Analysis/assets/75372182/7b347584-4745-4702-a555-d11b9258ffcf)
```

## Data Analysis

### Splicing Events and Isoform Usage Analysis

- **FLAIR Modules for Splicing Events and Isoform Usage:**
  - FLAIR modules were employed to analyze splicing events and isoform usage.

    ```
    flairMODULE6.withSRsupport.sh
    ```

    ```
    flair.diff.iso.usage.KNOWN.filteredDATA.bash.sh
    ```

- **Fisher Test for Splicing Events:**
  - To conduct Fisher tests on splicing events:

    ```
    diffsplice.fishers.exact.sh
    ```

### R Scripts for Analysis

In addition to the tools mentioned above, R scripts were utilized for further analysis and to generate visualizations of the data. These R scripts played a crucial role in providing insights and facilitating data interpretation.

Please refer to the specific R script files for detailed information on the analyses conducted and the visualizations produced.

## Splicing Data Visualization

- **Long-Read Data Visualization:**
  - For long-read sequencing data, the FLAIR align module was run independently for each sample to obtain separate BAM files for data visualization.

    ```
    flair.align.singlebams.sh
    ```

- **Short-Read Data Visualization:**
  - Separate BAM files from short-read sequencing data were obtained directly from Grape-nf output.

- **BAM File Merge:**
  - Samtools was used to merge sample replicate BAM files from different compartments prior to data visualization.

- **Splicing Events and Isoform Usage Visualization:**
  - USCS browser and IGV were used to inspect splicing events and isoform usage. The exact coordinates for plotting were extracted and utilized to run the sushi function on ggsashimi.

    ```
    run.ggsashimi.sh
    ```

## R-scripts:

- **PCA_scaled_centered.R**:
  - Purpose: Build scaled and centered PCAs for LRonly and SRsupported data.

- **Rscript_different_isoform_usage.R**:
  - Purpose: Perform different isoform usage analysis.

- **Splicing_events.R**:
  - Purpose: Analyze splicing events (IR, ES, Alt3, and Alt5).

- **Upset_plots_addingGENCODE.R**:
  - Purpose: Prepare upset plots of isoform counts and GENCODE for LRonly and SRsupported data.

- **PLOT_regression_plots_replicates.R**:
  - Purpose: Generate regression plots for LRonly and SRsupported data.

## Bash-scripts:

- **flairMODULE123.sh**:
  - Purpose: Bash script to run FLAIR along, correct, and quantify modules from independent long-read fastq by simultaneously processing all samples and building a unique set of transcript models.

- **flairMODULE4.sh**:
  - Purpose: Bash script to run FLAIR quantify module, which performs isoform quantifications matrix from transcript models file.

- **run.intronProspector.sh**:
  - Purpose: From BAM/SAM generates a list of splicing junctions (SJs).

- **sqanti3.withSRsupport.sh**:
  - Purpose: Bash script to run SQANTI3 (example for short-read supported data).

- **run.sqanti.filter.known.sh**:
  - Purpose: Filter known genes using SQANTI categories.

- **flairMODULE6.withSRsupport.sh**:
  - Purpose: Bash script to run FLAIR diffSplice - module 6 (example for short-read supported data).

- **flair.diff.iso.usage.KNOWN.filteredDATA.bash.sh**:
  - Purpose: Bash script to run FLAIR different isoform usage tool (example for known genes filtered data).

- **diffsplice.fishers.exact.sh**:
  - Purpose: Compute Fisher test on splicing events.

- **flair.align.singlebams.sh**:
  - Purpose: From long-read fastq, generates single BAMs for each sample using FLAIR align module 1.

- **run.ggsashimi.sh**:
  - Purpose: Uses sushi function to generate ggsashimi plots.
