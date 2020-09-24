# RNA pipelines - current status

Created: 2020/09/24 15:19:42
Last modified: 2020/09/24 17:48:28

- **Aim:** Evaluate the current pipelines available for processing RNA-seq data. This will help us decide if there is a pipeline currently available for our use, one we could adapt, or if we will need to create an RNA-seq pipeline from scratch

## Table of contents

- [RNA pipelines - current status](#rna-pipelines---current-status)
  - [Table of contents](#table-of-contents)
  - [Overview](#overview)
  - [RNA pipelines currently available](#rna-pipelines-currently-available)
    - [VIPER: Visualization Pipeline for RNA-seq](#viper-visualization-pipeline-for-rna-seq)
    - [TRAPLINE](#trapline)
    - [HppRNA](#hpprna)
    - [DRAGEN RNA Pipeline](#dragen-rna-pipeline)
    - [rna-seq-star-deseq2](#rna-seq-star-deseq2)
    - [rna-seq-kallisto-sleuth](#rna-seq-kallisto-sleuth)
    - [single-cell-rna-seq](#single-cell-rna-seq)
    - [single-cell-drop-seq](#single-cell-drop-seq)
    - [National Cancer Institute mRNA quantification analysis pipeline](#national-cancer-institute-mrna-quantification-analysis-pipeline)
    - [CloseCall](#closecall)
    - [snDrop_prep](#sndrop_prep)
    - [RNAsik-pipe](#rnasik-pipe)
    - [exceRpt](#excerpt)
    - [GEO2RNAseq](#geo2rnaseq)

## Overview


A summary of some of the RNA-seq data analysis tools by [Ruairi J MacKenzie, 2018](https://www.technologynetworks.com/genomics/articles/rna-seq-basics-applications-and-protocol-299461):

"Tools like Sailfish, RSEM and BitSeq12 will help you quantify your expression levels, whilst tools like MISO, which quantifies alternatively spliced genes, are available for more specialized analysis13. There is a library of these tools out there, and reading reviews and roundups are your best way to find the right tool for your research."

## RNA pipelines currently available

### VIPER: Visualization Pipeline for RNA-seq

- Paper [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2139-9) and github repo [here](https://github.com/hanfeisun/viper-rnaseq)
- Open source
- Snakemake
- Package management with Conda
- Hasn't been updated in four years - may not be actively maintained/supported
- [STAR](https://github.com/alexdobin/STAR) for alignment
- [Cufflinks](https://github.com/cole-trapnell-lab/cufflinks) for transcript assembly, normalisation of read counts etc.
- bams in BigWig format to allow visualisation in genome browser
- [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki) for fusion gene discovery for paired-end data
- Read quality metrics
- Differential expression and pathway analysis - can use [DEseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) or build in a preferred differential expression method
- Immunology analysis module available
- Whole-genome SNV calling (human and mouse) module available
- Viral analysis (human samples only) module available
- Batch effect correction module available
- This pipeline focuses on data visualisation, using Snakemake, parallelization/speed, ease of use (such as less options, easier to install)
- They say they differ from two other pipelines available by "the number of features included, package management software, and reporting functionalities"
- They say they focus on the best practice software rather than including many options (like HppRNA)
- See their comparison to other pipelines [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2139-9/tables/1)

### TRAPLINE

  - Package management with Galaxy

### HppRNA

- github repo [here](https://github.com/NextGenBioinformatics/hppRNA)
- RNA-Seq analysis of numerous samples
- They describe as "parameter-free"
- Hasn't been updated in two years - may not be actively maintained/supported

### DRAGEN RNA Pipeline

- Website [here](https://www.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/edico-genome-inc-dragen-rna-pipeline.html)
- Not open source?
- For research use only

### rna-seq-star-deseq2

- github repo [here](https://github.com/snakemake-workflows/rna-seq-star-deseq2)
- Open source
- Snakemake
- Package management with Conda
- Differential expression analysis with [STAR](https://github.com/alexdobin/STAR) and [Deseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

### rna-seq-kallisto-sleuth

- github repo [here](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth)
- Open source
- Snakemake
- Package management with Conda
- Differential expression analysis with [Kallisto](https://pachterlab.github.io/kallisto/) and [Sleuth](https://pachterlab.github.io/sleuth/).


### single-cell-rna-seq

- github repo [here](https://github.com/snakemake-workflows/single-cell-rna-seq)
- Open source
- Snakemake
- Package management with Conda
- Single cell RNA-seq workflow

### single-cell-drop-seq

- github repo [here](https://github.com/snakemake-workflows/single-cell-drop-seq)
- Open source
- Snakemake
- Package management with Conda
- Single cell RNA-seq workflow
- [STAR](https://github.com/alexdobin/STAR) for alignment

### [National Cancer Institute](https://www.cancer.gov/) mRNA quantification analysis pipeline

- Website [here](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#mrna-analysis-pipeline)
- [STAR](https://github.com/alexdobin/STAR) for alignment
- GDC gene fusion pipeline using [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki)
- Arriba gene fusion pipeline

### CloseCall

- github repo [here](https://github.com/StevenWingett/CloseCall) and paper [here](https://www.nature.com/articles/s41597-020-0372-3)
- Pipeline for processing RNA-RNA proximity data
- Open source
- Mapping and QC
- Monte Carlo Simulation to identify statistically significant RNA-RNA interactions
- Actively maintained (last commit 13 months ago)
- Written in perl, java and R (not in a workflow language?)

### snDrop_prep

- github repo [here](https://github.com/chensong611/snDrop_prep) and paper [here](https://www.nature.com/articles/s41467-019-10861-2)
- Open source
- Single-nucleus RNA-sequencing pipeline
- Python and bash (not in a workflow language)

### RNAsik-pipe

- github [here](https://github.com/MonashBioinformaticsPlatform/RNAsik-pipe) and website [here](https://monashbioinformaticsplatform.github.io/RNAsik-pipe/)
- Open source
- Written in BigDataScript (bds) (with underlying Java) (not in a workflow language)

### exceRpt

- github [here](https://github.com/gersteinlab/exceRpt), website [here](https://rna-seqblog.com/excerpt-a-comprehensive-analytic-platform-for-extracellular-rna-profiling/)
- Open source
- Extracellular RNA profiling
- Hasn't been updated in two years - may not be actively maintained/supported
- Written in java, R, bash, pl (not in a workflow language)

### GEO2RNAseq

- Hosted on Anaconda [here](https://anaconda.org/xentrics/r-geo2rnaseq), paper [here](https://www.biorxiv.org/content/10.1101/771063v1.full)
- Open source