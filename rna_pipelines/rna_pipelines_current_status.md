# RNA pipelines - current status

Created: 2020/09/24 15:19:42
Last modified: 2020/09/24 16:18:06

- **Aim:** Evaluate the current pipelines available for processing RNA-seq data. This will help us decide if there is a pipeline currently available for our use, one we could adapt, or if we will need to create an RNA-seq pipeline from scratch

## Table of contents

- [RNA pipelines - current status](#rna-pipelines---current-status)
  - [Table of contents](#table-of-contents)
  - [Overview](#overview)
  - [RNA pipelines currently available](#rna-pipelines-currently-available)
    - [Snakemake](#snakemake)

## Overview


A summary of some of the RNA-seq data analysis tools by [Ruairi J MacKenzie, 2018](https://www.technologynetworks.com/genomics/articles/rna-seq-basics-applications-and-protocol-299461):

"Tools like Sailfish, RSEM and BitSeq12 will help you quantify your expression levels, whilst tools like MISO, which quantifies alternatively spliced genes, are available for more specialized analysis13. There is a library of these tools out there, and reading reviews and roundups are your best way to find the right tool for your research."

## RNA pipelines currently available

### Snakemake

- VIPER: Visualization Pipeline for RNA-seq
  - Paper [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2139-9) and github repo [here](https://github.com/hanfeisun/viper-rnaseq)
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
  - Package management with Conda
  - This pipeline focuses on data visualisation, using Snakemake, parallelization/speed, ease of use (like less options, easier install)
  - They say they differ from other pipelines available by "the number of features included, package management software, and reporting functionalities"
  - They say they focus on the best practice software rather than including many options (like HppRNA)
  - See their comparison to other pipelines [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2139-9/tables/1)

- TRAPLINE
  - Package management with Galaxy

- HppRNA
