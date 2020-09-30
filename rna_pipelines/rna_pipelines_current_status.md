# RNA pipelines - current status

Created: 2020/09/24 15:19:42
Last modified: 2020/09/30 16:13:32

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
    - [GeneTEFlow](#geneteflow)
    - [Parabricks](#parabricks)
    - [DiMSum](#dimsum)
    - [scTyper](#sctyper)
    - [RiboMiner](#ribominer)
    - [SLFinder](#slfinder)
    - [Cell Ranger pipeline (10×Genomics)](#cell-ranger-pipeline-10genomics)
    - [Seurat](#seurat)
  - [Other](#other)
  - [Small non-coding RNA-Seq pipelines currently available](#small-non-coding-rna-seq-pipelines-currently-available)
    - [smrnaseq](#smrnaseq)
    - [sports1.1](#sports11)
    - [umitools](#umitools)
    - [smallseq](#smallseq)
    - [RNA-seq](#rna-seq)
    - [snakemake_RNA-seq](#snakemake_rna-seq)
    - [snakemake_fastqc](#snakemake_fastqc)
    - [snakemake_deeptools](#snakemake_deeptools)
    - [SnapT](#snapt)
    - [short-ncrna-annotation](#short-ncrna-annotation)
    - [ncRNA_Pipeline](#ncrna_pipeline)
    - [FlaiMapper](#flaimapper)
    - [irap](#irap)
    - [ymp](#ymp)
    - [rnaseq-pipeline](#rnaseq-pipeline)
    - [exceRNApipeline](#excernapipeline)
    - [sRNAflow](#srnaflow)
    - [gorap](#gorap)
  - [Other things of possible interest](#other-things-of-possible-interest)
  - [Notes](#notes)

## Overview


A summary of some of the RNA-seq data analysis tools by [Ruairi J MacKenzie, 2018](https://www.technologynetworks.com/genomics/articles/rna-seq-basics-applications-and-protocol-299461):

"Tools like Sailfish, RSEM and BitSeq12 will help you quantify your expression levels, whilst tools like MISO, which quantifies alternatively spliced genes, are available for more specialized analysis. There is a library of these tools out there, and reading reviews and roundups are your best way to find the right tool for your research."

Note that there also exists an R framework for comparing pipelines (see [here](chrome-extension://oemmndcbldboiebfnladdacbdfmadadm/https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-02136-7))

Papers that compare single-cell RNA pipelines (from [this paper](chrome-extension://oemmndcbldboiebfnladdacbdfmadadm/https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-02136-7)):

- Cobos FA, Alquicira-Hernandez J, Powell J, Mestdagh P, De Preter K. Comprehensive benchmarking of computational deconvolution of transcriptomics data. bio Rxiv. 2020.https://doi.org/10.1101/2020.01.10.897116.T.4.   
- Cole MB, Risso D, Wagner A, DeTomaso D, Ngai J, Purdom E, Dudoit S, Yosef N. Performance assessment and selection of normalization procedures for single-cell RNA-Seq. Cell Syst. 2019;8(4):315–28.https://doi.org/10.1016/j.cels.2019.03.010.5.   
- Dal Molin A, Baruzzo G, Di Camillo B. Single-cell RNA-sequencing: Assessment of differential expression analysis methods. Front Genet. 2017;8(62):.https://doi.org/10.3389/fgene.2017.00062.6.   
- Duo A, Robinson MD, Soneson C. A systematic performance evaluation of clustering methods for single-cell RNA-seq data. F1000Research. 2018;7:1141.https://doi.org/10.12688/f1000research.15666.2.
- Freytag S, Tian L, Lönnstedt I, Ng M, Bahlo M. Comparison of clustering tools in R for medium-sized 10x Genomics single-cell RNA-sequencing data. F1000Research. 2018;7(1297):1–29.https://doi.org/10.12688/f1000research.15809.2.8.   
- Gao M, Ling M, Tang X, Wang S, Xiao X, Qiao Y, Yang W, Yu R. Comparison of high-throughput single-cell RNA sequencing data processing pipelines. bioRxiv. 2020.https://doi.org/10.1101/2020.02.09.940221.9.   - Heiser CN, Lau KS. A quantitative framework for evaluating single-cell data structure preservation by dimensionality reduction techniques. bioRxiv. 2019684340.https://doi.org/10.1101/684340.10. 
- Hou W, Ji Z, Ji H, Hicks SC. A systematic evaluation of single-cell RNA-sequencing imputation methods. bioRxiv.2020.https://doi.org/10.1101/2020.01.29.925974.11.  
- Jaakkola MK, Seyednasrollah F, Mehmood A, Elo LL. Comparison of methods to detect differentially expressed genes between single-cell populations. Brief Bioinforma. 2017;18(5):735–43.https://doi.org/10.1093/bib/bbw057.12.  
- Krzak M, Raykov Y, Boukouvalas A, Cutillo L, Angelini C. Benchmark and parameter sensitivity analysis of single-cellRNA sequencing clustering methods. Front Genet. 2019;10:1253.https://doi.org/10.3389/fgene.2019.01253.13.  
- Soneson C, Robinson MD. Bias, robustness and scalability in single-cell differential expression analysis. Nat Methods.2018;15(4):255–61.https://doi.org/10.1038/nmeth.4612.14.  
- Sun S, Zhu J, Ma Y, Zhou X. Accuracy, robustness and scalability of dimensionality reduction methods for single-cell RNA-seq analysis. Genome Biol. 2019;20(269):1–21.https://doi.org/10.1186/s13059-019-1898-6.15.  
- Tian L, Dong X, Freytag S, Le Cao K-A, Su S, Amann-Zalcenstein D, Weber TS, Seidi A, Naik S, Ritchie ME.scRNA-seq mixology: towards better benchmarking of single cell RNA-seq protocols and analysis methods. bioRxiv.2018433102.https://doi.org/10.1101/433102.16.  
- Tran HTN, Ang KS, Chevrier M, Zhang X, Lee NYS, Goh M, Chen J. A benchmark of batch-effect correction methods for single-cell RNA sequencing data,. Genome Biol. 2020;21(1):1–32.https://doi.org/10.1186/s13059-019-1850-9.17.  
- Tsuyuzaki K, Sato H, Sato K, Nikaido I. Benchmarking principal component analysis for large-scale single-cellRNA-sequencing. Genome Biol. 2020;21(9):1–17.https://doi.org/10.1186/s13059-019-1900-3.18.  
- Vieth B, Parekh S, Ziegenhain C, Enard W, Hellmann I. A systematic evaluation of single cell RNA-seq analysis pipelines. Nat Commun. 2019;10(1):1–11.https://doi.org/10.1038/s41467-019-12266-7.19.  
- Wang T, Li B, Nelson CE, Nabavi S. Comparative analysis of differential gene expression analysis tools for single-cellRNA sequencing data. BMC Bioinformatics. 2019;20.1(40):1–16.https://doi.org/10.1186/s12859-019-2599-6.20.  
- Yip SH, Sham PC, Wang J. Evaluation of tools for highly variable gene discovery from single-cell RNA-seq data. Brief Bioinforma. 2018;20(4):1583–9.https://doi.org/10.1093/bib/bby011.21.  
- Zhang L, Zhang S. Comparison of computational methods for imputing single-cell RNA-sequencing data. IEEE/ACMTrans Comput Biol Bioinforma. 2018.https://doi.org/10.1109/TCBB.2018.2848633.22.

A [recent paper](https://www.embopress.org/doi/full/10.15252/msb.20188746) (2019) that describes the best practices in single-cell RNA-seq analysis

A [recent paper](AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAqUwggKhBgkqhkiG9w0BBwagggKSMIICjgIBADCCAocGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMdbpTwGltlCc7WTjbAgEQgIICWEvhDFT) (2018) that describes a website that acts as a referece set/database of human long non-coding RNAs

A [recent paper](chrome-extension://oemmndcbldboiebfnladdacbdfmadadm/https://www.nature.com/articles/s12276-018-0071-8.pdf) (2018) that describes Single-cell RNA sequencing technologiesand bioinformatics pipelines

## RNA pipelines currently available

### VIPER: Visualization Pipeline for RNA-seq

- Paper [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2139-9) and github repo [here](https://github.com/hanfeisun/viper-rnaseq)
- Open source
- Workflow language - snakemake
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
- Workflow language - snakemake
- Package management with Conda
- Differential expression analysis with [STAR](https://github.com/alexdobin/STAR) and [Deseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

### rna-seq-kallisto-sleuth

- github repo [here](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth)
- Open source
- Workflow language - snakemake
- Package management with Conda
- Differential expression analysis with [Kallisto](https://pachterlab.github.io/kallisto/) and [Sleuth](https://pachterlab.github.io/sleuth/)

### single-cell-rna-seq

- github repo [here](https://github.com/snakemake-workflows/single-cell-rna-seq)
- Open source
- Workflow language - snakemake
- Package management with Conda
- Single cell RNA-seq workflow

### single-cell-drop-seq

- github repo [here](https://github.com/snakemake-workflows/single-cell-drop-seq)
- Open source
- Workflow language - snakemake
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

### GeneTEFlow

- github [here](https://github.com/zhongw2/GeneTEFlow), paper [here](chrome-extension://oemmndcbldboiebfnladdacbdfmadadm/https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0232994&type=printable)
- Differential expression analysis of genes and locus-specific transposable elements from RNA sequencing
- Open source
- Nextflow
- [Trimmomatic](https://github.com/timflutre/trimmomatic) for fastq trimming
- Actively maintained (last commit 3 months ago)

### Parabricks

- Some RNA tools available in the GPU-accelerated toolkit
- Not open source
- [Star](https://www.nvidia.com/en-us/docs/parabricks/star/) for alignment
- [Star-Fusion](https://www.nvidia.com/en-us/docs/parabricks/star-fusion/) for fusion gene discovery for paired-end data

### DiMSum

- github [here](https://github.com/lehner-lab/DiMSum), paper [here](chrome-extension://oemmndcbldboiebfnladdacbdfmadadm/https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-02091-3)
- Deep mutational scanning (DMS) enabling the multiplexed measurement of the effects of thousands of variants of proteins, RNAs, and regulatory elements
- Open source
- Written in R (not in a workflow language)

### scTyper

- github [here](https://github.com/omicsCore/scTyper), paper [here](chrome-extension://oemmndcbldboiebfnladdacbdfmadadm/https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/s12859-020-03700-5)
- Cell typing analysis of single-cell RNA-seq data
- Written in R (not in a workflow language)

### RiboMiner

- github [here](https://github.com/xryanglab/RiboMiner), paper [here](chrome-extension://oemmndcbldboiebfnladdacbdfmadadm/https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/s12859-020-03670-8)
- Open source
- Mining multi-dimensional features of the translatome with ribosome profiling data
- Four function parts: Quality Control, Metagene Analysis, Feature Analysis, Enrichment Analysis:
- Written in python (not a pipeline, not in a workflow language)

### SLFinder

- github [here](https://github.com/LBC-Iriarte/SLFinder), paper [here](chrome-extension://oemmndcbldboiebfnladdacbdfmadadm/https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/s12859-020-03610-6)
- Open source
- Novel identification of splice-leader sequences
- Actively maintained (last commit 3 months ago)
- Written in bash (not in a workflow language)

### Cell Ranger pipeline (10×Genomics)
- website [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)
- Looks like it's open source
- A set of analysis pipelines that process Chromium single-cell RNA-seq output to align reads, generate feature-barcode matrices and perform clustering and gene expression analysis
- t-SNE is implemented (according to [this paper](chrome-extension://oemmndcbldboiebfnladdacbdfmadadm/https://www.nature.com/articles/s12276-018-0071-8.pdf))

### Seurat
- website [here](https://satijalab.org/seurat/)
- QC, analysis, and exploration of single-cell RNA-seq data
- Written in R (not in a workflow language)
- t-SNE is implemented (according to [this paper](chrome-extension://oemmndcbldboiebfnladdacbdfmadadm/https://www.nature.com/articles/s12276-018-0071-8.pdf))

## Other

- [HPC-REDItools - A tool for large-scale RNA-editing analysis](chrome-extension://oemmndcbldboiebfnladdacbdfmadadm/https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/s12859-020-03562-x) (RNA-editing is a molecular process through which some cells can make discrete changes to specific nucleotide sequences within an RNA molecule after it has been generated by RNA polymerase)
- Open source

## Small non-coding RNA-Seq pipelines currently available

### smrnaseq

- github [here](https://github.com/nf-core/smrnaseq) and paper [here](https://www.biorxiv.org/content/10.1101/610741v1)
- Best-practice analysis pipeline used for small RNA sequencing data
- Open source
- Nextflow
- Deployable to SLURM and AWS
- Actively maintained (last commit 12 months ago)

### sports1.1

- github [here](https://github.com/junchaoshi/sports1.1) and paper [here](https://www.sciencedirect.com/science/article/pii/S1672022918300445)
- Small non-coding RNA annotation Pipeline Optimized for rRNA- and tRNA-Derived Small RNAs
- Open source
- Written in R and perl (not in a workflow language)
- Actively maintained (last commit 2 months ago)

### umitools

- github [here](https://github.com/weng-lab/umitools) and paper [here](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4933-1)
- A toolset for handling sequencing data with unique molecular identifiers (UMIs)
- Unique molecular identifiers (UMIs) are a type of molecular barcoding that provides error correction and increased accuracy during sequencing. These molecular barcodes are short sequences used to uniquely tag each molecule in a sample library
- Open source
- Written in python (not in a workflow language)
- Hasn't been updated in two years - may not be actively maintained/supported

### smallseq

- github [here](https://github.com/eyay/smallseq)
- Analyze small RNAs from single-cells
- Open source
- Written in python (not in a workflow language)
- Hasn't been updated in three years - may not be actively maintained/supported

### RNA-seq

- github [here](https://github.com/biowdl/RNA-seq)
- A Biowdl workflows usable for processing RNA-seq data. This pipeline will performs QC (including adapter clipping), mapping, variant-calling and expression quantification
- Open source
- Workflow Description Language
- Actively maintained (last commit 20 days ago)

### snakemake_RNA-seq

- github [here](https://github.com/WilliamJeong2/snakemake_RNA-seq)
- A snakemake pipeline for the analysis of RNA-seq data that makes use of hisat2 and Stringtie
- Snakemake
- Open source
- Actively maintained (last commit 25 days ago)

### snakemake_fastqc

- github [here](https://github.com/AngryMaciek/snakemake_fastqc)
- A small snakemake workflow for FastQC
- Snakemake
- Open source
- Actively maintained (last commit 25 days ago)

### snakemake_deeptools

- github [here](https://github.com/AngryMaciek/snakemake_deeptools)
- Snakemake pipeline for RNA-Seq data analysis with deepTools
- Snakemake
- Open source
- Actively maintained (last commit 25 days ago)

### SnapT

- github [here](https://github.com/ursky/SnapT)
- A Small non-coding RNA annotation pipeline for Transcriptomic or metatranscriptomic data
- Written in python (not in a workflow language)
- Open source
- Actively maintained (last commit 6 months ago)

### short-ncrna-annotation

- github [here](https://github.com/SimonSchafferer/short-ncrna-annotation)
- Annotation of commonly used interval/range based data such as the browser extended display format (BED), or the general feature format (GFF). In addition, it provides sequence based annotation, by employing the NCBI blast+ software
- Written in R (not in a workflow language)
- Open source
- Hasn't been updated in six years - may not be actively maintained/supported

### ncRNA_Pipeline

- github [here](https://github.com/navygit/ncRNA_Pipeline)
- Open source
- Written in perl (not in a workflow language)
- Hasn't been updated in five years - may not be actively maintained/supported

### FlaiMapper

- github [here], paper [here](https://academic.oup.com/bioinformatics/article/31/5/665/2748143)
- Annotation of small ncRNA-derived fragments using RNA-seq high-throughput data
- Open source
- Written in bash (not in a workflow language)
- Hasn't been updated in four years - may not be actively maintained/supported

### irap

- github [here], paper [here](https://academic.oup.com/bioinformatics/article/31/5/665/2748143)
- Flexible RNA-seq analysis pipeline that allows the user to select and apply their preferred combination of existing tools for mapping reads, quantifying expression and testing for differential expression
- Open source
- Written in R, perl, bash (not in a workflow language)
- Hasn't been updated in two years - may not be actively maintained/supported

### ymp

- github [here](https://github.com/epruesse/ymp)
- Flexible omics pipeline (QC, trimming, contaminant removal), assemble metagenomes, annotate assemblies, or assemble and quantify RNA-Seq transcripts, offering a choice of tools for each of those procecssing stages.
- Written in python (not in a workflow language)
- Open source
- Actively maintained (last commit a month ago)

### rnaseq-pipeline

- github [here](https://github.com/PavlidisLab/rnaseq-pipeline)
- Written in python, R (not in a workflow language)
- Open source
- Actively maintained (last commit 3 months ago)

### exceRNApipeline

- github [here](https://github.com/zhuchcn/exceRNApipeline)
- Data processing pipeline for extracellular small RNA-seq from human specimen.The pipeline is designated for running on HPC with the job management system SLURM.
- Preprocess: remove adapters and trimming low quality nucleotides, using HTStream.
- UniVec: map to the NCBI's UniVec database to remove vector origins
- RiboRNA: map to the rRNA sequences
- Human Genome: map to human genome
- Repetitive Elements: map to RepeatMasker's repetitive elements sequences
- SILVA: map to SILVA's ribosomal rRNA gene of bacteria, archaea, and fungi
- Bacteria: map to all bacteria genomes in ensemble's
- Only supports single end sequencing data
- Snakemake
- Open source
- Actively maintained (last commit 6 months ago)

### sRNAflow

- github [here](https://github.com/zajakin/sRNAflow)
- Analysis of small RNA that fulfills the specific needs for samples derived from biofluids
- Nextflow
- Open source
- Actively maintained (last commit 28 days ago)

### gorap

- github [here](https://github.com/koriege/gorap), website [here](https://www.rna.uni-jena.de/research/software/)
- Screens genomic sequences for all non-coding RNAs present in the Rfam database
- Provides ncRNA based reconstruction of phylogenetic trees and is able to perform de novo predictions including TPM calculations from RNA-Seq experiments
- Open source
- Written in perl (not in a workflow language)
- Actively maintained (last commit 7 months ago)

## Other things of possible interest

- Prediction and annotation of tRNA-derived small RNAs [here](https://github.com/wangqinhu/tsRFinder)
- A JBrowse plugin to support small RNA visualization [here](https://github.com/bhofmei/jbplugin-smallrna)
- A set of tools related to the Rfam (database of non-coding RNA families https://rfam.org/) production pipeline [here](https://github.com/Rfam/rfam-production)
- Script to map annotated ncRNAs (and other elements) to a chromosome and visualise [here](https://github.com/fanagislab/draw_annotation/tree/master/bin)

## Notes

- There seems to be a number of RNA-seq pipelines currently available (open access) on github\
- There seems to very few pipelines dedicated to small non-coding RNA-seq
- 
