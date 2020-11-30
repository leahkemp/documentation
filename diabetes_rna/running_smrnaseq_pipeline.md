# Running smrnaseq pipeline

Created: 2020/11/13 12:29:25
Last modified: 2020/11/17 17:01:03

- **Aim:** In [this document](./rna_pipelines_current_status.md) I settled on using the [smrnaseq](https://github.com/nf-core/smrnaseq) nextflow pipeline to process our small non-coding RNA-seq data. This document documents/describes the process of trying this pipeline out on the data (extending on Miles Bentons work). Thi is part of the wider [small rnaseq hepatic portal project](./project_notes_small_rnaseq_hepatic_portal.md)
- **Prerequisite software:** [conda 4.9.0](https://docs.conda.io/en/latest/), [git 2.7.4](https://git-scm.com/)
- **OS:** Ubuntu 16.04 (Wintermute - research server)

## Table of contents

- [Running smrnaseq pipeline](#running-smrnaseq-pipeline)
  - [Table of contents](#table-of-contents)
  - [Overview of data](#overview-of-data)
  - [Run data through pipeline](#run-data-through-pipeline)
    - [Clone pipeline](#clone-pipeline)
    - [Install nextflow into a conda env](#install-nextflow-into-a-conda-env)
    - [Run the pipeline](#run-the-pipeline)
  - [Explore the outputs!](#explore-the-outputs)
    - [1. Raw read QC - FastQC](#1-raw-read-qc---fastqc)
      - [i. Insert Size calculation](#i-insert-size-calculation)
      - [ii. Collapse reads seqcsluter](#ii-collapse-reads-seqcsluter)
    - [2. Adapter trimming - Trim Galore!](#2-adapter-trimming---trim-galore)
  - [Re-run data through pipeline](#re-run-data-through-pipeline)
  - [Explore the outputs!](#explore-the-outputs-1)
    - [1. Raw read QC - FastQC](#1-raw-read-qc---fastqc-1)
      - [i. Insert Size calculation](#i-insert-size-calculation-1)
      - [ii. Collapse reads seqcsluter](#ii-collapse-reads-seqcsluter-1)
    - [2. Adapter trimming - Trim Galore!](#2-adapter-trimming---trim-galore-1)
    - [3. Alignment against miRBase mature miRNA - Bowtie1](#3-alignment-against-mirbase-mature-mirna---bowtie1)
    - [4. Alignment against miRBase hairpin](#4-alignment-against-mirbase-hairpin)
      - [Unaligned reads from step 3 (Bowtie1)](#unaligned-reads-from-step-3-bowtie1)
      - [Collapsed reads from step 2.2 (Bowtie1)](#collapsed-reads-from-step-22-bowtie1)
    - [5. Post-alignment processing of miRBase hairpin](#5-post-alignment-processing-of-mirbase-hairpin)
      - [i. Basic statistics from step 3 and step 4.1 (SAMtools)](#i-basic-statistics-from-step-3-and-step-41-samtools)
      - [ii. Analysis on miRBase hairpin counts (edgeR)](#ii-analysis-on-mirbase-hairpin-counts-edger)
      - [iii. miRNA and isomiR annotation from step 4.1 (mirtop)](#iii-mirna-and-isomir-annotation-from-step-41-mirtop)
    - [6. Alignment against host reference genome (Bowtie1)](#6-alignment-against-host-reference-genome-bowtie1)
      - [i. Post-alignment processing of alignment against host reference genome (SAMtools)](#i-post-alignment-processing-of-alignment-against-host-reference-genome-samtools)
    - [7. miRNA quality control (mirtrace)](#7-mirna-quality-control-mirtrace)
    - [8. Present QC for raw read, alignment, and expression results (MultiQC)](#8-present-qc-for-raw-read-alignment-and-expression-results-multiqc)

## Overview of data

```bash
cd /store/lkemp/smrnaseq_hps/
```

Have a look at the input files

```bash
ls -lh ./fastq/
```

Output:

```bash
total 18G
-rw-rw-r-- 1 lkemp lkemp 230M Oct 13 13:14 HPS004_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 212M Oct 13 13:14 HPS022_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 228M Oct 13 13:14 HPS030_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 192M Oct 13 13:14 HPS033_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 223M Oct 13 13:14 HPS041_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 170M Oct 13 13:14 HPS048_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 175M Oct 13 13:14 HPS061_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 200M Oct 13 13:14 HPS064_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 197M Oct 13 13:14 HPS070_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 235M Oct 13 13:14 HPS071_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 219M Oct 13 13:14 HPS072_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 218M Oct 13 13:14 HPS080_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 357M Oct 13 13:14 HPS081_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 398M Oct 13 13:14 HPS083_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 363M Oct 13 13:14 HPS086_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 350M Oct 13 13:14 HPS100_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 403M Oct 13 13:14 HPS108_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 386M Oct 13 13:14 HPS114_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 379M Oct 13 13:14 HPS116_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 333M Oct 13 13:14 HPS119_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 416M Oct 13 13:14 HPS121_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 376M Oct 13 13:14 HPS123_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 354M Oct 13 13:16 HPS124_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 411M Oct 13 13:14 HPS130_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 324M Oct 13 13:16 HPS139_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 293M Oct 13 13:17 HPS145_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 332M Oct 13 13:17 HPS146_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 343M Oct 13 13:17 HPS147_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 361M Oct 13 13:17 HPS150_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 351M Oct 13 13:17 HPS152_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 302M Oct 13 13:17 HPS155_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 357M Oct 13 13:17 HPS163_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 283M Oct 13 13:17 HPS172_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 298M Oct 13 13:17 HPS175_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 382M Oct 13 13:17 HPS183_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 328M Oct 13 13:17 HPS184_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 342M Oct 13 13:17 HPS186_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 334M Oct 13 13:17 HPS193_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 370M Oct 13 13:17 HPS199_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 381M Oct 13 13:17 HPS200_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 281M Oct 13 13:17 HPS203_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 393M Oct 13 13:17 HPS208_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 309M Oct 13 13:17 HPS210_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 330M Oct 13 13:17 HPS217_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 327M Oct 13 13:17 HPS223_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 380M Oct 13 13:17 HPS227_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 380M Oct 13 13:17 HPS230_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 329M Oct 13 13:17 HPS234_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 206M Oct 13 13:17 HPS240_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 214M Oct 13 13:17 HPS255_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 206M Oct 13 13:17 HPS257_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 194M Oct 13 13:17 HPS259_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 204M Oct 13 13:17 HPS263_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 192M Oct 13 13:17 HPS267_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 194M Oct 13 13:17 HPS269_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 201M Oct 13 13:17 HPS270_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 241M Oct 13 13:17 HPS281_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 226M Oct 13 13:17 HPS283_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 244M Oct 13 13:17 HPS286_combined.fastq.gz
-rw-rw-r-- 1 lkemp lkemp 187M Oct 13 13:17 HPS296_combined.fastq.gz
```

Get number of input fastq files for pipeline

```bash
ls ./fastq/ -1 | wc -l
```

Output:

```bash
60
```

## Run data through pipeline

### Clone pipeline

```bash
cd /store/lkemp/smrnaseq_hps/
git clone https://github.com/nf-core/smrnaseq.git
```

### Install nextflow into a conda env

```bash
conda create -n nextflow python=3.7.6
conda activate nextflow
mamba install -c bioconda nextflow=20.07.1
```

### Run the pipeline

`--min-length 17` set to avoid small smRNA's mapping to many locations

```bash
cd /store/lkemp/smrnaseq_hps/

nextflow run /store/lkemp/smrnaseq_hps/smrnaseq/main.nf --reads 'fastq/*_combined.fastq.gz' \
  -profile conda --protocol illumina --genome 'GRCh37' \
  --saveReference -resume --min_length 17
```

Worked great!

Print pipeline report

```bash
cat ./results/pipeline_info/pipeline_report.txt
```

Output:

```bash
========================================
 nf-core/smrnaseq v1.0.0
========================================
Run Name: chaotic_ptolemy

## nf-core/smrnaseq execution completed successfully! ##


The workflow was completed at 2020-11-13T12:14:45.729368+13:00 (duration: 28m 20s)

The command used to launch the workflow was as follows:

  nextflow run /store/lkemp/smrnaseq_hps/smrnaseq/main.nf --reads 'fastq/*_combined.fastq.gz' -profile conda --protocol illumina --genome GRCh37 --saveReference -resume --min_length 17



Pipeline Configuration:
-----------------------
 - Run Name: chaotic_ptolemy
 - Reads: fastq/*_combined.fastq.gz
 - Genome: GRCh37
 - Min Trimmed Length: 17
 - Trim 5' R1: 0
 - Trim 3' R1: 0
 - miRBase mature: s3://ngi-igenomes/igenomes//Homo_sapiens/Ensembl/GRCh37/Annotation/SmallRNA/mature.fa
 - miRBase hairpin: s3://ngi-igenomes/igenomes//Homo_sapiens/Ensembl/GRCh37/Annotation/SmallRNA/hairpin.fa
 - Bowtie Index for Ref: s3://ngi-igenomes/igenomes//Homo_sapiens/Ensembl/GRCh37/Sequence/BowtieIndex/genome
 - Save Reference: Yes
 - Protocol: illumina
 - miRTrace species: hsa
 - 3' adapter: TGGAATTCTCGGGTGCCAAGG
 - Output dir: ./results
 - Launch dir: /store/lkemp/smrnaseq_hps
 - Working dir: /store/lkemp/smrnaseq_hps/work
 - Current home: /home/lkemp
 - Current user: lkemp
 - Current path: /store/lkemp/smrnaseq_hps
 - Script dir: /store/lkemp/smrnaseq_hps/smrnaseq
 - Config Profile: conda
 - Max Resources: 128 GB memory, 16 cpus, 10d time per job
 - User: lkemp
 - Date Started: 2020-11-13T11:46:25.964109+13:00
 - Date Completed: 2020-11-13T12:14:45.729368+13:00
 - Pipeline script file path: /store/lkemp/smrnaseq_hps/smrnaseq/main.nf
 - Pipeline script hash ID: eba0f9b76b0f4757e88b4463a5c7ae1d
 - Docker image: nfcore/smrnaseq:1.0.0
 - Nextflow Version: 20.07.1
 - Nextflow Build: 5412
 - Nextflow Compile Timestamp: 24-07-2020 15:18 UTC


--
nf-core/smrnaseq is a bioinformatics best-practice analysis pipeline used for small RNA sequencing data at the National Genomics Infrastructure at SciLifeLab Stockholm, Sweden.
The pipeline uses Nextflow, a bioinformatics workflow tool. It pre-processes raw data from FastQ inputs, aligns the reads and performs extensive quality-control on the results.
For more information, please see the pipeline homepage: https://github.com/nf-core/smrnaseq
```

## Explore the outputs!

### 1. Raw read QC - FastQC

Get sequence lengths for all reads *before* trimming

```bash
mkdir ./results/fastqc_unzipped/

# Unzip fastqc files
for zip in ./results/fastqc/*_combined_fastqc.zip; do
unzip $zip -d ./results/fastqc_unzipped/
done

# Get sequence lengths before trimming for all samples from fastqc summary files
for summaryfile in ./results/fastqc_unzipped/HPS*_combined_fastqc/fastqc_data.txt; do
grep -e 'Filename' -e 'Sequence length' $summaryfile
done
```

Output:

```bash
Filename        HPS004_combined.fastq.gz
Sequence length 100
Filename        HPS022_combined.fastq.gz
Sequence length 100
Filename        HPS030_combined.fastq.gz
Sequence length 100
Filename        HPS033_combined.fastq.gz
Sequence length 100
Filename        HPS041_combined.fastq.gz
Sequence length 100
Filename        HPS048_combined.fastq.gz
Sequence length 100
Filename        HPS061_combined.fastq.gz
Sequence length 100
Filename        HPS064_combined.fastq.gz
Sequence length 100
Filename        HPS070_combined.fastq.gz
Sequence length 100
Filename        HPS071_combined.fastq.gz
Sequence length 100
Filename        HPS072_combined.fastq.gz
Sequence length 100
Filename        HPS080_combined.fastq.gz
Sequence length 100
Filename        HPS081_combined.fastq.gz
Sequence length 100
Filename        HPS083_combined.fastq.gz
Sequence length 100
Filename        HPS086_combined.fastq.gz
Sequence length 100
Filename        HPS100_combined.fastq.gz
Sequence length 100
Filename        HPS108_combined.fastq.gz
Sequence length 100
Filename        HPS114_combined.fastq.gz
Sequence length 100
Filename        HPS116_combined.fastq.gz
Sequence length 100
Filename        HPS119_combined.fastq.gz
Sequence length 100
Filename        HPS121_combined.fastq.gz
Sequence length 100
Filename        HPS123_combined.fastq.gz
Sequence length 100
Filename        HPS124_combined.fastq.gz
Sequence length 100
Filename        HPS130_combined.fastq.gz
Sequence length 100
Filename        HPS139_combined.fastq.gz
Sequence length 100
Filename        HPS145_combined.fastq.gz
Sequence length 100
Filename        HPS146_combined.fastq.gz
Sequence length 100
Filename        HPS147_combined.fastq.gz
Sequence length 100
Filename        HPS150_combined.fastq.gz
Sequence length 100
Filename        HPS152_combined.fastq.gz
Sequence length 100
Filename        HPS155_combined.fastq.gz
Sequence length 100
Filename        HPS163_combined.fastq.gz
Sequence length 100
Filename        HPS172_combined.fastq.gz
Sequence length 100
Filename        HPS175_combined.fastq.gz
Sequence length 100
Filename        HPS183_combined.fastq.gz
Sequence length 100
Filename        HPS184_combined.fastq.gz
Sequence length 100
Filename        HPS186_combined.fastq.gz
Sequence length 100
Filename        HPS193_combined.fastq.gz
Sequence length 100
Filename        HPS199_combined.fastq.gz
Sequence length 100
Filename        HPS200_combined.fastq.gz
Sequence length 100
Filename        HPS203_combined.fastq.gz
Sequence length 100
Filename        HPS208_combined.fastq.gz
Sequence length 100
Filename        HPS210_combined.fastq.gz
Sequence length 100
Filename        HPS217_combined.fastq.gz
Sequence length 100
Filename        HPS223_combined.fastq.gz
Sequence length 100
Filename        HPS227_combined.fastq.gz
Sequence length 100
Filename        HPS230_combined.fastq.gz
Sequence length 100
Filename        HPS234_combined.fastq.gz
Sequence length 100
Filename        HPS240_combined.fastq.gz
Sequence length 100
Filename        HPS255_combined.fastq.gz
Sequence length 100
Filename        HPS257_combined.fastq.gz
Sequence length 100
Filename        HPS259_combined.fastq.gz
Sequence length 100
Filename        HPS263_combined.fastq.gz
Sequence length 100
Filename        HPS267_combined.fastq.gz
Sequence length 100
Filename        HPS269_combined.fastq.gz
Sequence length 100
Filename        HPS270_combined.fastq.gz
Sequence length 100
Filename        HPS281_combined.fastq.gz
Sequence length 100
Filename        HPS283_combined.fastq.gz
Sequence length 100
Filename        HPS286_combined.fastq.gz
Sequence length 100
Filename        HPS296_combined.fastq.gz
Sequence length 100
```

All sequences 100bp in length as expected of the data

The report at `./results/pipeline_info/results_description.html`, mentions:

```bash
The FastQC plots displayed in the MultiQC report shows untrimmed reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the trim_galore directory
```

Therefore I'll look at the qc of the reads after trimming in the fastqc reports in the trim_galore directory

#### i. Insert Size calculation

```bash
/store/lkemp/smrnaseq_hps/results/trim_galore/insertsize/
```

#### ii. Collapse reads seqcsluter

### 2. Adapter trimming - Trim Galore!

For each sample, get the number and percentage of sequences removed because they became shorter than the length cutoff of 17 bp (`--min-length 17` flag in pipeline)

```bash
for summaryfile in ./results/trim_galore/HPS*_combined.fastq.gz_trimming_report.txt; do
grep -e 'Sequences removed because they became shorter than the length cutoff of 17 bp:' $summaryfile
done
```

Output:

```bash
Sequences removed because they became shorter than the length cutoff of 17 bp:  25580 (0.4%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  29988 (0.5%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  29263 (0.5%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  23180 (0.4%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  28554 (0.5%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  14396 (0.2%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  16588 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  19147 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  24621 (0.4%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  27172 (0.4%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  26108 (0.4%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  21945 (0.4%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  35793 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  87048 (0.8%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  46670 (0.5%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  36173 (0.4%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  77529 (0.8%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  56023 (0.6%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  39229 (0.4%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  27187 (0.2%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  68674 (0.7%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  40526 (0.4%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  40470 (0.4%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  40896 (0.4%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  32301 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  26202 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  28379 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  27855 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  34126 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  28716 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  26134 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  33863 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  25767 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  24703 (0.2%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  42221 (0.5%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  31456 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  33895 (0.4%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  30859 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  37795 (0.4%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  35930 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  34672 (0.4%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  40896 (0.4%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  30316 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  34424 (0.4%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  30531 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  40893 (0.4%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  35933 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  30342 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  19830 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  19640 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  19601 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  18808 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  21287 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  17970 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  18628 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  21733 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  28688 (0.4%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  22001 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  25078 (0.3%)
Sequences removed because they became shorter than the length cutoff of 17 bp:  17632 (0.3%)
```

For each sample, get the number and percentage of sequences removed because they were longer than the maximum length cutoff of 40 bp

```bash
for summaryfile in ./results/trim_galore/HPS*_combined.fastq.gz_trimming_report.txt; do
grep -e 'Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:' $summaryfile
done
```

Output:

```bash
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      6261929 (99.5%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      5762353 (99.4%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      5979309 (99.5%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      5689589 (99.6%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      6232197 (99.5%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      6138339 (99.7%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      6287192 (99.7%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      6143959 (99.7%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      5900113 (99.5%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      6457960 (99.5%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      6296118 (99.5%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      6236627 (99.6%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      10775840 (99.6%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      10290202 (99.1%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      9808112 (99.5%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      9765126 (99.6%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      9841520 (99.1%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      9596440 (99.4%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      10393168 (99.6%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      11037573 (99.7%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      10421673 (99.3%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      10943833 (99.6%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      11439897 (99.6%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      11348876 (99.6%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      10595898 (99.7%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      9951570 (99.7%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      10702668 (99.7%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      10665847 (99.7%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      10188208 (99.6%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      11439895 (99.7%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      10112447 (99.7%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      11061912 (99.7%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      8562172 (99.7%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      9976086 (99.7%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      9159838 (99.5%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      8993455 (99.6%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      9580386 (99.6%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      10153611 (99.6%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      9733691 (99.6%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      11021045 (99.6%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      9125361 (99.6%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      9967236 (99.5%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      9369937 (99.6%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      9208237 (99.6%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      9737807 (99.6%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      9995091 (99.5%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      10627321 (99.6%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      10621532 (99.7%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      6594950 (99.7%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      7284780 (99.7%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      6630465 (99.7%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      6280895 (99.7%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      6281218 (99.6%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      6555044 (99.7%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      5887612 (99.6%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      6198497 (99.6%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      6465871 (99.5%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      6795455 (99.6%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      7237398 (99.6%)
Sequences removed because after trimming they were longer than the maximum length cutoff of 40 bp:      6316988 (99.7%)
```

Seems like the majority of sequences are larger than 40bp in length and are being removed from the analysis. Before trimming, our sequences were 100bp in length. Is this what we expect? I haven't set the maximum length cutoff in the pipeline run, but it could be set somewhere internally/automatically in the pipeline.

Get sequence lengths for all reads *after* trimming

```bash
mkdir ./results/fastqc_trimmed_unzipped

# Unzip fastqc post trimming files 
for zip in ./results/trim_galore/*_combined_trimmed_fastqc.zip; do
unzip $zip -d ./results/fastqc_trimmed_unzipped/
done

# Get sequence lengths after trimming for all samples from fastqc summary files
for summaryfile in ./results/fastqc_trimmed_unzipped/HPS*_combined_trimmed_fastqc/fastqc_data.txt; do
grep -e 'Filename' -e 'Sequence length' $summaryfile
done
```

Output:

```bash
Filename        HPS004_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS022_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS030_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS033_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS041_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS048_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS061_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS064_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS070_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS071_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS072_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS080_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS081_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS083_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS086_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS100_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS108_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS114_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS116_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS119_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS121_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS123_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS124_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS130_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS139_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS145_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS146_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS147_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS150_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS152_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS155_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS163_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS172_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS175_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS183_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS184_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS186_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS193_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS199_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS200_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS203_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS208_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS210_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS217_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS223_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS227_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS230_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS234_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS240_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS255_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS257_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS259_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS263_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS267_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS269_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS270_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS281_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS283_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS286_combined_trimmed.fq.gz
Sequence length 17-40
Filename        HPS296_combined_trimmed.fq.gz
Sequence length 17-40
```

Now all the sequences are between the 17bp and 40bp length cutoffs

There is a `--max_length` flag for TrimGalore that I could play around with, but I'll have to see if the pipeline take this value.


## Re-run data through pipeline

It turns out there is a bug in the trimming in this pipeline that Miles fixed earlier this year (see the github issue [here](https://github.com/nf-core/smrnaseq/issues/41)). I hadn't realised that the current latest release of [nf-core/smrnaseq](https://github.com/nf-core/smrnaseq) is [v1.0.0](https://github.com/nf-core/smrnaseq/releases/tag/1.0.0) which I would have used above doesn't include this fix. I'll re-run the pipeline with a newer version of the pipeline (the [dev branch](https://github.com/nf-core/smrnaseq/tree/dev) that includes this fix).

```bash
cd /store/lkemp/smrnaseq_hps/

# Move old analysis to folder
mkdir master_branch_analysis
mv .nextflow.log ./master_branch_analysis
mv ./work/ ./master_branch_analysis
mv ./results/ ./master_branch_analysis
mv ./smrnaseq ./master_branch_analysis
mv .nextflow ./master_branch_analysis

# Setup to run pipeline on dev branch of pipeline
mkdir dev_branch_analysis
cd dev_branch_analysis

git clone https://github.com/nf-core/smrnaseq.git
cd smrnaseq
git checkout dev
cd ..

conda activate nextflow

nextflow run /store/lkemp/smrnaseq_hps/dev_branch_analysis/smrnaseq/main.nf --input '/store/lkemp/smrnaseq_hps/fastq/*.fastq.gz' -profile conda --protocol illumina --mature --hairpin
```

Getting the following error only when running the dev branch (works fine when running on the master branch):

```bash
Mature file not found: false
```

Error coming from https://github.com/nf-core/smrnaseq/blob/a115116a06ad274b475834316d0193e1b756b249/main.nf#L130

It looks as though the dev branch (in contrast to the master branch) is expecting you to pass these variables:

```bash
      --mature [file]               Path to the FASTA file of mature miRNAs
      --hairpin [file]              Path to the FASTA file of miRNA precursors
```

Get these files from [mirbase](http://mirbase.org/)

```bash
wget ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz
wget ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz
```

Run again

```bash
nextflow run /store/lkemp/smrnaseq_hps/dev_branch_analysis/smrnaseq/main.nf \
--input '/store/lkemp/smrnaseq_hps/fastq/*.fastq.gz' \
-profile conda \
--protocol illumina \
--mature /store/lkemp/smrnaseq_hps/dev_branch_analysis/mature.fa.gz \
--hairpin /store/lkemp/smrnaseq_hps/dev_branch_analysis/hairpin.fa.gz
```

Error

```bash
Reference genome file not found: false
```

Looks like it wants a reference genome file passed to this flag:

```bash
      --fasta [file]                Path to fasta reference
```

Before getting the reference genome file myself, I'll see if it's happy if I just specify the reference genome using the `--genome` flag since it looks like it might be able to grab it from the iGenomes reference? (will use GRCh38 since this is what `mature.fa.gz` and `hairpin.fa.gz` are based on)

```bash
nextflow run /store/lkemp/smrnaseq_hps/dev_branch_analysis/smrnaseq/main.nf \
--input '/store/lkemp/smrnaseq_hps/fastq/*.fastq.gz' \
-profile conda \
--protocol illumina \
--mature /store/lkemp/smrnaseq_hps/dev_branch_analysis/mature.fa.gz \
--hairpin /store/lkemp/smrnaseq_hps/dev_branch_analysis/hairpin.fa.gz \
--genome GRCh38
```

Error

```bash
Reference species for miRTrace is not defined.
```

Try specifying reference species for miRTrace with the `--mirtrace_species` flag

```bash
nextflow run /store/lkemp/smrnaseq_hps/dev_branch_analysis/smrnaseq/main.nf \
--input '/store/lkemp/smrnaseq_hps/fastq/*.fastq.gz' \
-profile conda \
--protocol illumina \
--mature /store/lkemp/smrnaseq_hps/dev_branch_analysis/mature.fa.gz \
--hairpin /store/lkemp/smrnaseq_hps/dev_branch_analysis/hairpin.fa.gz \
--genome GRCh38 \
--mirtrace_species hsa
```

It worked!! (launched the pipeline)

Although I'm still getting some warnings:

```bash
WARN: Access to undefined parameter `mirna_gtf` -- Initialise it to a default value eg. `params.mirna_gtf = some_value`
WARN: Access to undefined parameter `bt_indices` -- Initialise it to a default value eg. `params.bt_indices = some_value`
No GTF / Bowtie 1 index supplied - host reference genome analysis will be skipped.
WARN: Access to undefined parameter `input_paths` -- Initialise it to a default value eg. `params.input_paths = some_value`
```

Get GFF/GTF file with coordinates positions of precursor and miRNAs

```bash
wget ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3
```

The warning about the bt_indices looks like they want us to provide an index for the reference genome in order to run the host reference genome analysis (even though it's optional to pass this file). See the `main.nf` file:

```bash
      --bt_index [file]             Path to the bowtie 1 index files of the host reference genome. Optional.
```

The [bowtie manual](http://bowtie-bio.sourceforge.net/manual.shtml) indicates where to get these index files (they indicate to get it from NCBI)

```bash
# Get the index files for bowtie
wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_decoy_as.zip
# Unzip
unzip GRCh38_noalt_decoy_as.zip
```

Run again

```bash
nextflow run /store/lkemp/smrnaseq_hps/dev_branch_analysis/smrnaseq/main.nf \
--input '/store/lkemp/smrnaseq_hps/fastq/*_combined.fastq.gz' \
-profile conda \
--protocol illumina \
--mature /store/lkemp/smrnaseq_hps/dev_branch_analysis/mature.fa.gz \
--hairpin /store/lkemp/smrnaseq_hps/dev_branch_analysis/hairpin.fa.gz \
--genome GRCh38 \
--mirtrace_species hsa \
--mirna_gtf /store/lkemp/smrnaseq_hps/dev_branch_analysis/hsa.gff3 \
--bt_index /store/lkemp/smrnaseq_hps/dev_branch_analysis/GRCh38_noalt_decoy_as/ \
```

Just one remaining warning left: 

```bash
WARN: Access to undefined parameter `input_paths` -- Initialise it to a default value eg. `params.input_paths = some_value`
```

I won't worry about that for now, the input fastq files seem to be taken by the pipeline just fine. Now I want to add back some of the trimming and other parameters as well as switch from using conda to singularity (more robust and recommended by the software developers in the README to use conda as a last resort). I also want to reduce some of the resources the pipeline uses - I'll reduce the max_cpus variable in the pipeline config file (`/store/lkemp/smrnaseq_hps/dev_branch_analysis/smrnaseq/nextflow.config`) from 16 to 10. Run the pipeline:

```bash
nextflow run /store/lkemp/smrnaseq_hps/dev_branch_analysis/smrnaseq/main.nf \
--input '/store/lkemp/smrnaseq_hps/fastq/*_combined.fastq.gz' \
-profile singularity \
--protocol illumina \
--mature /store/lkemp/smrnaseq_hps/dev_branch_analysis/mature.fa.gz \
--hairpin /store/lkemp/smrnaseq_hps/dev_branch_analysis/hairpin.fa.gz \
--genome GRCh38 \
--mirtrace_species hsa \
--mirna_gtf /store/lkemp/smrnaseq_hps/dev_branch_analysis/hsa.gff3 \
--bt_index /store/lkemp/smrnaseq_hps/dev_branch_analysis/GRCh38_noalt_decoy_as/ \
--saveReference \
-resume \
--three_prime_adapter AGATCGGAAGAGCACACG \
--min_length 17 \
--email leah.kemp@esr.cri.nz \
--email_on_fail leah.kemp@esr.cri.nz
```

Failed with error:

```bash
Error executing process > 'bowtie_ref (HPS061_combined_trimmed.fq.gz)'

Caused by:
  Process `bowtie_ref (HPS061_combined_trimmed.fq.gz)` terminated with an error exit status (1)

Command executed:

  bowtie \
      GRCh38_noalt_decoy_as.4.bt2 \
      -q <(zcat HPS061_combined_trimmed.fq.gz) \
      -p 10 \
      -t \
      -k 50 \
      --best \
      --strata \
      -e 99999 \
      --chunkmbs 2048 \
      -S  \
      | samtools view -bS - > HPS061_combined.genome.bam

Command exit status:
  1

Command output:
  (empty)

Command error:
  Could not locate a Bowtie index corresponding to basename "GRCh38_noalt_decoy_as.4.bt2"
  Overall time: 00:00:00
  Command: /opt/conda/envs/nf-core-smrnaseq-1.1dev/bin/bowtie-align-s --wrapper basic-0 -q -p 10 -t -k 50 --best --strata -e 99999 --chunkmbs 2048 -S GRCh38_noalt_decoy_as.4.bt2 /dev/fd/63
```

Will try run passing the `GRCh38_noalt_decoy_as.4.bt2` file directly to the `--bt-index` flag

```bash
nextflow run /store/lkemp/smrnaseq_hps/dev_branch_analysis/smrnaseq/main.nf \
--input '/store/lkemp/smrnaseq_hps/fastq/*_combined.fastq.gz' \
-profile singularity \
--protocol illumina \
--mature /store/lkemp/smrnaseq_hps/dev_branch_analysis/mature.fa.gz \
--hairpin /store/lkemp/smrnaseq_hps/dev_branch_analysis/hairpin.fa.gz \
--genome GRCh38 \
--mirtrace_species hsa \
--mirna_gtf /store/lkemp/smrnaseq_hps/dev_branch_analysis/hsa.gff3 \
--bt_index /store/lkemp/smrnaseq_hps/dev_branch_analysis/GRCh38_noalt_decoy_as/GRCh38_noalt_decoy_as.4.bt2 \
--saveReference \
-resume \
--three_prime_adapter AGATCGGAAGAGCACACG \
--min_length 17 \
--email leah.kemp@esr.cri.nz \
--email_on_fail leah.kemp@esr.cri.nz
```

Error thrown again. Try other bowtie index file

```bash
wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
unzip GRCh38_noalt_as.zip
```

```bash
nextflow run /store/lkemp/smrnaseq_hps/dev_branch_analysis/smrnaseq/main.nf \
--input '/store/lkemp/smrnaseq_hps/fastq/*_combined.fastq.gz' \
-profile singularity \
--protocol illumina \
--mature /store/lkemp/smrnaseq_hps/dev_branch_analysis/mature.fa.gz \
--hairpin /store/lkemp/smrnaseq_hps/dev_branch_analysis/hairpin.fa.gz \
--genome GRCh38 \
--mirtrace_species hsa \
--mirna_gtf /store/lkemp/smrnaseq_hps/dev_branch_analysis/hsa.gff3 \
--bt_index /store/lkemp/smrnaseq_hps/dev_branch_analysis/GRCh38_noalt_as/ \
--saveReference \
-resume \
--three_prime_adapter AGATCGGAAGAGCACACG \
--min_length 17 \
--email leah.kemp@esr.cri.nz \
--email_on_fail leah.kemp@esr.cri.nz
```

Didn't work

From [this question on stackoverflow](https://stackoverflow.com/questions/38502194/could-not-locate-a-bowtie-index-corresponding-to-basename/38568771) it looks like I was passing the indexes incorrectly

```bash
nextflow run /store/lkemp/smrnaseq_hps/dev_branch_analysis/smrnaseq/main.nf \
--input '/store/lkemp/smrnaseq_hps/fastq/*_combined.fastq.gz' \
-profile singularity \
--protocol illumina \
--mature /store/lkemp/smrnaseq_hps/dev_branch_analysis/mature.fa.gz \
--hairpin /store/lkemp/smrnaseq_hps/dev_branch_analysis/hairpin.fa.gz \
--genome GRCh38 \
--mirtrace_species hsa \
--mirna_gtf /store/lkemp/smrnaseq_hps/dev_branch_analysis/hsa.gff3 \
--bt_index /store/lkemp/smrnaseq_hps/dev_branch_analysis/GRCh38_noalt_as/GRCh38_noalt_as \
--saveReference \
-resume \
--three_prime_adapter AGATCGGAAGAGCACACG \
--min_length 17 \
--email leah.kemp@esr.cri.nz \
--email_on_fail leah.kemp@esr.cri.nz
```

Didn't work

Try on GRCh37

```bash
wget https://genome-idx.s3.amazonaws.com/bt/GRCh37.zip
unzip GRCh37.zip
wget ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3
```

```bash
nextflow run /store/lkemp/smrnaseq_hps/dev_branch_analysis/smrnaseq/main.nf \
--input '/store/lkemp/smrnaseq_hps/fastq/*_combined.fastq.gz' \
-profile singularity \
--protocol illumina \
--mature /store/lkemp/smrnaseq_hps/dev_branch_analysis/mature.fa.gz \
--hairpin /store/lkemp/smrnaseq_hps/dev_branch_analysis/hairpin.fa.gz \
--genome GRCh38 \
--mirtrace_species hsa \
--mirna_gtf /store/lkemp/smrnaseq_hps/dev_branch_analysis/hsa.gff3 \
--bt_index /store/lkemp/smrnaseq_hps/dev_branch_analysis/GRCh37/GRCh37 \
--saveReference \
-resume \
--three_prime_adapter AGATCGGAAGAGCACACG \
--min_length 17 \
--email leah.kemp@esr.cri.nz \
--email_on_fail leah.kemp@esr.cri.nz
```

Didn't work, try with index files from the ftp site (ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/)

```bash
mkdir bowtie_indexes
cd bowtie_indexes
wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/GRCh38.*.ebwt
cd ..
```

```bash
nextflow run /store/lkemp/smrnaseq_hps/dev_branch_analysis/smrnaseq/main.nf \
--input '/store/lkemp/smrnaseq_hps/fastq/*_combined.fastq.gz' \
-profile singularity \
--protocol illumina \
--mature /store/lkemp/smrnaseq_hps/dev_branch_analysis/mature.fa.gz \
--hairpin /store/lkemp/smrnaseq_hps/dev_branch_analysis/hairpin.fa.gz \
--genome GRCh38 \
--mirtrace_species hsa \
--mirna_gtf /store/lkemp/smrnaseq_hps/dev_branch_analysis/hsa.gff3 \
--bt_index /store/lkemp/smrnaseq_hps/dev_branch_analysis/GRCh38/GRCh38 \
--saveReference \
-resume \
--three_prime_adapter AGATCGGAAGAGCACACG \
--min_length 17 \
--email leah.kemp@esr.cri.nz \
--email_on_fail leah.kemp@esr.cri.nz
```

Run without bowtie index (and therefore without the host reference genome analysis)

```bash
nextflow run /store/lkemp/smrnaseq_hps/dev_branch_analysis/smrnaseq/main.nf \
--input '/store/lkemp/smrnaseq_hps/fastq/*_combined.fastq.gz' \
-profile singularity \
--protocol illumina \
--mature /store/lkemp/smrnaseq_hps/dev_branch_analysis/mature.fa.gz \
--hairpin /store/lkemp/smrnaseq_hps/dev_branch_analysis/hairpin.fa.gz \
--genome GRCh38 \
--mirtrace_species hsa \
--mirna_gtf /store/lkemp/smrnaseq_hps/dev_branch_analysis/hsa.gff3 \
--saveReference \
-resume \
--three_prime_adapter AGATCGGAAGAGCACACG \
--min_length 17 \
--email leah.kemp@esr.cri.nz \
--email_on_fail leah.kemp@esr.cri.nz
```

Output:

```bash

```

## Explore the outputs!

### 1. Raw read QC - FastQC

Get sequence lengths for all reads *before* trimming

```bash
mkdir ./results/fastqc_unzipped/

# Unzip fastqc files
for zip in ./results/fastqc/*_combined_fastqc.zip; do
unzip $zip -d ./results/fastqc_unzipped/
done

# Get sequence lengths before trimming for all samples from fastqc summary files
for summaryfile in ./results/fastqc_unzipped/HPS*_combined_fastqc/fastqc_data.txt; do
grep -e 'Filename' -e 'Sequence length' $summaryfile
done
```

Output:

```bash

```

#### i. Insert Size calculation

```bash
/store/lkemp/smrnaseq_hps/results/trim_galore/insertsize/
```

#### ii. Collapse reads seqcsluter

### 2. Adapter trimming - Trim Galore!

### 3. Alignment against miRBase mature miRNA - Bowtie1

### 4. Alignment against miRBase hairpin

#### Unaligned reads from step 3 (Bowtie1)

#### Collapsed reads from step 2.2 (Bowtie1)

### 5. Post-alignment processing of miRBase hairpin

#### i. Basic statistics from step 3 and step 4.1 (SAMtools)

#### ii. Analysis on miRBase hairpin counts (edgeR)

TMM normalization and a table of top expression hairpin

MDS plot clustering samples

Heatmap of sample similarities

#### iii. miRNA and isomiR annotation from step 4.1 (mirtop)

### 6. Alignment against host reference genome (Bowtie1)

#### i. Post-alignment processing of alignment against host reference genome (SAMtools)

### 7. miRNA quality control (mirtrace)

### 8. Present QC for raw read, alignment, and expression results (MultiQC)

