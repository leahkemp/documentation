# Benchmarking genomic pipelines - resources

Created: 2020-04-22 13:37:04
Last modified: 2020/05/25 16:39:30

- **Aim:** Undertake benchmarking of genomics pipelines to optimise their resource use.
- **Prerequisite software:** [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html)
- **OS:** Ubuntu 16.04 (Wintermute - research server)

The idea is to run these pipelines against a reduced Genome In A Bottle (GIAB) sample [NA12878](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) and evaluate the real time and user time for each step in the pipelines.

## Table of contents

- [Benchmarking genomic pipelines - resources](#benchmarking-genomic-pipelines---resources)
  - [Table of contents](#table-of-contents)
  - [Create a test dataset](#create-a-test-dataset)
    - [Download and prepare the test WES data](#download-and-prepare-the-test-wes-data)
    - [Reduce the test dataset size for benchmarking](#reduce-the-test-dataset-size-for-benchmarking)
      - [Run the two samples (NIST7035 and NIST7086) through the bwa_map step of human_genomics_pipeline](#run-the-two-samples-nist7035-and-nist7086-through-the-bwamap-step-of-humangenomicspipeline)
  - [Testing](#testing)
    - [Setup](#setup)
    - [Wintermute](#wintermute)
      - [human_genomics_pipeline](#humangenomicspipeline)
    - [vcf_annotation_pipeline](#vcfannotationpipeline)
    - [Production](#production)
      - [human_genomics_pipeline](#humangenomicspipeline-1)
    - [vcf_annotation_pipeline](#vcfannotationpipeline-1)

## Create a test dataset

### Download and prepare the test WES data

Download

```bash
cd /store/lkemp/publicData/exomes/NA12878_exome/
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/Garvan_NA12878_HG001_HiSeq_Exome.README
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L002_R1_001.fastq.gz
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L002_R2_001.fastq.gz
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7086_CGTACTAG_L001_R1_001.fastq.gz
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7086_CGTACTAG_L001_R2_001.fastq.gz
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7086_CGTACTAG_L002_R1_001.fastq.gz
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7086_CGTACTAG_L002_R2_001.fastq.gz
```

Note: there can be some issues with correctly downloading the fastq files (likely due to the ESR proxy). The below can be run to check the integrity of the gunzip files, no error suggests the archive is OK.

```bash
for i in NIST*; do
  echo "$i" ;
  gunzip -t "$i"
  echo "...done..."
done
```

Collapse pooled runs

```bash
# NIST7035
cat NIST7035*_R1_001.fastq.gz > NIST7035_NIST_R1.fastq.gz
cat NIST7035*_R2_001.fastq.gz > NIST7035_NIST_R2.fastq.gz

# NIST7086
cat NIST7086*_R1_001.fastq.gz > NIST7086_NIST_R1.fastq.gz
cat NIST7086*_R2_001.fastq.gz > NIST7086_NIST_R2.fastq.gz
```

### Reduce the test dataset size for benchmarking

#### Run the two samples (NIST7035 and NIST7086) through the bwa_map step of human_genomics_pipeline

Create a local copy of NIST7035 and NIST7086 fastq files to run through pipeline (and manipulate later)

```bash
cd /store/lkemp/exome_project/resource_benchmarking/

mkdir fastq
cd fastq/

cp /store/lkemp/publicData/exomes/NA12878_exome/NIST7035_NIST_R1.fastq.gz .
cp /store/lkemp/publicData/exomes/NA12878_exome/NIST7035_NIST_R2.fastq.gz .
cp /store/lkemp/publicData/exomes/NA12878_exome/NIST7086_NIST_R1.fastq.gz .
cp /store/lkemp/publicData/exomes/NA12878_exome/NIST7086_NIST_R2.fastq.gz .
```

Clone [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline)

```bash
cd ..
git clone git@github.com:ESR-NZ/human_genomics_pipeline.git
cd human_genomics_pipeline
```

Create and activate a conda environment with pytho, snakemake and samtools installed

```bash
conda create -n pipeline_env python=3.7
conda activate pipeline_env
conda install -c bioconda snakemake=5.14.0
conda install -c bioconda samtools=1.10
```

Run the pipeline up to and including the gatk_apply_bqsr rule to get a full bam file mapped to the reference genome

```bash
# Dry run
snakemake -n -j 24 --use-conda --configfile config/resource_benchmarking.yml --until gatk_apply_bqsr

# Full run
snakemake -j 24 --use-conda --configfile config/resource_benchmarking.yml --until gatk_apply_bqsr
```

Extract chr1 from the bam files (retain paired reads)

```bash
cd /store/lkemp/exome_project/resource_benchmarking/human_genomics_pipeline/mapped/

# NIST7035
samtools view NIST7035_NIST_bwa_recal.bam chr1 -b > NIST7035_NIST_bwa_recal_chr1.bam

# NIST7086
samtools view NIST7086_NIST_bwa_recal.bam chr1 -b > NIST7086_NIST_bwa_recal_chr1.bam
```

Extract the fastq reads (R1 and R2) from the raw fastq files that correspond to chr1 (the reads that map to chr1 in the bam files)

```bash
# NIST7035
samtools bam2fq NIST7035_NIST_bwa_recal_chr1.bam > NIST7035_NIST_bwa_recal_chr1.fastq
cat NIST7035_NIST_bwa_recal_chr1.fastq | grep '^@.*/1$' -A 3 --no-group-separator > NIST7035_NIST_chr1_R1.fastq
cat NIST7035_NIST_bwa_recal_chr1.fastq | grep '^@.*/2$' -A 3 --no-group-separator > NIST7035_NIST_chr1_R2.fastq

# NIST7086
samtools bam2fq NIST7086_NIST_bwa_recal_chr1.bam > NIST7086_NIST_bwa_recal_chr1.fastq
cat NIST7086_NIST_bwa_recal_chr1.fastq | grep '^@.*/1$' -A 3 --no-group-separator > NIST7086_NIST_chr1_R1.fastq
cat NIST7086_NIST_bwa_recal_chr1.fastq | grep '^@.*/2$' -A 3 --no-group-separator > NIST7086_NIST_chr1_R2.fastq
```

Zip fastq files

```bash
bgzip NIST7035_NIST_chr1_R1.fastq
bgzip NIST7035_NIST_chr1_R2.fastq
bgzip NIST7086_NIST_chr1_R1.fastq
bgzip NIST7086_NIST_chr1_R2.fastq
```

This will reduce the size of the dataset (width) without reducing the depth of reads that could influence the performance of the downstream rules (which could skew the resource benchmarking results for these steps).

## Testing

Approach:

- Run each rule separately
- Run with doubling threads: 1, 2, 4, 8, 16 etc. (until runtime plateaus)
- Compare real time vs. user time (minimise any divergence between them)
- Re-run benchmarking on new machines (eg. production) to fine-tune resource allocation

### Setup

Replace the full fastq input file with the reduced ones

```bash
cd /store/lkemp/exome_project/resource_benchmarking/fastq/
rm -r NIST70*
cp ../human_genomics_pipeline/mapped/NIST7035_NIST_chr1* .
cp ../human_genomics_pipeline/mapped/NIST7086_NIST_chr1* .
```

Clone a fresh pipeline and checkout the branch for resource benchmarking on Wintermute

```bash
cd ..
rm -r human_genomic_pipeline
git clone git@github.com:ESR-NZ/human_genomics_pipeline.git
cd human_genomics_pipeline
git checkout resource_benchmarking_wintermute
```

### Wintermute

#### human_genomics_pipeline

Set the number of threads in each rule to the maximum number that we will test (32), therefore we can control the number of threads used for each test with the `j` flag passed to snakemake on the command line (if the number of threads for a given rule are larger that the threads passed to this flag, they will be scaled down)

Run rules in order shown

*note. I removed all files created by a rule before re-running the same rule with a new number of threads. I also made sure a rules conda environment was created prior to running benhcmarking for a rule to ensure the time represents on;y activatea conda environment and not creating one from scratch. One sample was used NIST7035_NIST_chr1*

- fastqc

```bash
time snakemake -j 1 --use-conda --configfile config/resource_benchmarking.yml --until fastqc
time snakemake -j 2 --use-conda --configfile config/resource_benchmarking.yml --until fastqc
time snakemake -j 4 --use-conda --configfile config/resource_benchmarking.yml --until fastqc
time snakemake -j 8 --use-conda --configfile config/resource_benchmarking.yml --until fastqc
time snakemake -j 16 --use-conda --configfile config/resource_benchmarking.yml --until fastqc
time snakemake -j 32 --use-conda --configfile config/resource_benchmarking.yml --until fastqc
```

- multiqc

```bash
time snakemake -j 1 --use-conda --configfile config/resource_benchmarking.yml --until multiqc_pre_trim
time snakemake -j 2 --use-conda --configfile config/resource_benchmarking.yml --until multiqc_pre_trim
time snakemake -j 4 --use-conda --configfile config/resource_benchmarking.yml --until multiqc_pre_trim
time snakemake -j 8 --use-conda --configfile config/resource_benchmarking.yml --until multiqc_pre_trim
time snakemake -j 16 --use-conda --configfile config/resource_benchmarking.yml --until multiqc_pre_trim
time snakemake -j 32 --use-conda --configfile config/resource_benchmarking.yml --until multiqc_pre_trim
```

- trim_galore

```bash
time snakemake -j 1 --use-conda --configfile config/resource_benchmarking.yml --until trim_galore_pe
time snakemake -j 2 --use-conda --configfile config/resource_benchmarking.yml --until trim_galore_pe
time snakemake -j 4 --use-conda --configfile config/resource_benchmarking.yml --until trim_galore_pe
time snakemake -j 8 --use-conda --configfile config/resource_benchmarking.yml --until trim_galore_pe
time snakemake -j 16 --use-conda --configfile config/resource_benchmarking.yml --until trim_galore_pe
time snakemake -j 32 --use-conda --configfile config/resource_benchmarking.yml --until trim_galore_pe
```

- multiqc_post_trim

```bash
time snakemake -j 1 --use-conda --configfile config/resource_benchmarking.yml --until multiqc_post_trim
time snakemake -j 2 --use-conda --configfile config/resource_benchmarking.yml --until multiqc_post_trim
time snakemake -j 4 --use-conda --configfile config/resource_benchmarking.yml --until multiqc_post_trim
time snakemake -j 8 --use-conda --configfile config/resource_benchmarking.yml --until multiqc_post_trim
time snakemake -j 16 --use-conda --configfile config/resource_benchmarking.yml --until multiqc_post_trim
time snakemake -j 32 --use-conda --configfile config/resource_benchmarking.yml --until multiqc_post_trim
```

- bwa_map

```bash
time snakemake -j 1 --use-conda --configfile config/resource_benchmarking.yml --until bwa_map
time snakemake -j 2 --use-conda --configfile config/resource_benchmarking.yml --until bwa_map
time snakemake -j 4 --use-conda --configfile config/resource_benchmarking.yml --until bwa_map
time snakemake -j 8 --use-conda --configfile config/resource_benchmarking.yml --until bwa_map
time snakemake -j 16 --use-conda --configfile config/resource_benchmarking.yml --until bwa_map
time snakemake -j 32 --use-conda --configfile config/resource_benchmarking.yml --until bwa_map
```

- sambamba_sort

```bash
time snakemake -j 1 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_sort
time snakemake -j 2 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_sort
time snakemake -j 4 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_sort
time snakemake -j 8 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_sort
time snakemake -j 16 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_sort
time snakemake -j 32 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_sort
```

- sambamba_mkdups

```bash
time snakemake -j 1 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_mkdups
time snakemake -j 2 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_mkdups
time snakemake -j 4 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_mkdups
time snakemake -j 8 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_mkdups
time snakemake -j 16 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_mkdups
time snakemake -j 32 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_mkdups
```

- sambamba_index

```bash
time snakemake -j 1 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_index
time snakemake -j 2 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_index
time snakemake -j 4 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_index
time snakemake -j 8 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_index
time snakemake -j 16 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_index
time snakemake -j 32 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_index
```

- gatk_add_replace_read_groups

```bash
time snakemake -j 1 --use-conda --configfile config/resource_benchmarking.yml --until gatk_add_replace_read_groups
time snakemake -j 2 --use-conda --configfile config/resource_benchmarking.yml --until gatk_add_replace_read_groups
time snakemake -j 4 --use-conda --configfile config/resource_benchmarking.yml --until gatk_add_replace_read_groups
time snakemake -j 8 --use-conda --configfile config/resource_benchmarking.yml --until gatk_add_replace_read_groups
time snakemake -j 16 --use-conda --configfile config/resource_benchmarking.yml --until gatk_add_replace_read_groups
time snakemake -j 32 --use-conda --configfile config/resource_benchmarking.yml --until gatk_add_replace_read_groups
```

- sambamba_index_rgadd

```bash
time snakemake -j 1 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_index_rgadd
time snakemake -j 2 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_index_rgadd
time snakemake -j 4 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_index_rgadd
time snakemake -j 8 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_index_rgadd
time snakemake -j 16 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_index_rgadd
time snakemake -j 32 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_index_rgadd
```

- gatk_base_recalibrator

```bash
time snakemake -j 1 --use-conda --configfile config/resource_benchmarking.yml --until gatk_base_recalibrator
time snakemake -j 2 --use-conda --configfile config/resource_benchmarking.yml --until gatk_base_recalibrator
time snakemake -j 4 --use-conda --configfile config/resource_benchmarking.yml --until gatk_base_recalibrator
time snakemake -j 8 --use-conda --configfile config/resource_benchmarking.yml --until gatk_base_recalibrator
time snakemake -j 16 --use-conda --configfile config/resource_benchmarking.yml --until gatk_base_recalibrator
time snakemake -j 32 --use-conda --configfile config/resource_benchmarking.yml --until gatk_base_recalibrator
```

- gatk_apply_qgsr

- gatk_haplotype_caller_single

```bash
time snakemake -j 24 --use-conda --configfile config/resource_benchmarking.yml --until bwa_map
```

### vcf_annotation_pipeline

resource_benchmarking_wintermute branch?

### Production

#### human_genomics_pipeline

### vcf_annotation_pipeline