# Benchmarking genomic pipelines - quality

Created: 2020-04-22 13:37:04
Last modified: 2020/09/15 13:30:25

- **Aim:** Undertake best practice benchmarking of genomics pipelines to test their quality for clinical use.
- **Prerequisite software:** [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html), [bgzip](http://www.htslib.org/doc/bgzip.html), [tabix](http://www.htslib.org/doc/tabix.html)
- **OS:** Ubuntu 16.04 (Wintermute - research server) for human_genomics_pipeline + vcf_annotation_pipeline runs and CentOS-7 (ORAC) for parabricks germline pipeline runs

The idea is to run these pipelines against the Genome In A Bottle (GIAB) sample [NA12878](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) and compare the quality of the variant calls to the truth vcf. Benchmarking will follow the best practices described in [this paper](https://www.nature.com/articles/s41587-019-0054-x).

This document aims to document the code used in each quality benchmarking run, whereas [this document](https://github.com/leahkemp/documentation/blob/master/benchmarking_pipelines_quality/benchmarking_pipelines_quality_results.md) aims to describe a few key run parameters and record the results of each quality benchmarking run. The full parameters for each benchmarking run can be found [here](https://github.com/ESR-NZ/human_genomics_pipeline/tree/quality_benchmarking), [here](https://github.com/ESR-NZ/vcf_annotation_pipeline/tree/quality_benchmarking) and [here](https://github.com/ESR-NZ/ESR-Parabricks/quality_benchmarking).

## Table of contents

- [Benchmarking genomic pipelines - quality](#benchmarking-genomic-pipelines---quality)
  - [Table of contents](#table-of-contents)
  - [Setup](#setup)
    - [Install benchmarking software](#install-benchmarking-software)
    - [Download and prepare data](#download-and-prepare-data)
      - [Fastq data to run through pipelines](#fastq-data-to-run-through-pipelines)
      - [Truth vcf to compare output to](#truth-vcf-to-compare-output-to)
    - [Other setup](#other-setup)
      - [sdf files for hap.py](#sdf-files-for-happy)
  - [Benchmarking](#benchmarking)
    - [intra_truth_comparison](#intra_truth_comparison)
    - [quality_bench1.0](#quality_bench10)
      - [human_genomics_pipeline + minimal vcf_annotation_pipeline](#human_genomics_pipeline--minimal-vcf_annotation_pipeline)
      - [parabricks germline pipeline](#parabricks-germline-pipeline)
    - [quality_bench1.1](#quality_bench11)
      - [human_genomics_pipeline + minimal vcf_annotation_pipeline](#human_genomics_pipeline--minimal-vcf_annotation_pipeline-1)
    - [quality_bench1.2](#quality_bench12)
      - [human_genomics_pipeline + minimal vcf_annotation_pipeline](#human_genomics_pipeline--minimal-vcf_annotation_pipeline-2)
    - [quality_bench1.3](#quality_bench13)
      - [human_genomics_pipeline + minimal vcf_annotation_pipeline](#human_genomics_pipeline--minimal-vcf_annotation_pipeline-3)
    - [quality_bench1.4](#quality_bench14)
      - [human_genomics_pipeline + minimal vcf_annotation_pipeline](#human_genomics_pipeline--minimal-vcf_annotation_pipeline-4)
      - [parabricks germline pipeline](#parabricks-germline-pipeline-1)
    - [quality_bench1.5](#quality_bench15)
      - [no pipeline (intra_truth_comparison re-run)](#no-pipeline-intra_truth_comparison-re-run)
      - [human_genomics_pipeline + minimal vcf_annotation_pipeline (quality_bench1.0 re-run)](#human_genomics_pipeline--minimal-vcf_annotation_pipeline-quality_bench10-re-run)
      - [parabricks germline pipeline (quality_bench1.0 re-run)](#parabricks-germline-pipeline-quality_bench10-re-run)
      - [human_genomics_pipeline + minimal vcf_annotation_pipeline (quality_bench1.1 re-run)](#human_genomics_pipeline--minimal-vcf_annotation_pipeline-quality_bench11-re-run)
      - [human_genomics_pipeline + minimal vcf_annotation_pipeline (quality_bench1.2 re-run)](#human_genomics_pipeline--minimal-vcf_annotation_pipeline-quality_bench12-re-run)
      - [human_genomics_pipeline + minimal vcf_annotation_pipeline (quality_bench1.3 re-run)](#human_genomics_pipeline--minimal-vcf_annotation_pipeline-quality_bench13-re-run)
      - [human_genomics_pipeline + minimal vcf_annotation_pipeline (quality_bench1.4 re-run)](#human_genomics_pipeline--minimal-vcf_annotation_pipeline-quality_bench14-re-run)
      - [parabricks germline pipeline (quality_bench1.4 re-run)](#parabricks-germline-pipeline-quality_bench14-re-run)
    - [quality_bench1.6](#quality_bench16)
      - [human_genomics_pipeline + minimal vcf_annotation_pipeline](#human_genomics_pipeline--minimal-vcf_annotation_pipeline-5)
    - [quality_bench1.7](#quality_bench17)
      - [human_genomics_pipeline + minimal vcf_annotation_pipeline](#human_genomics_pipeline--minimal-vcf_annotation_pipeline-6)

## Setup

### Install benchmarking software

Create a conda environment for installing hap.py (we need python 2 in order to get htslib2 that is required for installation of hap.py)

```bash
conda create -n happy_install_env python=2.7.15
conda activate happy_install_env
```

Install dependencies for [hap.py](https://github.com/Illumina/hap.py) and [RTG Tools](https://github.com/RealTimeGenomics/rtg-tools/tree/eb13bbb82d2fbeab7d54a92e8493ddd2acf0d349)

```bash
conda install -c conda-forge cython=0.29.15
conda install -c conda-forge scipy=1.2.1 # Also installs numpy=1.16.5 dependency
conda install -c conda-forge pandas=0.24.2
conda install -c bioconda pybedtools=0.8.1 # Also installs pysam=0.15.3 dependency
conda install -c bioconda bx-python=0.8.8
```

Install hap.py and RTG Tools (using the python helper script)

```bash
cd /store/lkemp/exome_project/quality_benchmarking/
git clone git@github.com:Illumina/hap.py.git
cd hap.py
python install.py /store/lkemp/exome_project/quality_benchmarking/hap.py-install \
--with-rtgtools \
--no-tests
```

Install other benchmarking software in a new conda environment

```bash
conda deactivate
conda create -n benchmarking_env python=3.7
conda activate benchmarking_env
conda install -c bioconda bedtools=2.29.2
conda install -c bioconda bcftools=1.10.2
conda install -c bioconda bedops=2.4.39
conda install -c bioconda gatk=3.8
```

### Download and prepare data

#### Fastq data to run through pipelines

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

*Note. NIST7035 and NIST7086 represent 2 vials of the same sample (NA12878)*

#### Truth vcf to compare output to

```bash
# VCF
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST.hc.snps.indels.vcf
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST.hc.snps.indels.vcf.idx
# regions bed
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz
```

Create separate vcf files for columns in known vcf (NIST7035 and NIST7086)

```bash
bcftools view \
-s NIST7035 \
project.NIST.hc.snps.indels.vcf \
-o project.NIST.hc.snps.indels.NIST7035.vcf \
--exclude-private

bcftools view \
-s NIST7086 \
project.NIST.hc.snps.indels.vcf \
-o project.NIST.hc.snps.indels.NIST7086.vcf \
--exclude-private
```

See [here](https://github.com/ga4gh/benchmarking-tools/blob/master/resources/high-confidence-sets/giab.md) for more info on GIAB resources

### Other setup

#### sdf files for hap.py

Create an sdf file for the reference human genome (that was used in the benchmarking pipeline run). Create this with the RTG-tools format function:

- ucsc.hg19.fasta

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/rtg format \
--output /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
/store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta
```

## Benchmarking

human_genomic_pipeline will undertake pre-processing and variant calling. Because variant filtering occurs with vcf annotation_pipeline, we will benchmark the vcf files that have gone through both pipelines. However, we will use a minimal version of the vcf_annotation_pipeline since the annotation rules are not required for benchmarking.

### intra_truth_comparison

Compare the two truth sample vcfs (NIST7035 and NIST7086)

bgzip and index vcfs for comparison

```bash
bgzip < project.NIST.hc.snps.indels.NIST7035.vcf > project.NIST.hc.snps.indels.NIST7035.vcf.gz
tabix project.NIST.hc.snps.indels.NIST7035.vcf.gz
bgzip < project.NIST.hc.snps.indels.NIST7086.vcf > project.NIST.hc.snps.indels.NIST7086.vcf.gz
tabix project.NIST.hc.snps.indels.NIST7086.vcf.gz
```

NIST7035 ('baseline') compared to NIST7086 ('truth')

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/intra_truth_comparison/
mkdir happy_project.NIST.hc.snps.indels.NIST7035_v_project.NIST.hc.snps.indels.NIST7086
cd happy_project.NIST.hc.snps.indels.NIST7035_v_project.NIST.hc.snps.indels.NIST7086
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_project.NIST.hc.snps.indels.NIST7035_v_project.NIST.hc.snps.indels.NIST7086 \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

NIST7086 ('baseline') compared to NIST7035 ('truth')

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/intra_truth_comparison/
mkdir happy_project.NIST.hc.snps.indels.NIST7086_v_project.NIST.hc.snps.indels.NIST7035
cd happy_project.NIST.hc.snps.indels.NIST7086_v_project.NIST.hc.snps.indels.NIST7035
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_project.NIST.hc.snps.indels.NIST7086_v_project.NIST.hc.snps.indels.NIST7035 \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

### quality_bench1.0

#### human_genomics_pipeline + minimal vcf_annotation_pipeline

- NIST7035

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.0/
mkdir happy_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035
cd happy_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

- NIST7086

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.0/
mkdir happy_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086
cd happy_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
../vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

#### parabricks germline pipeline

*Output vcfs were transferred to Wintermute*

bgzip and index vcfs for comparison

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.0/

bgzip < ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf > ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz
bgzip < ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf > ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz
```

- NIST7035

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.0/
mkdir happy_NIST7035_NIST_v_project.NIST.hc.snps.indels.NIST7035
cd happy_NIST7035_NIST_v_project.NIST.hc.snps.indels.NIST7035
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../parabricks/NIST7035_NIST.vcf \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7035_NIST_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

- NIST7086

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.0/
mkdir happy_NIST7086_NIST_v_project.NIST.hc.snps.indels.NIST7086
cd happy_NIST7086_NIST_v_project.NIST.hc.snps.indels.NIST7086
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
../parabricks/NIST7086_NIST.vcf \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7086_NIST_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

### quality_bench1.1

bgzip and index vcfs for comparison

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.1/
# Pipeline output vcf (NIST7035)
bgzip < ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf > ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz
# Pipeline output vcf (NIST7086)
bgzip < ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf > ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz
```

#### human_genomics_pipeline + minimal vcf_annotation_pipeline

- NIST7035

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.1/
mkdir happy_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035
cd happy_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

- NIST7086

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.1/
mkdir happy_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086
cd happy_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
../vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

### quality_bench1.2

#### human_genomics_pipeline + minimal vcf_annotation_pipeline

Extract snps and indels based on filter tranche levels in vcf pipeline output (stored in the 'FILTER' column) so they can be compared separately (see [here](https://gatkforums.broadinstitute.org/gatk/discussion/1255/using-jexl-to-apply-hard-filters-or-select-variants-based-on-annotation-values) and [here](https://gatkforums.broadinstitute.org/gatk/discussion/12406/selectvariants-from-filter-column-gatk4) for info on passing selection conditions to gatk SelectVariants)

NIST7035

- SNPS

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.2/vcf_annotation_pipeline/filtered/

# Extract all snps
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered.vcf \
--select-type-to-include SNP \
-O NIST7035_NIST_filtered_SNP.vcf

# Extract all PASS
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered_SNP.vcf \
-select "vc.isNotFiltered()" \
-O NIST7035_NIST_filtered_SNP_PASS.vcf

# Extract all FILTERED
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered_SNP.vcf \
-select "vc.isFiltered()" \
-O NIST7035_NIST_filtered_SNP_FILTERED.vcf

# Extract for each tranche
for i in 'CNN_2D_SNP_Tranche_99.00_100.00' 'CNN_2D_SNP_Tranche_98.00_99.00' 'CNN_2D_SNP_Tranche_97.00_98.00' 'CNN_2D_SNP_Tranche_96.00_97.00' 'CNN_2D_SNP_Tranche_95.00_96.00' 'CNN_2D_SNP_Tranche_94.00_95.00' 'CNN_2D_SNP_Tranche_93.00_94.00' 'CNN_2D_SNP_Tranche_92.00_93.00' 'CNN_2D_SNP_Tranche_91.00_92.00' 'CNN_2D_SNP_Tranche_90.00_91.00' 'CNN_2D_SNP_Tranche_89.00_90.00' 'CNN_2D_SNP_Tranche_88.00_89.00' 'CNN_2D_SNP_Tranche_87.00_88.00' 'CNN_2D_SNP_Tranche_86.00_87.00' 'CNN_2D_SNP_Tranche_85.00_86.00' 'CNN_2D_SNP_Tranche_84.00_85.00' 'CNN_2D_SNP_Tranche_83.00_84.00' 'CNN_2D_SNP_Tranche_82.00_83.00' 'CNN_2D_SNP_Tranche_81.00_82.00' 'CNN_2D_SNP_Tranche_80.00_81.00'
do gatk SelectVariants \
  -R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
  -V NIST7035_NIST_filtered.vcf \
  --select "vc.getFilters().contains('$i')" \
  -O NIST7035_NIST_filtered_$i.vcf
done

# Get number of variants in each tranche
for i in 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_99.00_100.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_98.00_99.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_97.00_98.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_96.00_97.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_95.00_96.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_94.00_95.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_93.00_94.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_92.00_93.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_91.00_92.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_90.00_91.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_89.00_90.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_88.00_89.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_87.00_88.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_86.00_87.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_85.00_86.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_84.00_85.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_83.00_84.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_82.00_83.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_81.00_82.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_80.00_81.00.vcf'
do cat $i | grep -v '##' | wc -l
done
```

- INDELS

```bash
# Extract all indels
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered.vcf \
--select-type-to-include INDEL \
-O NIST7035_NIST_filtered_INDEL.vcf

# Extract all PASS
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered_INDEL.vcf \
-select "vc.isNotFiltered()" \
-O NIST7035_NIST_filtered_INDEL_PASS.vcf

# Extract all FILTERED
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered_INDEL.vcf \
-select "vc.isFiltered()" \
-O NIST7035_NIST_filtered_INDEL_FILTERED.vcf

# Extract for each tranche
for i in 'CNN_2D_INDEL_Tranche_99.00_100.00' 'CNN_2D_INDEL_Tranche_98.00_99.00' 'CNN_2D_INDEL_Tranche_97.00_98.00' 'CNN_2D_INDEL_Tranche_96.00_97.00' 'CNN_2D_INDEL_Tranche_95.00_96.00' 'CNN_2D_INDEL_Tranche_94.00_95.00' 'CNN_2D_INDEL_Tranche_93.00_94.00' 'CNN_2D_INDEL_Tranche_92.00_93.00' 'CNN_2D_INDEL_Tranche_91.00_92.00' 'CNN_2D_INDEL_Tranche_90.00_91.00' 'CNN_2D_INDEL_Tranche_89.00_90.00' 'CNN_2D_INDEL_Tranche_88.00_89.00' 'CNN_2D_INDEL_Tranche_87.00_88.00' 'CNN_2D_INDEL_Tranche_86.00_87.00' 'CNN_2D_INDEL_Tranche_85.00_86.00' 'CNN_2D_INDEL_Tranche_84.00_85.00' 'CNN_2D_INDEL_Tranche_83.00_84.00' 'CNN_2D_INDEL_Tranche_82.00_83.00' 'CNN_2D_INDEL_Tranche_81.00_82.00' 'CNN_2D_INDEL_Tranche_80.00_81.00'
do gatk SelectVariants \
  -R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
  -V NIST7035_NIST_filtered.vcf \
  --select "vc.getFilters().contains('$i')" \
  -O NIST7035_NIST_filtered_$i.vcf
done

# Get number of variants in each tranche
for i in 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_99.00_100.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_98.00_99.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_97.00_98.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_96.00_97.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_95.00_96.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_94.00_95.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_93.00_94.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_92.00_93.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_91.00_92.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_90.00_91.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_89.00_90.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_88.00_89.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_87.00_88.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_86.00_87.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_85.00_86.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_84.00_85.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_83.00_84.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_82.00_83.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_81.00_82.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_80.00_81.00.vcf'
do cat $i | grep -v '##' | wc -l
done
```

NIST7086

- SNPS

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.2/vcf_annotation_pipeline/filtered/

# Extract all snps
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered.vcf \
--select-type-to-include SNP \
-O NIST7086_NIST_filtered_SNP.vcf

# Extract all PASS
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered_SNP.vcf \
-select "vc.isNotFiltered()" \
-O NIST7086_NIST_filtered_SNP_PASS.vcf

# Extract all FILTERED
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered_SNP.vcf \
-select "vc.isFiltered()" \
-O NIST7086_NIST_filtered_SNP_FILTERED.vcf

# Extract for each tranche
for i in 'CNN_2D_SNP_Tranche_99.00_100.00' 'CNN_2D_SNP_Tranche_98.00_99.00' 'CNN_2D_SNP_Tranche_97.00_98.00' 'CNN_2D_SNP_Tranche_96.00_97.00' 'CNN_2D_SNP_Tranche_95.00_96.00' 'CNN_2D_SNP_Tranche_94.00_95.00' 'CNN_2D_SNP_Tranche_93.00_94.00' 'CNN_2D_SNP_Tranche_92.00_93.00' 'CNN_2D_SNP_Tranche_91.00_92.00' 'CNN_2D_SNP_Tranche_90.00_91.00' 'CNN_2D_SNP_Tranche_89.00_90.00' 'CNN_2D_SNP_Tranche_88.00_89.00' 'CNN_2D_SNP_Tranche_87.00_88.00' 'CNN_2D_SNP_Tranche_86.00_87.00' 'CNN_2D_SNP_Tranche_85.00_86.00' 'CNN_2D_SNP_Tranche_84.00_85.00' 'CNN_2D_SNP_Tranche_83.00_84.00' 'CNN_2D_SNP_Tranche_82.00_83.00' 'CNN_2D_SNP_Tranche_81.00_82.00' 'CNN_2D_SNP_Tranche_80.00_81.00'
do gatk SelectVariants \
  -R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
  -V NIST7086_NIST_filtered.vcf \
  --select "vc.getFilters().contains('$i')" \
  -O NIST7086_NIST_filtered_$i.vcf
done

# Get number of variants in each tranche
for i in 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_99.00_100.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_98.00_99.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_97.00_98.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_96.00_97.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_95.00_96.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_94.00_95.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_93.00_94.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_92.00_93.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_91.00_92.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_90.00_91.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_89.00_90.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_88.00_89.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_87.00_88.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_86.00_87.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_85.00_86.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_84.00_85.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_83.00_84.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_82.00_83.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_81.00_82.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_80.00_81.00.vcf'
do cat $i | grep -v '##' | wc -l
done
```

- INDELS

```bash
# Extract all indels
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered.vcf \
--select-type-to-include INDEL \
-O NIST7086_NIST_filtered_INDEL.vcf

# Extract all PASS
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered_INDEL.vcf \
-select "vc.isNotFiltered()" \
-O NIST7086_NIST_filtered_INDEL_PASS.vcf

# Extract all FILTERED
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered_INDEL.vcf \
-select "vc.isFiltered()" \
-O NIST7086_NIST_filtered_INDEL_FILTERED.vcf

# Extract for each tranche
for i in 'CNN_2D_INDEL_Tranche_99.00_100.00' 'CNN_2D_INDEL_Tranche_98.00_99.00' 'CNN_2D_INDEL_Tranche_97.00_98.00' 'CNN_2D_INDEL_Tranche_96.00_97.00' 'CNN_2D_INDEL_Tranche_95.00_96.00' 'CNN_2D_INDEL_Tranche_94.00_95.00' 'CNN_2D_INDEL_Tranche_93.00_94.00' 'CNN_2D_INDEL_Tranche_92.00_93.00' 'CNN_2D_INDEL_Tranche_91.00_92.00' 'CNN_2D_INDEL_Tranche_90.00_91.00' 'CNN_2D_INDEL_Tranche_89.00_90.00' 'CNN_2D_INDEL_Tranche_88.00_89.00' 'CNN_2D_INDEL_Tranche_87.00_88.00' 'CNN_2D_INDEL_Tranche_86.00_87.00' 'CNN_2D_INDEL_Tranche_85.00_86.00' 'CNN_2D_INDEL_Tranche_84.00_85.00' 'CNN_2D_INDEL_Tranche_83.00_84.00' 'CNN_2D_INDEL_Tranche_82.00_83.00' 'CNN_2D_INDEL_Tranche_81.00_82.00' 'CNN_2D_INDEL_Tranche_80.00_81.00'
do gatk SelectVariants \
  -R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
  -V NIST7086_NIST_filtered.vcf \
  --select "vc.getFilters().contains('$i')" \
  -O NIST7086_NIST_filtered_$i.vcf
done

# Get number of variants in each tranche
for i in 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_99.00_100.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_98.00_99.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_97.00_98.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_96.00_97.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_95.00_96.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_94.00_95.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_93.00_94.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_92.00_93.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_91.00_92.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_90.00_91.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_89.00_90.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_88.00_89.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_87.00_88.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_86.00_87.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_85.00_86.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_84.00_85.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_83.00_84.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_82.00_83.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_81.00_82.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_80.00_81.00.vcf'
do cat $i | grep -v '##' | wc -l
done
```

This data was then explored in R - see [here](evaluation_of_tranch_filtering/evaluation_of_tranch_filtering.md)

From these results, I wanted to filter the data from a tranche of around 82.00. I'll extract variants (snps and indels) that are in a tranche lower than 82.00 or are a pass to compare against the truth vcf

```bash
# NIST7035
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_81.00_82.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_80.00_81.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_81.00_82.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_80.00_81.00') || vc.isNotFiltered()" \
-O NIST7035_NIST_filtered_less_than_82.00.vcf

# NIST7086
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_81.00_82.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_80.00_81.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_81.00_82.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_80.00_81.00') || vc.isNotFiltered()" \
-O NIST7086_NIST_filtered_less_than_82.00.vcf
```

bgzip and index vcfs for comparison

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.2/
# Pipeline output vcf (NIST7035)
bgzip < ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf > ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz
bgzip < ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_82.00.vcf > ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_82.00.vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_82.00.vcf.gz
# Pipeline output vcf (NIST7086)
bgzip < ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf > ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz
bgzip < ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered_less_than_82.00.vcf > ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered_less_than_82.00.vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered_less_than_82.00.vcf.gz
```

- NIST7035

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.2/
mkdir happy_NIST7035_NIST_filtered_less_than_82.00_v_project.NIST.hc.snps.indels.NIST7035
cd happy_NIST7035_NIST_filtered_less_than_82.00_v_project.NIST.hc.snps.indels.NIST7035
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_82.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7035_NIST_filtered_less_than_82.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

- NIST7086

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.2/
mkdir happy_NIST7086_NIST_filtered_less_than_82.00_v_project.NIST.hc.snps.indels.NIST7086
cd happy_NIST7086_NIST_filtered_less_than_82.00_v_project.NIST.hc.snps.indels.NIST7086
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
../vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered_less_than_82.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7086_NIST_filtered_less_than_82.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

### quality_bench1.3

#### human_genomics_pipeline + minimal vcf_annotation_pipeline

Extract snps and indels based on filter tranche levels in vcf pipeline output (stored in the 'FILTER' column) so they can be compared separately (see [here](https://gatkforums.broadinstitute.org/gatk/discussion/1255/using-jexl-to-apply-hard-filters-or-select-variants-based-on-annotation-values) and [here](https://gatkforums.broadinstitute.org/gatk/discussion/12406/selectvariants-from-filter-column-gatk4) for info on passing selection conditions to gatk SelectVariants)

NIST7035

- SNPS

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.3/vcf_annotation_pipeline/filtered/

# Extract all snps
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered.vcf \
--select-type-to-include SNP \
-O NIST7035_NIST_filtered_SNP.vcf

# Extract all PASS
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered_SNP.vcf \
-select "vc.isNotFiltered()" \
-O NIST7035_NIST_filtered_SNP_PASS.vcf

# Extract all FILTERED
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered_SNP.vcf \
-select "vc.isFiltered()" \
-O NIST7035_NIST_filtered_SNP_FILTERED.vcf

# Extract for each tranche
for i in 'CNN_2D_SNP_Tranche_79.00_80.00' 'CNN_2D_SNP_Tranche_78.00_79.00' 'CNN_2D_SNP_Tranche_77.00_78.00' 'CNN_2D_SNP_Tranche_76.00_77.00' 'CNN_2D_SNP_Tranche_75.00_76.00' 'CNN_2D_SNP_Tranche_74.00_75.00' 'CNN_2D_SNP_Tranche_73.00_74.00' 'CNN_2D_SNP_Tranche_72.00_73.00' 'CNN_2D_SNP_Tranche_71.00_72.00' 'CNN_2D_SNP_Tranche_70.00_71.00'
do gatk SelectVariants \
  -R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
  -V NIST7035_NIST_filtered.vcf \
  --select "vc.getFilters().contains('$i')" \
  -O NIST7035_NIST_filtered_$i.vcf
done

# Get number of variants in each tranche
for i in 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_79.00_80.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_78.00_79.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_77.00_78.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_76.00_77.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_75.00_76.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_74.00_75.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_73.00_74.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_72.00_73.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_71.00_72.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_SNP_Tranche_70.00_71.00.vcf'
do cat $i | grep -v '##' | wc -l
done
```

- INDELS

```bash
# Extract all indels
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered.vcf \
--select-type-to-include INDEL \
-O NIST7035_NIST_filtered_INDEL.vcf

# Extract all PASS
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered_INDEL.vcf \
-select "vc.isNotFiltered()" \
-O NIST7035_NIST_filtered_INDEL_PASS.vcf

# Extract all FILTERED
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered_INDEL.vcf \
-select "vc.isFiltered()" \
-O NIST7035_NIST_filtered_INDEL_FILTERED.vcf

# Extract for each tranche
for i in 'CNN_2D_INDEL_Tranche_79.00_80.00' 'CNN_2D_INDEL_Tranche_78.00_79.00' 'CNN_2D_INDEL_Tranche_77.00_78.00' 'CNN_2D_INDEL_Tranche_76.00_77.00' 'CNN_2D_INDEL_Tranche_75.00_76.00' 'CNN_2D_INDEL_Tranche_74.00_75.00' 'CNN_2D_INDEL_Tranche_73.00_74.00' 'CNN_2D_INDEL_Tranche_72.00_73.00' 'CNN_2D_INDEL_Tranche_71.00_72.00' 'CNN_2D_INDEL_Tranche_70.00_71.00'
do gatk SelectVariants \
  -R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
  -V NIST7035_NIST_filtered.vcf \
  --select "vc.getFilters().contains('$i')" \
  -O NIST7035_NIST_filtered_$i.vcf
done

# Get number of variants in each tranche
for i in 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_79.00_80.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_78.00_79.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_77.00_78.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_76.00_77.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_75.00_76.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_74.00_75.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_73.00_74.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_72.00_73.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_71.00_72.00.vcf' 'NIST7035_NIST_filtered_CNN_2D_INDEL_Tranche_70.00_71.00.vcf'
do cat $i | grep -v '##' | wc -l
done
```

NIST7086

- SNPS

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.3/vcf_annotation_pipeline/filtered/

# Extract all snps
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered.vcf \
--select-type-to-include SNP \
-O NIST7086_NIST_filtered_SNP.vcf

# Extract all PASS
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered_SNP.vcf \
-select "vc.isNotFiltered()" \
-O NIST7086_NIST_filtered_SNP_PASS.vcf

# Extract all FILTERED
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered_SNP.vcf \
-select "vc.isFiltered()" \
-O NIST7086_NIST_filtered_SNP_FILTERED.vcf

# Extract for each tranche
for i in 'CNN_2D_SNP_Tranche_79.00_80.00' 'CNN_2D_SNP_Tranche_78.00_79.00' 'CNN_2D_SNP_Tranche_77.00_78.00' 'CNN_2D_SNP_Tranche_76.00_77.00' 'CNN_2D_SNP_Tranche_75.00_76.00' 'CNN_2D_SNP_Tranche_74.00_75.00' 'CNN_2D_SNP_Tranche_73.00_74.00' 'CNN_2D_SNP_Tranche_72.00_73.00' 'CNN_2D_SNP_Tranche_71.00_72.00' 'CNN_2D_SNP_Tranche_70.00_71.00'
do gatk SelectVariants \
  -R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
  -V NIST7086_NIST_filtered.vcf \
  --select "vc.getFilters().contains('$i')" \
  -O NIST7086_NIST_filtered_$i.vcf
done

# Get number of variants in each tranche
for i in 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_79.00_80.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_78.00_79.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_77.00_78.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_76.00_77.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_75.00_76.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_74.00_75.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_73.00_74.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_72.00_73.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_71.00_72.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_SNP_Tranche_70.00_71.00.vcf'
do cat $i | grep -v '##' | wc -l
done
```

- INDELS

```bash
# Extract all indels
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered.vcf \
--select-type-to-include INDEL \
-O NIST7086_NIST_filtered_INDEL.vcf

# Extract all PASS
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered_INDEL.vcf \
-select "vc.isNotFiltered()" \
-O NIST7086_NIST_filtered_INDEL_PASS.vcf

# Extract all FILTERED
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered_INDEL.vcf \
-select "vc.isFiltered()" \
-O NIST7086_NIST_filtered_INDEL_FILTERED.vcf

# Extract for each tranche
for i in 'CNN_2D_INDEL_Tranche_79.00_80.00' 'CNN_2D_INDEL_Tranche_78.00_79.00' 'CNN_2D_INDEL_Tranche_77.00_78.00' 'CNN_2D_INDEL_Tranche_76.00_77.00' 'CNN_2D_INDEL_Tranche_75.00_76.00' 'CNN_2D_INDEL_Tranche_74.00_75.00' 'CNN_2D_INDEL_Tranche_73.00_74.00' 'CNN_2D_INDEL_Tranche_72.00_73.00' 'CNN_2D_INDEL_Tranche_71.00_72.00' 'CNN_2D_INDEL_Tranche_70.00_71.00'
do gatk SelectVariants \
  -R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
  -V NIST7086_NIST_filtered.vcf \
  --select "vc.getFilters().contains('$i')" \
  -O NIST7086_NIST_filtered_$i.vcf
done

# Get number of variants in each tranche
for i in 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_79.00_80.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_78.00_79.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_77.00_78.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_76.00_77.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_75.00_76.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_74.00_75.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_73.00_74.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_72.00_73.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_71.00_72.00.vcf' 'NIST7086_NIST_filtered_CNN_2D_INDEL_Tranche_70.00_71.00.vcf'
do cat $i | grep -v '##' | wc -l
done
```

This data was then explored in R - see [here](evaluation_of_tranch_filtering/evaluation_of_tranch_filtering.md)

From these results, I wanted to filter the data from a tranche of around 82.00. I'll extract variants (snps and indels) that are in a tranche lower than 82.00 or are a pass to compare against the truth vcf. I'll also I'll extract variants (snps and indels) that are in a tranche lower than 70.00 or are a pass.

```bash
# NIST7035
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_81.00_82.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_80.00_81.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_79.00_80.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_78.00_79.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_77.00_78.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_76.00_77.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_75.00_76.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_74.00_75.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_73.00_74.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_72.00_73.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_71.00_72.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_70.00_71.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_81.00_82.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_80.00_81.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_79.00_80.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_78.00_79.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_77.00_78.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_76.00_77.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_75.00_76.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_74.00_75.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_73.00_74.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_72.00_73.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_71.00_72.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_70.00_71.00') || vc.isNotFiltered()" \
-O NIST7035_NIST_filtered_less_than_82.00.vcf

gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_70.00_71.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_70.00_71.00') || vc.isNotFiltered()" \
-O NIST7035_NIST_filtered_less_than_70.00.vcf

# NIST7086
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_81.00_82.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_80.00_81.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_79.00_80.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_78.00_79.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_77.00_78.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_76.00_77.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_75.00_76.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_74.00_75.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_73.00_74.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_72.00_73.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_71.00_72.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_70.00_71.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_81.00_82.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_80.00_81.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_79.00_80.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_78.00_79.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_77.00_78.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_76.00_77.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_75.00_76.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_74.00_75.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_73.00_74.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_72.00_73.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_71.00_72.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_70.00_71.00') || vc.isNotFiltered()" \
-O NIST7086_NIST_filtered_less_than_82.00.vcf

gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_70.00_71.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_70.00_71.00') || vc.isNotFiltered()" \
-O NIST7086_NIST_filtered_less_than_70.00.vcf
```

bgzip and index vcfs for comparison

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.3/
# Pipeline output vcf (NIST7035)
bgzip < ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf > ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz
bgzip < ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_82.00..vcf > ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_82.00..vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_82.00..vcf.gz
bgzip < ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_70.00.vcf > ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_70.00.vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_70.00.vcf.gz
# Pipeline output vcf (NIST7086)
bgzip < ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf > ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz
bgzip < ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered_less_than_82.00..vcf > ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered_less_than_82.00..vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered_less_than_82.00..vcf.gz
bgzip < ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered_less_than_70.00.vcf > ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered_less_than_70.00.vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered_less_than_70.00.vcf.gz
```

I also want to compare the filtering of the other tranches, I'll extract them so they can be compared to the known vcf with hap.py

```bash
# NIST7035
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_71.00_72.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_70.00_71.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_71.00_72.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_70.00_71.00') || vc.isNotFiltered()" \
-O NIST7035_NIST_filtered_less_than_72.00.vcf

gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_73.00_74.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_72.00_73.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_71.00_72.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_70.00_71.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_73.00_74.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_72.00_73.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_71.00_72.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_70.00_71.00') || vc.isNotFiltered()" \
-O NIST7035_NIST_filtered_less_than_74.00.vcf

gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_75.00_76.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_74.00_75.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_73.00_74.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_72.00_73.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_71.00_72.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_70.00_71.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_75.00_76.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_74.00_75.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_73.00_74.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_72.00_73.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_71.00_72.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_70.00_71.00') || vc.isNotFiltered()" \
-O NIST7035_NIST_filtered_less_than_76.00.vcf

gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_77.00_78.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_76.00_77.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_75.00_76.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_74.00_75.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_73.00_74.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_72.00_73.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_71.00_72.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_70.00_71.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_77.00_78.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_76.00_77.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_75.00_76.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_74.00_75.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_73.00_74.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_72.00_73.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_71.00_72.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_70.00_71.00') || vc.isNotFiltered()" \
-O NIST7035_NIST_filtered_less_than_78.00.vcf

gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_79.00_80.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_78.00_79.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_77.00_78.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_76.00_77.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_75.00_76.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_74.00_75.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_73.00_74.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_72.00_73.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_71.00_72.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_70.00_71.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_79.00_80.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_78.00_79.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_77.00_78.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_76.00_77.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_75.00_76.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_74.00_75.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_73.00_74.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_72.00_73.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_71.00_72.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_70.00_71.00') || vc.isNotFiltered()" \
-O NIST7035_NIST_filtered_less_than_80.00.vcf

# NIST7086
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_71.00_72.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_70.00_71.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_71.00_72.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_70.00_71.00') || vc.isNotFiltered()" \
-O NIST7086_NIST_filtered_less_than_72.00.vcf

gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_73.00_74.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_72.00_73.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_71.00_72.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_70.00_71.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_73.00_74.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_72.00_73.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_71.00_72.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_70.00_71.00') || vc.isNotFiltered()" \
-O NIST7086_NIST_filtered_less_than_74.00.vcf

gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_75.00_76.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_74.00_75.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_73.00_74.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_72.00_73.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_71.00_72.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_70.00_71.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_75.00_76.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_74.00_75.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_73.00_74.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_72.00_73.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_71.00_72.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_70.00_71.00') || vc.isNotFiltered()"  \
-O NIST7086_NIST_filtered_less_than_76.00.vcf

gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_77.00_78.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_76.00_77.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_75.00_76.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_74.00_75.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_73.00_74.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_72.00_73.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_71.00_72.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_70.00_71.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_77.00_78.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_76.00_77.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_75.00_76.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_74.00_75.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_73.00_74.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_72.00_73.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_71.00_72.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_70.00_71.00') || vc.isNotFiltered()" \
-O NIST7086_NIST_filtered_less_than_78.00.vcf

gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_79.00_80.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_78.00_79.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_77.00_78.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_76.00_77.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_75.00_76.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_74.00_75.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_73.00_74.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_72.00_73.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_71.00_72.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_70.00_71.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_79.00_80.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_78.00_79.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_77.00_78.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_76.00_77.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_75.00_76.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_74.00_75.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_73.00_74.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_72.00_73.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_71.00_72.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_70.00_71.00') || vc.isNotFiltered()" \
-O NIST7086_NIST_filtered_less_than_80.00.vcf
```

bgzip and index vcfs for comparison

```bash
for i in 'NIST7035_NIST_filtered_less_than_72.00.vcf' 'NIST7035_NIST_filtered_less_than_74.00.vcf' 'NIST7035_NIST_filtered_less_than_76.00.vcf' 'NIST7035_NIST_filtered_less_than_78.00.vcf' 'NIST7035_NIST_filtered_less_than_80.00.vcf' 'NIST7086_NIST_filtered_less_than_72.00.vcf' 'NIST7086_NIST_filtered_less_than_74.00.vcf' 'NIST7086_NIST_filtered_less_than_76.00.vcf' 'NIST7086_NIST_filtered_less_than_78.00.vcf' 'NIST7086_NIST_filtered_less_than_80.00.vcf'
do bgzip $i
done

for i in 'NIST7035_NIST_filtered_less_than_72.00.vcf.gz' 'NIST7035_NIST_filtered_less_than_74.00.vcf.gz' 'NIST7035_NIST_filtered_less_than_76.00.vcf.gz' 'NIST7035_NIST_filtered_less_than_78.00.vcf.gz' 'NIST7035_NIST_filtered_less_than_80.00.vcf.gz' 'NIST7086_NIST_filtered_less_than_72.00.vcf.gz' 'NIST7086_NIST_filtered_less_than_74.00.vcf.gz' 'NIST7086_NIST_filtered_less_than_76.00.vcf.gz' 'NIST7086_NIST_filtered_less_than_78.00.vcf.gz' 'NIST7086_NIST_filtered_less_than_80.00.vcf.gz'
do tabix $i
done
```

- NIST7035

```bash
# Less than 82
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.3/
mkdir happy_NIST7035_NIST_filtered_less_than_82.00_v_project.NIST.hc.snps.indels.NIST7035
cd happy_NIST7035_NIST_filtered_less_than_82.00_v_project.NIST.hc.snps.indels.NIST7035

/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_82.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7035_NIST_filtered_less_than_82.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

```bash
# Less than 80
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.3/
mkdir happy_NIST7035_NIST_filtered_less_than_80.00_v_project.NIST.hc.snps.indels.NIST7035
cd happy_NIST7035_NIST_filtered_less_than_80.00_v_project.NIST.hc.snps.indels.NIST7035

/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_80.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7035_NIST_filtered_less_than_80.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

```bash
# Less than 78
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.3/
mkdir happy_NIST7035_NIST_filtered_less_than_78.00_v_project.NIST.hc.snps.indels.NIST7035
cd happy_NIST7035_NIST_filtered_less_than_78.00_v_project.NIST.hc.snps.indels.NIST7035

/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_78.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7035_NIST_filtered_less_than_78.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

```bash
# Less than 76
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.3/
mkdir happy_NIST7035_NIST_filtered_less_than_76.00_v_project.NIST.hc.snps.indels.NIST7035
cd happy_NIST7035_NIST_filtered_less_than_76.00_v_project.NIST.hc.snps.indels.NIST7035

/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_76.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7035_NIST_filtered_less_than_76.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

```bash
# Less than 74
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.3/
mkdir happy_NIST7035_NIST_filtered_less_than_74.00_v_project.NIST.hc.snps.indels.NIST7035
cd happy_NIST7035_NIST_filtered_less_than_74.00_v_project.NIST.hc.snps.indels.NIST7035

/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_74.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7035_NIST_filtered_less_than_74.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

```bash
# Less than 72
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.3/
mkdir happy_NIST7035_NIST_filtered_less_than_72.00_v_project.NIST.hc.snps.indels.NIST7035
cd happy_NIST7035_NIST_filtered_less_than_72.00_v_project.NIST.hc.snps.indels.NIST7035

/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_72.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7035_NIST_filtered_less_than_72.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

- NIST7086

```bash
# Less than 82
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.3/
mkdir happy_NIST7086_NIST_filtered_less_than_82.00_v_project.NIST.hc.snps.indels.NIST7086
cd happy_NIST7086_NIST_filtered_less_than_82.00_v_project.NIST.hc.snps.indels.NIST7086

/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
../vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered_less_than_82.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7086_NIST_filtered_less_than_82.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

```bash
# Less than 80
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.3/
mkdir happy_NIST7086_NIST_filtered_less_than_80.00_v_project.NIST.hc.snps.indels.NIST7086
cd happy_NIST7086_NIST_filtered_less_than_80.00_v_project.NIST.hc.snps.indels.NIST7086

/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
../vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered_less_than_80.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7086_NIST_filtered_less_than_80.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

```bash
# Less than 78
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.3/
mkdir happy_NIST7086_NIST_filtered_less_than_78.00_v_project.NIST.hc.snps.indels.NIST7086
cd happy_NIST7086_NIST_filtered_less_than_78.00_v_project.NIST.hc.snps.indels.NIST7086

/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
../vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered_less_than_78.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7086_NIST_filtered_less_than_78.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

```bash
# Less than 76
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.3/
mkdir happy_NIST7086_NIST_filtered_less_than_76.00_v_project.NIST.hc.snps.indels.NIST7086
cd happy_NIST7086_NIST_filtered_less_than_76.00_v_project.NIST.hc.snps.indels.NIST7086

/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
../vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered_less_than_76.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7086_NIST_filtered_less_than_76.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

```bash
# Less than 76
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.3/
mkdir happy_NIST7086_NIST_filtered_less_than_74.00_v_project.NIST.hc.snps.indels.NIST7086
cd happy_NIST7086_NIST_filtered_less_than_74.00_v_project.NIST.hc.snps.indels.NIST7086

/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
../vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered_less_than_74.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7086_NIST_filtered_less_than_74.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

```bash
# Less than 76
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.3/
mkdir happy_NIST7086_NIST_filtered_less_than_72.00_v_project.NIST.hc.snps.indels.NIST7086
cd happy_NIST7086_NIST_filtered_less_than_72.00_v_project.NIST.hc.snps.indels.NIST7086

/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
../vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered_less_than_72.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7086_NIST_filtered_less_than_72.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

### quality_bench1.4

bgzip and index vcfs for comparison

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.4/
# Pipeline output vcf (NIST7035)
bgzip < ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf > ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz
# Pipeline output vcf (NIST7086)
bgzip < ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf > ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz
```

#### human_genomics_pipeline + minimal vcf_annotation_pipeline

- NIST7035

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.4/
mkdir happy_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035
cd happy_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

- NIST7086

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.4/
mkdir happy_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086
cd happy_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
../vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

#### parabricks germline pipeline

*Output vcfs were transferred to Wintermute*

- NIST7035

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.4/
mkdir happy_NIST7035_NIST_v_project.NIST.hc.snps.indels.NIST7035
cd happy_NIST7035_NIST_v_project.NIST.hc.snps.indels.NIST7035
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../parabricks/NIST7035_NIST.vcf \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7035_NIST_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

- NIST7086

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.4/
mkdir happy_NIST7086_NIST_v_project.NIST.hc.snps.indels.NIST7086
cd happy_NIST7086_NIST_v_project.NIST.hc.snps.indels.NIST7086
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
../parabricks/NIST7086_NIST.vcf \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7086_NIST_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

### quality_bench1.5

#### no pipeline (intra_truth_comparison re-run)

- NIST7035 ('baseline') compared to NIST7086 ('truth')

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.5/
mkdir happy_intra_truth_comparison_re-run_project.NIST.hc.snps.indels.NIST7035_v_project.NIST.hc.snps.indels.NIST7086
cd happy_intra_truth_comparison_re-run_project.NIST.hc.snps.indels.NIST7035_v_project.NIST.hc.snps.indels.NIST7086
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_intra_truth_comparison_re-run_project.NIST.hc.snps.indels.NIST7035_v_project.NIST.hc.snps.indels.NIST7086 \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--type ga4gh \
--threads 16
```

- NIST7086 ('baseline') compared to NIST7035 ('truth')

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.5/
mkdir happy_intra_truth_comparison_re-run_project.NIST.hc.snps.indels.NIST7086_v_project.NIST.hc.snps.indels.NIST7035
cd happy_intra_truth_comparison_re-run_project.NIST.hc.snps.indels.NIST7086_v_project.NIST.hc.snps.indels.NIST7035
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_intra_truth_comparison_re-run_project.NIST.hc.snps.indels.NIST7086_v_project.NIST.hc.snps.indels.NIST7035 \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--type ga4gh \
--threads 16
```

#### human_genomics_pipeline + minimal vcf_annotation_pipeline (quality_bench1.0 re-run)

- NIST7035

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.5/
mkdir happy_quality_bench1.0_re-run_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035
cd happy_quality_bench1.0_re-run_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../../quality_bench1.0/vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_quality_bench1.0_re-run_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--type ga4gh \
--threads 16
```

- NIST7086

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.5/
mkdir happy_quality_bench1.0_re-run_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086
cd happy_quality_bench1.0_re-run_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
../../quality_bench1.0/vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_quality_bench1.0_re-run_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--type ga4gh \
--threads 16
```

#### parabricks germline pipeline (quality_bench1.0 re-run)

- NIST7035

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.5/
mkdir happy_quality_bench1.0_re-run_NIST7035_NIST_v_project.NIST.hc.snps.indels.NIST7035
cd happy_quality_bench1.0_re-run_NIST7035_NIST_v_project.NIST.hc.snps.indels.NIST7035
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../../quality_bench1.0/parabricks/NIST7035_NIST.vcf \
-f /store/lkemp/publicData/exomes/NA12878_exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_quality_bench1.0_re-run_NIST7035_NIST_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--type ga4gh \
--threads 16
```

- NIST7086

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.5/
mkdir happy_quality_bench1.0_re-run_NIST7086_NIST_v_project.NIST.hc.snps.indels.NIST7086
cd happy_quality_bench1.0_re-run_NIST7086_NIST_v_project.NIST.hc.snps.indels.NIST7086
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
../../quality_bench1.0/parabricks/NIST7086_NIST.vcf \
-f /store/lkemp/publicData/exomes/NA12878_exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_quality_bench1.0_re-run_NIST7086_NIST_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--type ga4gh \
--threads 16
```

#### human_genomics_pipeline + minimal vcf_annotation_pipeline (quality_bench1.1 re-run)

- NIST7035

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.5/
mkdir happy_quality_bench1.1_re-run_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035
cd happy_quality_bench1.1_re-run_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../../quality_bench1.1/vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_quality_bench1.1_re-run_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--type ga4gh \
--threads 16
```

- NIST7086

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.5/
mkdir happy_quality_bench1.1_re-run_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086
cd happy_quality_bench1.1_re-run_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
../../quality_bench1.1/vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_quality_bench1.1_re-run_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--type ga4gh \
--threads 16
```

#### human_genomics_pipeline + minimal vcf_annotation_pipeline (quality_bench1.2 re-run)

- NIST7035

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.5/
mkdir happy_quality_bench1.2_re-run_NIST7035_NIST_filtered_less_than_82.00_v_project.NIST.hc.snps.indels.NIST7035
cd happy_quality_bench1.2_re-run_NIST7035_NIST_filtered_less_than_82.00_v_project.NIST.hc.snps.indels.NIST7035
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../../quality_bench1.2/vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_82.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_quality_bench1.2_re-run_NIST7035_NIST_filtered_less_than_82.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--type ga4gh \
--threads 16
```

- NIST7086

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.5/
mkdir happy_quality_bench1.2_re-run_NIST7086_NIST_filtered_less_than_82.00_v_project.NIST.hc.snps.indels.NIST7086
cd happy_quality_bench1.2_re-run_NIST7086_NIST_filtered_less_than_82.00_v_project.NIST.hc.snps.indels.NIST7086
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
../../quality_bench1.2/vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered_less_than_82.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_quality_bench1.2_re-run_NIST7086_NIST_filtered_less_than_82.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--type ga4gh \
--threads 16
```

#### human_genomics_pipeline + minimal vcf_annotation_pipeline (quality_bench1.3 re-run)

- NIST7035

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.5/
mkdir happy_quality_bench1.3_re-run_NIST7035_NIST_filtered_less_than_70.00_v_project.NIST.hc.snps.indels.NIST7035
cd happy_quality_bench1.3_re-run_NIST7035_NIST_filtered_less_than_70.00_v_project.NIST.hc.snps.indels.NIST7035
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../../quality_bench1.3/vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_70.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_quality_bench1.3_re-run_NIST7035_NIST_filtered_less_than_70.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--type ga4gh \
--threads 16
```

- NIST7086

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.5/
mkdir happy_quality_bench1.3_re-run_NIST7086_NIST_filtered_less_than_70.00_v_project.NIST.hc.snps.indels.NIST7086
cd happy_quality_bench1.3_re-run_NIST7086_NIST_filtered_less_than_70.00_v_project.NIST.hc.snps.indels.NIST7086
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
../../quality_bench1.3/vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered_less_than_70.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_quality_bench1.3_re-run_NIST7086_NIST_filtered_less_than_70.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--type ga4gh \
--threads 16
```

#### human_genomics_pipeline + minimal vcf_annotation_pipeline (quality_bench1.4 re-run)

- NIST7035

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.5/
mkdir happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035
cd happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../../quality_bench1.4/vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--type ga4gh \
--threads 16
```

- NIST7086

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.5/
mkdir happy_quality_bench1.4_re-run_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086
cd happy_quality_bench1.4_re-run_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
../../quality_bench1.4/vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_quality_bench1.4_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--type ga4gh \
--threads 16
```

#### parabricks germline pipeline (quality_bench1.4 re-run)

- NIST7035

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.5/
mkdir happy_quality_bench1.4_re-run_NIST7035_NIST_v_project.NIST.hc.snps.indels.NIST7035
cd happy_quality_bench1.4_re-run_NIST7035_NIST_v_project.NIST.hc.snps.indels.NIST7035
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../../quality_bench1.4/parabricks/NIST7035_NIST.vcf \
-f /store/lkemp/publicData/exomes/NA12878_exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_quality_bench1.4_NIST7035_NIST_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--type ga4gh \
--threads 16
```

- NIST7086

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.5/
mkdir happy_quality_bench1.4_re-run_NIST7086_NIST_v_project.NIST.hc.snps.indels.NIST7086
cd happy_quality_bench1.4_re-run_NIST7086_NIST_v_project.NIST.hc.snps.indels.NIST7086
```

```bash
/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
../../quality_bench1.4/parabricks/NIST7086_NIST.vcf \
-f /store/lkemp/publicData/exomes/NA12878_exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_quality_bench1.4_NIST7086_NIST_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--type ga4gh \
--threads 16
```

### quality_bench1.6

Extract snps and indels based on filter tranche levels in vcf pipeline output (stored in the 'FILTER' column) so they can be compared separately (see [here](https://gatkforums.broadinstitute.org/gatk/discussion/1255/using-jexl-to-apply-hard-filters-or-select-variants-based-on-annotation-values) and [here](https://gatkforums.broadinstitute.org/gatk/discussion/12406/selectvariants-from-filter-column-gatk4) for info on passing selection conditions to gatk SelectVariants)

#### human_genomics_pipeline + minimal vcf_annotation_pipeline

- NIST7035

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.4/vcf_annotation_pipeline/filtered/

# less than 99.0
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_93.00_94.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_94.00_95.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_95.00_96.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_96.00_97.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_97.00_98.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_98.00_99.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_93.00_94.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_94.00_95.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_95.00_96.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_96.00_97.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_97.00_98.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_98.00_99.00') || vc.isNotFiltered()" \
-O NIST7035_NIST_filtered_less_than_99.00.vcf

# less than 98.0
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_93.00_94.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_94.00_95.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_95.00_96.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_96.00_97.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_97.00_98.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_93.00_94.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_94.00_95.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_95.00_96.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_96.00_97.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_97.00_98.00') || vc.isNotFiltered()"  \
-O NIST7035_NIST_filtered_less_than_98.00.vcf

# less than 97.0
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_93.00_94.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_94.00_95.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_95.00_96.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_96.00_97.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_93.00_94.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_94.00_95.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_95.00_96.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_96.00_97.00') || vc.isNotFiltered()"  \
-O NIST7035_NIST_filtered_less_than_97.00.vcf

# less than 96.0
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_93.00_94.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_94.00_95.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_95.00_96.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_93.00_94.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_94.00_95.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_95.00_96.00') || vc.isNotFiltered()"  \
-O NIST7035_NIST_filtered_less_than_96.00.vcf

# less than 95.0
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_93.00_94.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_94.00_95.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_93.00_94.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_94.00_95.00') || vc.isNotFiltered()"  \
-O NIST7035_NIST_filtered_less_than_95.00.vcf

# less than 94.0
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_93.00_94.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_93.00_94.00') || vc.isNotFiltered()"  \
-O NIST7035_NIST_filtered_less_than_94.00.vcf

# less than 93.0
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_92.00_93.00') || vc.isNotFiltered()"  \
-O NIST7035_NIST_filtered_less_than_93.00.vcf

# less than 92.0
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_91.00_92.00') || vc.isNotFiltered()"  \
-O NIST7035_NIST_filtered_less_than_92.00.vcf

# less than 91.0
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_90.00_91.00') || vc.isNotFiltered()"  \
-O NIST7035_NIST_filtered_less_than_91.00.vcf

# less than 90.0
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7035_NIST_filtered.vcf \
--select "vc.isNotFiltered()"  \
-O NIST7035_NIST_filtered_less_than_90.00.vcf
```

- NIST7086

```bash
# less than 99.0
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_93.00_94.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_94.00_95.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_95.00_96.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_96.00_97.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_97.00_98.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_98.00_99.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_93.00_94.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_94.00_95.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_95.00_96.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_96.00_97.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_97.00_98.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_98.00_99.00') || vc.isNotFiltered()" \
-O NIST7086_NIST_filtered_less_than_99.00.vcf

# less than 98.0
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_93.00_94.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_94.00_95.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_95.00_96.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_96.00_97.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_97.00_98.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_93.00_94.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_94.00_95.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_95.00_96.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_96.00_97.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_97.00_98.00') || vc.isNotFiltered()"  \
-O NIST7086_NIST_filtered_less_than_98.00.vcf

# less than 97.0
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_93.00_94.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_94.00_95.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_95.00_96.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_96.00_97.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_93.00_94.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_94.00_95.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_95.00_96.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_96.00_97.00') || vc.isNotFiltered()"  \
-O NIST7086_NIST_filtered_less_than_97.00.vcf

# less than 96.0
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_93.00_94.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_94.00_95.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_95.00_96.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_93.00_94.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_94.00_95.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_95.00_96.00') || vc.isNotFiltered()"  \
-O NIST7086_NIST_filtered_less_than_96.00.vcf

# less than 95.0
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_93.00_94.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_94.00_95.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_93.00_94.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_94.00_95.00') || vc.isNotFiltered()"  \
-O NIST7086_NIST_filtered_less_than_95.00.vcf

# less than 94.0
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_93.00_94.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_93.00_94.00') || vc.isNotFiltered()"  \
-O NIST7086_NIST_filtered_less_than_94.00.vcf

# less than 93.0
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_92.00_93.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_92.00_93.00') || vc.isNotFiltered()"  \
-O NIST7086_NIST_filtered_less_than_93.00.vcf

# less than 92.0
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_SNP_Tranche_91.00_92.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_91.00_92.00') || vc.isNotFiltered()"  \
-O NIST7086_NIST_filtered_less_than_92.00.vcf

# less than 91.0
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered.vcf \
--select "vc.getFilters().contains('CNN_2D_SNP_Tranche_90.00_91.00') || vc.getFilters().contains('CNN_2D_INDEL_Tranche_90.00_91.00') || vc.isNotFiltered()"  \
-O NIST7086_NIST_filtered_less_than_91.00.vcf

# less than 90.0
gatk SelectVariants \
-R /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-V NIST7086_NIST_filtered.vcf \
--select "vc.isNotFiltered()"  \
-O NIST7086_NIST_filtered_less_than_90.00.vcf
```

bgzip and index vcfs for comparison

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.4/vcf_annotation_pipeline/filtered/


for i in 'NIST7035_NIST_filtered_less_than_99.00.vcf' 'NIST7035_NIST_filtered_less_than_98.00.vcf' 'NIST7035_NIST_filtered_less_than_97.00.vcf' 'NIST7035_NIST_filtered_less_than_96.00.vcf' 'NIST7035_NIST_filtered_less_than_95.00.vcf' 'NIST7035_NIST_filtered_less_than_94.00.vcf' 'NIST7035_NIST_filtered_less_than_93.00.vcf' 'NIST7035_NIST_filtered_less_than_92.00.vcf' 'NIST7035_NIST_filtered_less_than_91.00.vcf' 'NIST7035_NIST_filtered_less_than_90.00.vcf' 'NIST7086_NIST_filtered_less_than_99.00.vcf' 'NIST7086_NIST_filtered_less_than_98.00.vcf' 'NIST7086_NIST_filtered_less_than_97.00.vcf' 'NIST7086_NIST_filtered_less_than_96.00.vcf' 'NIST7086_NIST_filtered_less_than_95.00.vcf' 'NIST7086_NIST_filtered_less_than_94.00.vcf' 'NIST7086_NIST_filtered_less_than_93.00.vcf' 'NIST7086_NIST_filtered_less_than_92.00.vcf' 'NIST7086_NIST_filtered_less_than_91.00.vcf' 'NIST7086_NIST_filtered_less_than_90.00.vcf'
do bgzip $i
done

for i in 'NIST7035_NIST_filtered_less_than_99.00.vcf.gz' 'NIST7035_NIST_filtered_less_than_98.00.vcf.gz' 'NIST7035_NIST_filtered_less_than_97.00.vcf.gz' 'NIST7035_NIST_filtered_less_than_96.00.vcf.gz' 'NIST7035_NIST_filtered_less_than_95.00.vcf.gz' 'NIST7035_NIST_filtered_less_than_94.00.vcf.gz' 'NIST7035_NIST_filtered_less_than_93.00.vcf.gz' 'NIST7035_NIST_filtered_less_than_92.00.vcf.gz' 'NIST7035_NIST_filtered_less_than_91.00.vcf.gz' 'NIST7035_NIST_filtered_less_than_90.00.vcf.gz' 'NIST7086_NIST_filtered_less_than_99.00.vcf.gz' 'NIST7086_NIST_filtered_less_than_98.00.vcf.gz' 'NIST7086_NIST_filtered_less_than_97.00.vcf.gz' 'NIST7086_NIST_filtered_less_than_96.00.vcf.gz' 'NIST7086_NIST_filtered_less_than_95.00.vcf.gz' 'NIST7086_NIST_filtered_less_than_94.00.vcf.gz' 'NIST7086_NIST_filtered_less_than_93.00.vcf.gz' 'NIST7086_NIST_filtered_less_than_92.00.vcf.gz' 'NIST7086_NIST_filtered_less_than_91.00.vcf.gz' 'NIST7086_NIST_filtered_less_than_90.00.vcf.gz'
do tabix $i
done
```

- NIST7035

```bash
# less than 99.0
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.6/
mkdir happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_99.00_v_project.NIST.hc.snps.indels.NIST7035
cd happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_99.00_v_project.NIST.hc.snps.indels.NIST7035

/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../../quality_bench1.4/vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_99.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_99.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--type ga4gh \
--threads 16
```

```bash
# less than 98.0
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.6/
mkdir happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_98.00_v_project.NIST.hc.snps.indels.NIST7035
cd happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_98.00_v_project.NIST.hc.snps.indels.NIST7035

/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../../quality_bench1.4/vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_98.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_98.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--type ga4gh \
--threads 16
```

```bash
# less than 97.0
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.6/
mkdir happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_97.00_v_project.NIST.hc.snps.indels.NIST7035
cd happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_97.00_v_project.NIST.hc.snps.indels.NIST7035

/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../../quality_bench1.4/vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_97.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_97.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--type ga4gh \
--threads 16
```

```bash
# less than 96.0
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.6/
mkdir happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_96.00_v_project.NIST.hc.snps.indels.NIST7035
cd happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_96.00_v_project.NIST.hc.snps.indels.NIST7035

/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../../quality_bench1.4/vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_96.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_96.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--type ga4gh \
--threads 16
```

```bash
# less than 95.0
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.6/
mkdir happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_95.00_v_project.NIST.hc.snps.indels.NIST7035
cd happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_95.00_v_project.NIST.hc.snps.indels.NIST7035

/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../../quality_bench1.4/vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_95.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_95.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--type ga4gh \
--threads 16
```

```bash
# less than 94.0
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.6/
mkdir happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_94.00_v_project.NIST.hc.snps.indels.NIST7035
cd happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_94.00_v_project.NIST.hc.snps.indels.NIST7035

/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../../quality_bench1.4/vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_94.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_94.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--type ga4gh \
--threads 16
```

```bash
# less than 93.0
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.6/
mkdir happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_93.00_v_project.NIST.hc.snps.indels.NIST7035
cd happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_93.00_v_project.NIST.hc.snps.indels.NIST7035

/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../../quality_bench1.4/vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_93.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_93.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--type ga4gh \
--threads 16
```

```bash
# less than 92.0
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.6/
mkdir happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_92.00_v_project.NIST.hc.snps.indels.NIST7035
cd happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_92.00_v_project.NIST.hc.snps.indels.NIST7035

/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../../quality_bench1.4/vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_92.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_92.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--type ga4gh \
--threads 16
```

```bash
# less than 91.0
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.6/
mkdir happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_91.00_v_project.NIST.hc.snps.indels.NIST7035
cd happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_91.00_v_project.NIST.hc.snps.indels.NIST7035

/store/lkemp/exome_project/quality_benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../../quality_bench1.4/vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered_less_than_91.00.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_quality_bench1.4_re-run_NIST7035_NIST_filtered_less_than_91.00_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/quality_benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--type ga4gh \
--threads 16
```

- NIST7086

### quality_bench1.7

Create vcf files of false negatives (fn), false positives (fp) so they can be viewed in IGV to try an understand why they might be incorrectly called

#### human_genomics_pipeline + minimal vcf_annotation_pipeline

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.4/vcf_annotation_pipeline/filtered/

# False negatives
bedtools intersect \
-a /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
-b NIST7035_NIST_filtered.vcf \
-v \
-header \
> ../../fn_project.NIST.hc.snps.indels.NIST7035_v_NIST7035_NIST_filtered.vcf

bedtools intersect \
-a /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
-b NIST7086_NIST_filtered.vcf \
-v \
-header \
> ../../fn_project.NIST.hc.snps.indels.NIST7086_v_NIST7086_NIST_filtered.vcf

# False positives
bedtools intersect \
-a NIST7035_NIST_filtered.vcf \
-b /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
-v \
-header \
> ../../fp_project.NIST.hc.snps.indels.NIST7035_v_NIST7035_NIST_filtered.vcf

bedtools intersect \
-a NIST7086_NIST_filtered.vcf \
-b /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
-v \
-header \
> ../../fp_project.NIST.hc.snps.indels.NIST7086_v_NIST7086_NIST_filtered.vcf
```

Index the output vcfs

```bash
cd /store/lkemp/exome_project/quality_benchmarking/NA12878_exome/quality_bench1.4/

for i in 'fn_project.NIST.hc.snps.indels.NIST7035_v_NIST7035_NIST_filtered.vcf' 'fn_project.NIST.hc.snps.indels.NIST7086_v_NIST7086_NIST_filtered.vcf' 'fp_project.NIST.hc.snps.indels.NIST7035_v_NIST7035_NIST_filtered.vcf' 'fp_project.NIST.hc.snps.indels.NIST7086_v_NIST7086_NIST_filtered.vcf'
do bgzip $i
done

for i in 'fn_project.NIST.hc.snps.indels.NIST7035_v_NIST7035_NIST_filtered.vcf.gz' 'fn_project.NIST.hc.snps.indels.NIST7086_v_NIST7086_NIST_filtered.vcf.gz' 'fp_project.NIST.hc.snps.indels.NIST7035_v_NIST7035_NIST_filtered.vcf.gz' 'fp_project.NIST.hc.snps.indels.NIST7086_v_NIST7086_NIST_filtered.vcf.gz'
do tabix $i
done
```

Find variants in the vcf file to look at in IGV

```bash
# False negatives (missed variants) for NIST7035
zcat fn_project.NIST.hc.snps.indels.NIST7035_v_NIST7035_NIST_filtered.vcf.gz | zgrep -v '#' | head -n 20

# False positives for NIST7035
zcat fp_project.NIST.hc.snps.indels.NIST7035_v_NIST7035_NIST_filtered.vcf.gz | zgrep -v '#' | head -n 20
```
