# Create test dataset for pipelines

Created: 2020/10/29 10:07:39
Last modified: 2020/10/29 16:13:39

- **Aim:** Create a test dataset for [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) and [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline) (a trio and a singleton)
- **Prerequisite software:**  [Conda 4.8.5](https://docs.conda.io/projects/conda/en/latest/index.html)
- **OS:** CentOS-7 (ORAC - ESR cluster)

## Table of contents

- [Create test dataset for pipelines](#create-test-dataset-for-pipelines)
  - [Table of contents](#table-of-contents)
  - [Setup](#setup)
  - [Get bams and vcf for publicly available trio](#get-bams-and-vcf-for-publicly-available-trio)
  - [Manually create pedigree file](#manually-create-pedigree-file)
  - [Reduce dataset - subset by exome capture regions](#reduce-dataset---subset-by-exome-capture-regions)
  - [Prepare files and file directory structure to run through pipeline](#prepare-files-and-file-directory-structure-to-run-through-pipeline)
  - [Run through vcf_annotation_pipeline](#run-through-vcf_annotation_pipeline)
  - [Randomly sub-sample variants](#randomly-sub-sample-variants)
  - [Create a bed file from vcf](#create-a-bed-file-from-vcf)
  - [Pull out fastq reads from bam](#pull-out-fastq-reads-from-bam)
  - [Create singleton](#create-singleton)
  - [Other](#other)


## Setup

```bash
mkdir create_test_dataset
cd create_test_dataset
```

## Get bams and vcf for publicly available trio

ChineseTrio

```bash
# vcf
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/analysis/RTG_RTGJointTrio_06062019/GRCh37/family.merged.vcf.gz

# bams
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG005_NA24631_son/NIST_Stanford_Illumina_6kb_matepair/bams/HG005.mate_pair.sorted.bam
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG005_NA24631_son/NIST_Stanford_Illumina_6kb_matepair/bams/HG005.mate_pair.sorted.bam.bai
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG006_NA24694-huCA017E_father/NIST_Stanford_Illumina_6kb_matepair/bams/HG006.mate_pair.sorted.bam
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG006_NA24694-huCA017E_father/NIST_Stanford_Illumina_6kb_matepair/bams/HG006.mate_pair.sorted.bam.bai
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG007_NA24695-hu38168_mother/NIST_Stanford_Illumina_6kb_matepair/bams/HG007.mate_pair.sorted.bam
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG007_NA24695-hu38168_mother/NIST_Stanford_Illumina_6kb_matepair/bams/HG007.mate_pair.sorted.bam.bai
```

See [here](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/analysis/RTG_RTGJointTrio_06062019/GRCh37/README.md) for information on this data

## Manually create pedigree file

Label manually created pedigree file as `NA24694_pedigree.ped` and put in the `.pedigrees/` directory

## Reduce dataset - subset by exome capture regions

Create conda env with GATK4 installed

```bash
conda create -n gatk4 python=3.7.6
conda install -c bioconda gatk4=4.1.9.0
conda activate gatk4
```

Subset vcf by custom exome capture regions

```bash
gatk SelectVariants \
-V family.merged.vcf.gz \
-O family.merged.subset.vcf \
-L custom.bed.gz
```

## Prepare files and file directory structure to run through pipeline

```bash
mkdir -p human_genomics_pipeline/results/called/ pedigrees/
mv family.merged.subset.vcf human_genomics_pipeline/results/called/NA24694_raw_snps_indels.g.vcf
```

## Run through vcf_annotation_pipeline

Get pipeline

```bash
git clone https://github.com/ESR-NZ/vcf_annotation_pipeline.git
cd vcf_annotation_pipeline
git checkout 75774b4ffef3258d0c7a5cff13072971e8553104
```

Configure pipeline with the following files:

- ./vcf_annotation_pipeline/config/cluster.json
- ./vcf_annotation_pipeline/config/config.yaml
- ./vcf_annotation_pipeline/workflow/run_hpc.sh

Run pipeline

```bash
cd workflow

# Create screen to run pipeline in
screen -S create_test_dataset

# Create environment to run pipeline in
conda env create -f pipeline_run_env.yml
conda activate pipeline_run_env

# Dryrun
bash dryrun_hpc.sh

# Run
bash run_hpc.sh
```

## Randomly sub-sample variants

```bash

```

## Create a bed file from vcf

```bash

```

## Pull out fastq reads from bam

```bash

```

## Create singleton

## Other

Provide custom exome capture file and bed file

Random seed for random sampling of variants - re-generating this file might produce different results due to randomisation

Provide pipeline configuration files
