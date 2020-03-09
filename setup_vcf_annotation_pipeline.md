# Set up and run vcf_annotation_pipeline

- **Aim:** Set up and run the [vcf_annotation_pipeline](https://github.com/leahkemp/vcf_annotation_pipeline.git)
- **Prerequisite software:** Conda 4.8.2
- **Prerequisite data** Reference human genome
- **OS:** Ubuntu 16.04 (Wintermute - research server)

The vcf_annotation_pipeline is designed to accept data outputted by the ['human_genomics_pipeline'](https://github.com/ESR-NZ/human_genomics_pipeline). Therefore it assumes that a reference human genome has been downloaded (see the associated documentation [here](https://github.com/leahkemp/documentation/blob/master/setup_human_genomics_pipeline.md)).

## Steps

1. Download data/repository
2. Set up the working environment
3. Run the pipeline

## 1. Download data/repository

### VEP database

Create a conda environment

```bash
conda create --name vep python=3.7
conda activate vep
```

Create a VEP folder in the publicData directory

```bash
cd /store/lkemp
mkdir vep
```

Install conda package of VEP

```bash
conda install -c bioconda ensembl-vep=99.2
```

This installs only the variant effect predictor (VEP) library code. We will need to install data libraries.

Download VEP database

[Download using a command installed with the ensembl-vep conda package](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/ensembl-vep/meta.yaml)

*These are very large files and will likely take some time to download*

```bash
# GRCh37/hg19
vep_install -a cf -s homo_sapiens -y GRCh37 -c /store/lkemp/publicData/vep/GRCh37 --CONVERT
# GRCh38/hg38
vep_install -a cf -s homo_sapiens -y GRCh38 -c /store/lkemp/publicData/vep/GRCh38 --CONVERT
```

### Clone repository

Clone the [vcf_annotation_pipeline](https://github.com/leahkemp/vcf_annotation_pipeline.git) repository

```bash
git clone https://github.com/leahkemp/vcf_annotation_pipeline.git
```

## 2. Set up the working environment

### Set the working directories

### Create a conda environment

```bash
conda create --name annot_pipeline_env python=3.7
conda activate annot_pipeline_env
```

Install snakemake in your conda environment

```bash
conda install --channel bioconda snakemake
```

## 3. Run the pipeline

Start a dry run

```bash
snakemake -n --use-conda
```

If there are no issues, start a full run

```bash
snakemake --use-conda
```
