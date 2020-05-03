# Install VEP locally

Created: 2020/03/11 11:25:43
Last modified: 2020/03/11 11:26:45

- **Aim:** Install [VEP](https://asia.ensembl.org/info/docs/tools/vep/index.html) locally (using the [conda ensembl-vep package](https://anaconda.org/bioconda/ensembl-vep)) to ensure it works before integrating it into a snakemake workflow
- **Prerequisite software:** [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html)
- **OS:** Ubuntu 16.04 (Wintermute - research server)

This document describes how I installed VEP locally. VEP is difficult to install locally ([instructions here](https://asia.ensembl.org/info/docs/tools/vep/script/vep_download.html)) largely due to perl dependencies. There is a conda package available for VEP [ensembl-vep](https://anaconda.org/bioconda/ensembl-vep) available through bioconda, yet we haven't yet tested it's utility alone, or within a snakemake workflow. Here I will attempt to run the ensembl-vep package locally before attempting to incorporate it into a snakemake workflow.

*Note - there is also a [docker container available for VEP](https://hub.docker.com/r/ensemblorg/ensembl-vep/)*

## Table of contents

- [Install VEP locally](#install-vep-locally)
  - [Table of contents](#table-of-contents)
  - [Create environment](#create-environment)
  - [Install VEP](#install-vep)
  - [Download VEP database](#download-vep-database)
  - [Run VEP](#run-vep)

## Create environment

Create a conda environment

```bash
conda create --name vep python=3.7 # note. there were dependency conflicts when I used python version 3.8
conda activate vep
```

Create a VEP folder in the publicData directory

```bash
cd /store/lkemp
mkdir vep
```

## Install VEP

Install conda package of VEP

```bash
conda install -c bioconda ensembl-vep=99.2
```

This installs only the variant effect predictor (VEP) library code. We will need to install data libraries.

## Download VEP database

[Download using a command installed with the ensembl-vep conda package](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/ensembl-vep/meta.yaml)

*These are very large files and will likely take some time to download*

```bash
# GRCh37/hg19
vep_install -a cf -s homo_sapiens -y GRCh37 -c /store/lkemp/publicData/vep/GRCh37 --CONVERT
# GRCh38/hg38
vep_install -a cf -s homo_sapiens -y GRCh38 -c /store/lkemp/publicData/vep/GRCh38 --CONVERT
```

## Run VEP

See the 'cache options' section cache options at [emsembl-vep](http://asia.ensembl.org/info/docs/tools/vep/script/vep_options.html)

```bash
vep --assembly GRCh37 --cache --dir /store/lkemp/publicData/vep/GRCh37/ --fasta /store/lkemp/publicData/referenceGenome/hg19/ucsc.hg19.fasta -i /store/lkemp/exome_project/human_genomics_pipeline_official/vcf/MS_16BL0795_S5.raw.snps.indels.AS.g.vcf --stats_text --everything -o /store/lkemp/exome_project/human_genomics_pipeline_official/vcf/MS_16BL0795_S5.raw.snps.indels.AS.g.VEP.vcf --vcf --force_overwrite --offline
```
