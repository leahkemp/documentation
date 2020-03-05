# Install VEP locally

This document describes how I installed VEP locally. VEP is difficult to install locally ([instructions here](https://asia.ensembl.org/info/docs/tools/vep/script/vep_download.html)) largely due to perl dependencies. There is a conda package available for VEP [ensembl-vep](https://anaconda.org/bioconda/ensembl-vep) available through bioconda, yet we haven't yet tested it's utility alone, or within a snakemake workflow. Here I will attempt to run the ensembl-vep package locally before attempting to incorporate it into a snakemake workflow.

*Note - there is also a [docker container available for VEP](https://hub.docker.com/r/ensemblorg/ensembl-vep/)*

- **Aim:** Install VEP locally (using the conda ensembl-vep package) to ensure it works before integrating it into a snakemake workflow
- **Prerequisite software:** Conda 4.8.2
- **OS:** Ubuntu 16.04 (Wintermute - research server)

## Create environment

Create a conda environment

```bash
conda create --name vepLocalTest python=3.8
conda activate vepLocalTest
```

Create a VEP folder in the publicData directory

```bash
cd /store/lkemp
mkdir vep
```

## Install VEP

Install conda package of VEP

```bash
conda install -c bioconda ensembl-vep
```

This installed ensembl-vep version 92.4

## Download VEP database

See instructions [here](http://asia.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache)

*Make sure you install the VEP cache data that is the same version as the VEP you downloaded*

```bash
cd /store/lkemp/publicData/vep/
# GRCh37
wget ftp://ftp.ensembl.org:21/pub/release-92/variation/VEP/homo_sapiens_merged_vep_92_GRCh37.tar.gz
# GRCh38
wget ftp://ftp.ensembl.org:21/pub/release-92/variation/VEP/homo_sapiens_merged_vep_92_GRCh38.tar.gz
```
