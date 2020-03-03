# Install VEP locally

This document describes how I installed VEP locally on wintermute (ubuntu 16.04). VEP is difficult to install locally ([instructions here](https://asia.ensembl.org/info/docs/tools/vep/script/vep_download.html)) largely due to perl dependancies. There is a conda package available for vep [ensembl-vep](https://anaconda.org/bioconda/ensembl-vep) available through bioconda, yet we haven't yet tested it's utility within a snakemake workflow. Here I will attempt to run the ensembl-vep package locally before attempting to incorperate it into a snakemake workflow.

**Aim:** Install VEP locally (using the conda emsembl-vep package) to ensure it works before integrating it into a snakemake workflow.
**Prerequisite software:** Conda 4.8.2
**OS:** Ubuntu 16.04 (Wintermute - research server)

## Create environment

Create a conda environment

```bash
conda create --name vepLocalTest python=3.8
conda activate vepLocalTest
```

Create a vep folder in the publicData directory

```bash
cd /store/lkemp
mkdir vep
```