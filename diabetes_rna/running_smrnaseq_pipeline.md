# Running smrnaseq pipeline

Created: 2020/11/13 12:29:25
Last modified: 2020/11/13 12:35:25

- **Aim:** In [this document](./rna_pipelines_current_status.md) I settled on using the [smrnaseq](https://github.com/nf-core/smrnaseq) nextflow pipeline to process our small non-coding RNA-seq data. This document documents/describes the process of trying this pipeline out on the data (extending on Miles Bentons work)
- **Prerequisite software:** [conda 4.9.0](https://docs.conda.io/en/latest/), [git 2.7.4](https://git-scm.com/)
- **OS:** Ubuntu 16.04 (Wintermute - research server)

## Table of contents

- [Running smrnaseq pipeline](#running-smrnaseq-pipeline)
  - [Table of contents](#table-of-contents)
  - [Clone pipeline](#clone-pipeline)
  - [Install nextflow into pipeline run env](#install-nextflow-into-pipeline-run-env)
  - [Run the pipeline](#run-the-pipeline)

## Clone pipeline

```bash
cd /store/lkemp/smrnaseq_hps/
git clone https://github.com/nf-core/smrnaseq.git
```

## Install nextflow into pipeline run env

```bash
conda create -n nextflow python=3.7.6
conda activate nextflow
mamba install -c bioconda nextflow=20.07.1
```

## Run the pipeline

```bash
cd /store/lkemp/smrnaseq_hps/

nextflow run /store/lkemp/smrnaseq_hps/smrnaseq/main.nf --reads 'fastq/*_combined.fastq.gz' \
  -profile conda --protocol illumina --genome 'GRCh37' \
  --saveReference -resume --min_length 17
```

Worked great!