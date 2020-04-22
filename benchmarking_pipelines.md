# Benchmarking genomic pipelines

Created: 2020-04-22 13:37:04
Last modified: 2020/04/22 16:22:14

- **Aim:** Undertake benchmarking of genomics pipelines to test their quality for clinical use
- **Prerequisite software:** [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html)
- **OS:** Ubuntu 16.04 (Wintermute - research server)

## Table of contents

- [Benchmarking genomic pipelines](#benchmarking-genomic-pipelines)
  - [Table of contents](#table-of-contents)
  - [Setup](#setup)
    - [Create conda environment](#create-conda-environment)
    - [Install software](#install-software)
    - [Download data](#download-data)
  - [Benchmarking](#benchmarking)
    - [human_genomics_pipeline](#humangenomicspipeline)
    - [Nvidia parabricks](#nvidia-parabricks)

## Setup

### Create conda environment

```bash
conda create -n benchmarking_env python=2.7.5
conda activate benchmarking_env
```

(Need python 2 in order to get htslib2 that is required for installation)

### Install software

Install [hap.py](https://github.com/Illumina/hap.py) dependencies

```bash
conda install -c conda-forge cython=0.29.15
conda install -c conda-forge scipy=1.2.1 # Also installs numpy=1.16.5 dependency
conda install -c conda-forge pandas=0.24.2
conda install -c bioconda pybedtools=0.8.1 # Also installs pysam=0.15.3 dependency
conda install -c bioconda bx-python=0.8.8
```

Install [hap.py](https://github.com/Illumina/hap.py) (using the python and the helper script)

```bash
git clone git@github.com:Illumina/hap.py.git
cd hap.py
python install.py ~/hap.py-install --no-tests
```

### Download data

## Benchmarking

### [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline)

With [hap.py](https://github.com/Illumina/hap.py)

```bash
python /home/lkemp/hap.py-install/bin/hap.py /store/lkemp/exome_project/benchmarking/NA12878_exome/known/project.NIST.hc.snps.indels.vcf /store/lkemp/exome_project/benchmarking/NA12878_exome/human_genomics_pipeline/vcf/NA12878_NIST.raw.snps.indels.AS.g.vcf -f /store/lkemp/exome_project/benchmarking/NA12878_exome/known/nexterarapidcapture_expandedexome_targetedregions.bed.gz -o benchmarking -r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh38/Homo_sapiens_assembly38.fasta
```

### [Nvidia parabricks](https://github.com/ESR-NZ/ESR-Parabricks)