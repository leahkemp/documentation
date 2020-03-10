# Set up and run vcf_annotation_pipeline

- **Aim:** Set up and run the [vcf_annotation_pipeline](https://github.com/leahkemp/vcf_annotation_pipeline.git)
- **Prerequisite software:** [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html), [tabix](http://www.htslib.org/doc/tabix.html), [bgzip](http://www.htslib.org/doc/bgzip.html), [gunzip](https://linux.die.net/man/1/gunzip)
- **Prerequisite data** Reference human genome and dbSNP database
- **OS:** Ubuntu 16.04 (Wintermute - research server)
- **Date**: 11/03/2020

The vcf_annotation_pipeline is designed to accept data outputted by the ['human_genomics_pipeline'](https://github.com/ESR-NZ/human_genomics_pipeline). Therefore it assumes that a reference human genome and dbSNP database has been downloaded. The ["Setup and run human_genomics_pipeline"](https://github.com/leahkemp/documentation/blob/master/setup_human_genomics_pipeline.md) documentation provides instructions on how to download this data.

## Table of contents

* [Download data/repository](#download-data/repository)
* [Set up the working environment](#set-up-the-working-environment)
* [Run the pipeline](#run-the-pipeline)

## Download data/repository

### VEP database

Create a conda environment

```bash
conda create --name download_data_env python=3.7
conda activate download_data_env
```

Install conda package of VEP and gatk4 (will use later)

```bash
conda install -c bioconda ensembl-vep=99.2
conda install -c bioconda gatk4=4.1.5.0
```

Download VEP database with with the [ensembl-vep conda package](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/ensembl-vep/meta.yaml)

*These are very large files and will likely take some time to download*

```bash
# GRCh37/hg19
vep_install -a cf -s homo_sapiens -y GRCh37 -c /store/lkemp/publicData/vep/GRCh37 --CONVERT
# GRCh38/hg38
vep_install -a cf -s homo_sapiens -y GRCh38 -c /store/lkemp/publicData/vep/GRCh38 --CONVERT
```

### Other databases

GRCh37

dbNSFP database was custom built by Miles?

```bash
# Mills
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz

# 1000G indel
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.gz

# 1000G snp
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz

# Omni
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/1000G_omni2.5.hg19.sites.vcf.gz

# Hapmap3
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/hapmap_3.3.hg19.sites.vcf.gz

# Axiom exome plus
## Not available for GRCh37?
```

If an index file (.tbi) is not available, one will need to be created using tabix (or [IndexFeatureFile by gatk](https://gatk.broadinstitute.org/hc/en-us/articles/360036899892-IndexFeatureFile)). To do so, the files will need to be bgzipped. To check if they are the correct format, run:

```bash
htsfile Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
```

If the file is not bgzipped, you will need to unzip it, then rezip it to bgzip format

```bash
gunzip Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
bgzip Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
```

Then the files can be indexed with tabix

```bash
tabix Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
```

Repeat with the other files

```bash
gunzip 1000G_phase1.indels.hg19.sites.vcf.gz
bgzip 1000G_phase1.indels.hg19.sites.vcf
tabix 1000G_phase1.indels.hg19.sites.vcf.gz

gunzip 1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
bgzip 1000G_phase1.snps.high_confidence.hg19.sites.vcf
tabix 1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz

gunzip 1000G_omni2.5.hg19.sites.vcf.gz
bgzip 1000G_omni2.5.hg19.sites.vcf
tabix 1000G_omni2.5.hg19.sites.vcf.gz

gunzip hapmap_3.3.hg19.sites.vcf.gz
bgzip hapmap_3.3.hg19.sites.vcf
tabix hapmap_3.3.hg19.sites.vcf.gz
```

GRCh38

Get dbNSFP database for GRCh38?

```bash
# Mills
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

# 1000G indel

### Not available?

# 1000G snp
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi

# Omni
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/1000G_omni2.5.hg38.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/1000G_omni2.5.hg38.vcf.gz.tbi

# Hapmap3
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/hapmap_3.3.hg38.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/hapmap_3.3.hg38.vcf.gz.tbi

# Axiom exome plus
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi
```

Since the index files (.tbi) are available for download for GRCh38, we don't need to create them with tabix

### Clone repository

Clone the [vcf_annotation_pipeline](https://github.com/leahkemp/vcf_annotation_pipeline.git) repository

```bash
git clone https://github.com/leahkemp/vcf_annotation_pipeline.git
```

## Set up the working environment

### Set the working directories

Set the working directories in the first section of the Snakefile so that the pipeline can access the vcf data outputted by 'human_genomics_pipeline' (or other vcf data to be annotated) as well as the databases we downloaded (eg. Mills, 1000G snp, Omni etc.)

### Create a conda environment

```bash
conda create --name annot_pipeline_env python=3.7
conda activate annot_pipeline_env
```

Install snakemake in your conda environment

```bash
conda install --channel bioconda snakemake
```

## Run the pipeline

Start a dry run

```bash
snakemake -n --use-conda
```

If there are no issues, start a full run

```bash
snakemake --use-conda
```
