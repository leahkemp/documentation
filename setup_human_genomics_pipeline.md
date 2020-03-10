# Set up and run human_genomics_pipeline

- **Aim:** Set up and run the [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline)
- **Prerequisite software:**  [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html), [tabix](http://www.htslib.org/doc/tabix.html), [bgzip](http://www.htslib.org/doc/bgzip.html), [gunzip](https://linux.die.net/man/1/gunzip), [bwa](http://bio-bwa.sourceforge.net/), [samtools](http://www.htslib.org/), [gatk](https://gatk.broadinstitute.org/hc/en-us)
- **OS:** Ubuntu 16.04 (Wintermute - research server)
- **Date:** 2020-03-11

## Table of contents

- [Set up and run human_genomics_pipeline](#set-up-and-run-humangenomicspipeline)
  - [Table of contents](#table-of-contents)
  - [Download data/repository](#download-datarepository)
    - [Clone repository](#clone-repository)
    - [Reference human genome](#reference-human-genome)
      - [Option one: download from the GATK resource bundle (recommended)](#option-one-download-from-the-gatk-resource-bundle-recommended)
      - [Option two: download from UCSC](#option-two-download-from-ucsc)
    - [dbSNP database](#dbsnp-database)
      - [Option one: download from NCBI (recommended)](#option-one-download-from-ncbi-recommended)
      - [Option two: download from the GATK resource bundle](#option-two-download-from-the-gatk-resource-bundle)
    - [Example WGS data](#example-wgs-data)
  - [Set up the working environment](#set-up-the-working-environment)
    - [Set the working directories](#set-the-working-directories)
    - [Create a conda environment](#create-a-conda-environment)
  - [Run the pipeline](#run-the-pipeline)

## Download data/repository

Note. it's a good idea to use data such as reference human genome, it's associated files and dbSNP database downloaded from the same source since different sources may label data differently (eg. [chromosome labeling and length](https://gatkforums.broadinstitute.org/gatk/discussion/11359/input-files-reference-and-features-have-incompatible-contigs)

### Clone repository

Clone the [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) repository

```bash
git clone https://github.com/ESR-NZ/human_genomics_pipeline.git
```

### Reference human genome

There are many places you can download the reference human genome (and many ways to download it). Here I will describe two methods I used to download and prepare two releases of reference human genome.

#### Option one: download from the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360036212652-Resource-Bundle) (recommended)

Downloading the reference human genome from the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360036212652-Resource-Bundle) allows us to also download their associated fasta sequence dictionary file (.dict) and fasta index file (.fai) files that would otherwise [need to be created](https://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference)

- Download over the ftp server and unzip files

```bash
# GRCh37/hg19
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.fasta.gz
gunzip ucsc.hg19.fasta.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.dict.gz
gunzip ucsc.hg19.dict.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.fasta.fai.gz
gunzip ucsc.hg19.fasta.fai.gz

# GRCh38/hg38
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.gz
gunzip Homo_sapiens_assembly38.fasta.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.dict
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.fai
```

- Create index files for the genome sequence (.amb, .ann, .bwt, .pac, .sa)

```bash
# GRCh37/hg19
bwa index -a bwtsw ucsc.hg19.fasta
# GRCh38/hg38
bwa index -a bwtsw Homo_sapiens_assembly38.fasta
```

['bwtsw' is required](http://seqanswers.com/forums/showthread.php?t=3547) so that bwa uses the correct algorithm to handle a large whole genome sequence

#### Option two: download from [UCSC](https://hgdownload.soe.ucsc.edu/downloads.html)

Downloading the reference human genome from [UCSC](https://hgdownload.soe.ucsc.edu/downloads.html) will provide only the fasta file with the genome sequence. We will [need to create]((https://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference)) it's associated fasta sequence dictionary file (.dict) and fasta index file (.fai) files.

- Download over the ftp server and unzip files

```bash
# GRCh37/hg19
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
# GRCh38/hg38
wget ftp://hgdownload.soe.ucsc.edu:21/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
```

- Create index files for the genome sequence (.amb, .ann, .bwt, .pac, .sa)

```bash
# GRCh37/hg19
bwa index -a bwtsw hg19.fa
# GRCh38/hg38
bwa index -a bwtsw hg38.fa
```

['bwtsw' is required](http://seqanswers.com/forums/showthread.php?t=3547) so that bwa uses the correct algorithm to handle a large whole genome sequence

- Make the fasta sequence dictionary file (.dict) using picard within GATK

```bash
# GRCh37/hg19
java -jar /store/mbenton/software/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar CreateSequenceDictionary -R /store/lkemp/publicData/referenceGenome/GRCh37/hg19.fa
# GRCh38/hg38
java -jar /store/mbenton/software/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar CreateSequenceDictionary -R /store/lkemp/publicData/referenceGenome/GRCh38/hg38.fa
```

- Make the fasta index file (.fai) with samtools

```bash
# GRCh37/hg19
samtools faidx /store/lkemp/publicData/referenceGenome/GRCh37/hg19.fa
# GRCh38/hg38
samtools faidx /store/lkemp/publicData/referenceGenome/GRCh38/hg38.fa
```

### dbSNP database

See [here](https://gatkforums.broadinstitute.org/gatk/discussion/1247/what-should-i-use-as-known-variants-sites-for-running-tool-x) for information on what files you will need for GATK portions of the pipeline. Below is a summary of the files needed for two GATK functions called by the pipeline:

#### Option one: download from [NCBI](https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/) (recommended)

This option is recommended since the latest release of the dbSNP database was available on this site (current latest release: build 153)

*These are large files and make take some time to download*

Information on dbSNP files can be found on the [NCBI website](https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/) and how to download them is described [here](https://bioinformatics.stackexchange.com/questions/4578/how-to-download-dbsnp-database).

- Download the latest dbSNP database

```bash
# Build 153
# GRCh37/hg19
wget ftp://ftp.ncbi.nlm.nih.gov:21/snp/latest_release/VCF/GCF_000001405.25.gz
wget ftp://ftp.ncbi.nlm.nih.gov:21/snp/latest_release/VCF/GCF_000001405.25.gz.tbi

# GRCh38/hg38
wget ftp://ftp.ncbi.nlm.nih.gov:21/snp/latest_release/VCF/GCF_000001405.38.gz
wget ftp://ftp.ncbi.nlm.nih.gov:21/snp/latest_release/VCF/GCF_000001405.38.gz.tbi
```

#### Option two: download from the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360036212652-Resource-Bundle)

*These are large files and make take some time to download*

These are old releases of the dbSNP database so downloading these is not recommended

```bash
# GRCh37/hg19
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/dbsnp_138.hg19.vcf.gz

# GRCh38/hg38
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/dbsnp_146.hg38.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/dbsnp_146.hg38.vcf.gz.tbi
```

- Create tabix file (.tbi) if not downloaded. In order to do this, our dbsnp vcf file [needs to be bgzf-compressed file](https://github.com/samtools/bcftools/issues/668) (also see [here](https://www.biostars.org/p/138514/)). Check the format of your dbsnp file with:

```bash
htsfile dbsnp_138.hg19.vcf.gz
```

Output should say something like:

```bash
dbsnp_138.hg19.vcf.gz:    VCF version 4.0 BGZf-compressed variant calling data
```

If it doesn't (for example it says you have a gzip-compressed file), you can convert it to a bgzf-compressed file by unzipping it, then recompressing it with bgzip

```bash
gunzip dbsnp_138.hg19.vcf.gz
bgzip dbsnp_138.hg19.vcf.gz
```

Once you have the correct file format, you can create a tabix file (.tbi) with tabix

```bash
# GRCh37/hg19
tabix dbsnp_138.hg19.vcf.gz
```

### Example WGS data

- Download raw WGS data from [NIST genome in a bottle](https://www.nist.gov/programs-projects/genome-bottle).

```bash
wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/140528_D00360_0018_AH8VC6ADXX/Project_RM8391_RM8392/Sample_2A1/2A1_CGATGT_L001_R1_001.fastq.gz
wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/140528_D00360_0018_AH8VC6ADXX/Project_RM8391_RM8392/Sample_2A1/2A1_CGATGT_L001_R2_001.fastq.gz
```

This downloads WGS data for the mother of the AshkenazimTrio ([sample HG004](https://github.com/genome-in-a-bottle/giab_data_indexes))

## Set up the working environment

### Set the working directories

Set the working directories of the human_genomics_pipeline by manually editing the first section of 'Snakefile' and 'Merge_QC.snakemake'. Ensure that the pipeline can find the:

- reference human genome
- dbSNP database
- WGS or WES data

Also, make sure that the global wildcard function (page 23 of the snakefile) that finds your sample names will find/capture your sample name. For example, I needed to change this...

```python
SAMPLES, = glob_wildcards("fastq/{sample}_1.fastq.gz")
```

...to this...

```python
SAMPLES, = glob_wildcards("..fastq/{sample}_R1.fastq.gz")
```

...in order for it to find my WES data labelled 'CH_13BL2450_S1_R1.fastq.gz' and 'CH_13BL2450_S1_R2.fastq.gz'

### Create a conda environment

Create a conda environment including python, then activate it

```bash
conda create --name pipeline_env python=3.7
conda activate pipeline_env
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

To run the multiqc step, direct to the qc/fastqc directory created by the pipeline and run:

```bash
# install multiqc and run on files in the current directory
conda install channel --bioconda multiqc
multiqc .
```
