# Set up and run human_genomics_pipeline

**Aim:** Set up and run the [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline)
**Prerequisite software:** Conda 4.8.2
**OS:** Ubuntu 16.04 (Wintermute - research server)

## Download data/repository

### Clone repository

Clone the human_genomics_pipeline git repository

```bash
git clone https://github.com/ESR-NZ/human_genomics_pipeline.git
```

### Reference human genome

There are many places you can download the reference human genome (and many ways to download it).

TODO: add info on how I initially produced the reference human genome

Make dict file using picard (within GATK)

```bash
java -jar /store/mbenton/software/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar CreateSequenceDictionary -R /store/lkemp/publicData/referenceGenome/GRCh37/hg19.fa
```

Make index file using samtools

```bash
samtools faidx /store/lkemp/publicData/referenceGenome/GRCh37/hg19.fa
```

We will download it from the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360036212652-Resource-Bundle) over the ftp server.

*Downloading the fasta files from the GATK bundle allows us to also download their associated fasta sequence dictionary file (.dict) and fasta index file (.fai) files that would otherwise [need to be created](https://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference) for GATK to use in the pipeline*

Download and unzip

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

Index the reference human genomes

```bash
bwa index -a bwtsw ucsc.hg19.fasta
bwa index -a bwtsw Homo_sapiens_assembly38.fasta
```

### dbSNP database

Information on dbSNP files can be found on the [NCBI website](https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/). You can also download small index files of the ftp site using the command line so that you can see what is available for download. For example:

```bash
wget ftp://ftp.ncbi.nih.gov:21/snp/organisms/
wget ftp://ftp.ncbi.nih.gov:21/snp/organisms/human_9606_b151_GRCh37p13/VCF/
wget ftp://ftp.ncbi.nih.gov:21/snp/organisms/human_9606_b151_GRCh38p7/VCF/
```

Then download the [appropriate](https://bioinformatics.stackexchange.com/questions/4578/how-to-download-dbsnp-database) dbSNP database. b151 is the current newest version of the database.

```bash
# GRCh37/hg19
wget ftp://ftp.ncbi.nih.gov:21/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz
wget ftp://ftp.ncbi.nih.gov:21/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz.tbi

# GRCh38/hg38
wget ftp://ftp.ncbi.nih.gov:21/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz
wget ftp://ftp.ncbi.nih.gov:21/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz.tbi
```

Note: *these are large files and make take some time to download*

### Example WGS data

Getting example WGS raw data from [NIST genome in a bottle](https://www.nist.gov/programs-projects/genome-bottle).

```bash
wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/140528_D00360_0018_AH8VC6ADXX/Project_RM8391_RM8392/Sample_2A1/2A1_CGATGT_L001_R1_001.fastq.gz
wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/140528_D00360_0018_AH8VC6ADXX/Project_RM8391_RM8392/Sample_2A1/2A1_CGATGT_L001_R2_001.fastq.gz
```

This downloads WGS data for the mother of the AshkenazimTrio ([sample HG004](https://github.com/genome-in-a-bottle/giab_data_indexes))

## Set up the working environment

### Set the working directories

Set the working directories of the human_genomics_pipeline by manually editing the first section of 'Snakefile' and 'Merge_QC.snakemake'.

Ensure that your fastq files are within a fastq directory within your working directory. eg. ../yourWorkingDirectory/fastq/4001_1.fastq.gz

### Create a conda environment

Make sure conda is installed

```bash
which conda
```

If not, install anaconda, directions for doing this on Ubuntu 16.04 can be found [here](https://www.digitalocean.com/community/tutorials/how-to-install-the-anaconda-python-distribution-on-ubuntu-16-04)

Update conda if necessary

```bash
conda update -n base -c defaults conda
```

Add conda channels (note. the order you run these commands is important)

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Create a conda environment including python, then activate it

```bash
cd /home/lkemp/human_genomics_pipeline
conda create --name pipeline_env python=3.7
activate pipeline_env
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
# check mutliqc is installed
which multiqc
# install if neccesary
conda install channel --bioconda multiqc
multiqc .
```
