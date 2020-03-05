# Set up and run human_genomics_pipeline

**Aim:** Set up and run the [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline)
**Prerequisite software:** Conda 4.8.2
**OS:** Ubuntu 16.04 (Wintermute - research server)

## Steps

[1. Download data/repository](#step-1)
[2. Set up the working environment](#step-2)
[3. Run the pipeline](#step-3)

## 1. Download data/repository {: #step-1 }

### Clone repository

Clone the [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) git repository

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

---

**BaseRecalibrator** requires known SNPs and indels passed with the -knownSites argument to function properly. We use all the following files:

- The most recent dbSNP release (build ID > 132)
- Mills_and_1000G_gold_standard.indels.b37.vcf
- 1000G_phase1.indels.b37.vcf (currently from the 1000 Genomes Phase I indel calls)

**HaplotypeCaller**
These tools do NOT require known sites, but if SNPs are provided with the -dbsnp argument they will use them for variant annotation. We use this file:

- The most recent dbSNP release (build ID > 132)

---

#### Option one: download from the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360036212652-Resource-Bundle) (recommended)

*These are large files and make take some time to download*

```bash
# GRCh37/hg19
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/dbsnp_138.hg19.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/dbsnp_138.hg19.vcf.idx.gz # Might not need?
# Additional data required for BaseRecalibrator
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.gz

# GRCh38/hg38
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/dbsnp_146.hg38.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/dbsnp_146.hg38.vcf.gz.tbi
# Additional data required for BaseRecalibrator
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
# No indel file available?
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
# GRCh38/hg38 (already downloaded the .tbi file, don;t need to create one)
```

#### Option two: download from [NCBI](https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/)

*These are large files and make take some time to download*

Information on dbSNP files can be found on the [NCBI website](https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/).

- Download the [appropriate](https://bioinformatics.stackexchange.com/questions/4578/how-to-download-dbsnp-database) dbSNP database. b151 is the current newest version of the database.

```bash
# GRCh37/hg19
wget ftp://ftp.ncbi.nih.gov:21/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz
wget ftp://ftp.ncbi.nih.gov:21/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz.tbi

# GRCh38/hg38
wget ftp://ftp.ncbi.nih.gov:21/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz
wget ftp://ftp.ncbi.nih.gov:21/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz.tbi
```

TODO: find out if you need to create the additional files (Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz and 1000G_phase1.indels.hg19.sites.vcf.gz if you download the dbsnp database from NCBI)

### Example WGS data

- Download raw WGS data from [NIST genome in a bottle](https://www.nist.gov/programs-projects/genome-bottle).

```bash
wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/140528_D00360_0018_AH8VC6ADXX/Project_RM8391_RM8392/Sample_2A1/2A1_CGATGT_L001_R1_001.fastq.gz
wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/140528_D00360_0018_AH8VC6ADXX/Project_RM8391_RM8392/Sample_2A1/2A1_CGATGT_L001_R2_001.fastq.gz
```

This downloads WGS data for the mother of the AshkenazimTrio ([sample HG004](https://github.com/genome-in-a-bottle/giab_data_indexes))

## 2. Set up the working environment {: #step-2 }

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

## 3. Run the pipeline  {: #step-3 }

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
# check multiqc is installed
which multiqc
# install if necessary
conda install channel --bioconda multiqc
multiqc .
```
