# How to build VCF-DART from GitHub repositories

- **Aim:** Build and run [VCF-DART](https://www.sciencedirect.com/science/article/abs/pii/S1525157819303538?via%3Dihub) from two github repositories ([here](https://github.com/sirselim/diagnostics_exome_reporting.git) and [here](https://github.com/sirselim/WES_ShinyDiscover.git))
- **Prerequisite software:** [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html)
- **OS:** Ubuntu 16.04 (Wintermute - research server)
- **Date**: 2020-03-11

## Table of contents

- [How to build VCF-DART from GitHub repositories](#how-to-build-vcf-dart-from-github-repositories)
  - [Table of contents](#table-of-contents)
  - [Clone repositories](#clone-repositories)
  - [Install software](#install-software)
    - [USING CONDA (much easier that doing this manually)](#using-conda-much-easier-that-doing-this-manually)
    - [WITHOUT CONDA (MANUAL INSTALLATION)](#without-conda-manual-installation)
    - [snpEFF](#snpeff)
    - [bcftools, tabix and bgzip](#bcftools-tabix-and-bgzip)
    - [parallel](#parallel)
    - [bedops](#bedops)
    - [VEP](#vep)
  - [Install R libraries](#install-r-libraries)
  - [Download data](#download-data)
    - [dbSNP](#dbsnp)
    - [VEP](#vep-1)
    - [dbNSFP](#dbnsfp)
  - [Set file directories](#set-file-directories)
  - [Run RShiny app](#run-rshiny-app)
  - [Further information](#further-information)

## Clone repositories

VCF-DART is composed of two modules publicly available on GitHub

- diagnostics_exome_reporting
- WES_ShinyDiscover

Create a local project folder

```bash
mkdir VCF-DART
cd VCF-DART
```

Clone both repositories into this folder using git bash

```bash
git clone https://github.com/sirselim/diagnostics_exome_reporting.git
git clone https://github.com/sirselim/WES_ShinyDiscover.git
```

Example data can also be cloned

```bash
git clone https://github.com/sirselim/VCF-DART_example_exomes.git
```

## Install software

The following software is required to run VCF-DART successfully:

- [VEP](https://www.ensembl.org/vep) version 92.5
- [snpEFF](snpeff.sourceforge.net/) version 4.3t (for SNPSift dbNSFP annotation)
- [tabix](www.htslib.org/doc/tabix.html) version 1.7-2 (compression and indexing)
- parallel version 20161222
- [bedops](https://bedops.readthedocs.io/) version 2.4.26 (for vcf-sort) (vcf-sort is a function in the package vcftools?)
- [bcftools](https://samtools.github.io/bcftools/bcftools.html) version 1.7-2
- [R](https://www.r-project.org/) version 3.4.4
- [Shiny Server](https://www.rstudio.com/products/shiny/shiny-server/) version 1.5.7.907
- Java (open-jdk 11)

### USING CONDA (much easier that doing this manually)

Create a conda environment

```bash
conda create --name vcfdart_env
```

Use conda to install the software dependencies

```bash
conda install --channel bioconda ensembl-vep # VEP
conda install --channel bioconda snpeff # snpEff
conda install --channel bioconda snpsift # SnpSift
conda install --channel bioconda ensembl-vep # VEP
conda install --channel bioconda bedops # bedops
conda install --channel bioconda bcftools=1.7
conda install --channel bioconda vcftools #vcf-sort ???
```

### WITHOUT CONDA (MANUAL INSTALLATION)

Create a local bin folder to hold all the programs required by VCF-DART

```bash
mkdir bin
```

Add this folder to your path variable so these programs can found and executed by VCF-DART

```bash
export PATH="$home/lkemp/bin:$PATH" # Change 'lkemp' to your folder name
```

### snpEFF

Note that SnpEff and SnpSift are bundled together in the same download

Download and unzip using the command line

```bash
wget https://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip
unzip snpEff_v4_3t_core.zip
```

There shouldn't be any configuration required, however, you may need to change the location of the database.

Download SnpEff databases

```bash
wget
```

See http://snpeff.sourceforge.net/download.html for more information on how to download and install snpEFF

### bcftools, tabix and bgzip

tabix and bgzip is bundled within HTSlib which is bundled within BCFtools. Therefore, you can download bcftools and install HTSlib from it (which provides the tabix and bgzip utility)

Download bcftools and unzip

```bash
wget https://sourceforge.net/projects/samtools/files/samtools/1.7/bcftools-1.7.tar.bz2
tar xjf bcftools-1.7.tar.bz2
```

Install HTSlib

```bash
cd bcftools-1.7
cd htslib-1.7
./configure --prefix=/home/lkemp # Will place it in the bin subdirectory we have already created
make
make install
```

Install bcftools

```bash
cd ..
./configure --prefix=/home/lkemp
make
make install
```

See [htslib](http://www.htslib.org/download/) for more information on how to download and install bcftools/HTSlib/tabix/bgzip

### parallel

Download and unzip

```bash
cd ..
wget http://ftp.gnu.org/gnu/parallel/parallel-20160222.tar.bz2
tar xjf parallel-20160222.tar.bz2
```

Install

```bash
cd parallel-20160222
./configure --prefix=/home/lkemp
make
make install
```

### bedops

Download and unzip

```bash
cd ..
wget https://github.com/bedops/bedops/releases/download/v2.4.26/bedops_linux_x86_64-v2.4.26.tar.bz2
tar jxvf bedops_linux_x86_64-v2.4.26.tar.bz2
```

Copy the extracted binaries to bin folder

```bash
cp bin/bedops /home/lkemp/bin
```

NOTE: BEDOPS DOESN'T SEEM TO CONTAIN VCF-SORT AS EXPECTED/STATED BY PAPER

Therefore I will install VCFtools

```bash
git clone https://github.com/vcftools/vcftools.git
```

Install

```bash
cd vcftools
./autogen.sh
./configure --prefix=/home/lkemp
make
make install
```

[VCFtools](https://vcftools.github.io/perl_module.html#vcf-sort)

### VEP

Note: *Apparently this software is very difficult to install due to dependencies. Apparently there is also a containerized verison of this software available. Miles hasn't used the containerized version when he set up VCF-DART, but in theory it's a good idea to use*

Clone the git repository and checkout an branch containing the right verison of the software

```bash
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
git checkout release/92.5
perl INSTALL.pl
```

- - -
The VEP package requires Perl (>=5.10 recommended, tested on 5.8, 5.10, 5.14, 5.18, 5.22) and the DBI package installed. The remaining dependencies can be installed using the included INSTALL.pl script. The installer may also be used to check for updates and co-dependent packages, simply re-run INSTALL.pl.

```bash
perl INSTALL.pl
```

- - -

*My example workflow*

- I was missing a perl module, so the following error message was returned when I ran the installation

```bash
ERROR: DBI module not found. VEP requires the DBI perl module to function
```

- So I installed the missing module:

```bash
# Check that perl and cpan are installed
which perl
which cpan
# Intall DBI perl module
cpan App::cpanminus
cpanm DBI
```

(Click [here](http://www.cpan.org/modules/INSTALL.html) for more information on installing perl modules)

- And re-run the installation of VEP

```bash
perl INSTALL.pl
```

- Another error

```bash
WARNING: DBD::mysql module not found. VEP can only run in offline (--offline) mode without DBD::mysql installed
```

```bash
perl INSTALL.pl --NO_HTSLIB
```

- My response to the prompts:

```bash
Do you wish to exit so you can get updates (y) or continue (n):
n
```

```bash
Destination directory ./Bio already exists. Do you want to overwrite it (if updating VEP this is probably OK) (y/n)?
y
```

```bash
The VEP can either connect to remote or local databases, or use local cache files.
Using local cache files is the fastest and most efficient way to run the VEP
Cache files will be stored in /home/lkemp/.vep
Do you want to install any cache files (y/n)?
n
```

```bash
The VEP can use FASTA files to retrieve sequence data for HGVS notations and reference sequence checks.
FASTA files will be stored in /home/lkemp/.vep
n
```

```bash
The VEP can use plugins to add functionality and data.
Plugins will be installed in /home/lkemp/.vep/Plugins
Do you want to install any plugins (y/n)?
n
```

Copy software to 'bin'

```bash
cd ..
cp -R /home/lkemp/scratch/ensembl-vep /home/lkemp/bin
```

(See https://asia.ensembl.org/info/docs/tools/vep/script/vep_download.html for more information on how to download and install VEP)

## Install R libraries

The following R packages are required to run VCF-DART successfully:

- magrittr version 1.5
- shiny version 1.1.0
- rmarkdown version 1.10
- pander version 0.62
- DT version 0.4

Within R/RStudio, install the required packages

```r
install.packages("https://cran.r-project.org/src/contrib/magrittr_1.5.tar.gz", repos=NULL, type="source")
install.packages("https://cran.r-project.org/src/contrib/Archive/shiny/shiny_1.1.0.tar.gz", repos=NULL, type="source")
install.packages("https://cran.r-project.org/src/contrib/Archive/rmarkdown/rmarkdown_1.10.tar.gz", repos=NULL, type="source")
install.packages("https://cran.r-project.org/src/contrib/Archive/pander/pander_0.6.2.tar.gz", repos=NULL, type="source")
install.packages("https://cran.r-project.org/src/contrib/Archive/DT/DT_0.4.tar.gz", repos=NULL, type="source")
```

Note. you may need to download the package 'crosstalk' to ensure the final package (DT) installs correctly

```r
install.packages("crosstalk")
# Then repeat the installation of the DT package
install.packages("https://cran.r-project.org/src/contrib/Archive/DT/DT_0.4.tar.gz", repos=NULL, type="source")
```

Done!

## Download data

### dbSNP

[Single nucleotide polymorphism database](https://www.ncbi.nlm.nih.gov/snp/)

```bash
# GRCh37/hg19
wget ftp://ftp.ncbi.nih.gov:21/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz

# GRCh38/hg38
wget ftp://ftp.ncbi.nih.gov:21/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz
```

*See the file "setup_human_genomics_pipeline" for a more detailed description of this step*

### VEP

[Variant effect predictor database]()

### dbNSFP

## Set file directories

Change the file paths of all the software, databases and raw data in diagnostics_exome_reporting/WESdiag_pipeline_dev.sh

## Run RShiny app

- Open RStudio
- Open server.R from diagnostics_exome_reporting
- Run app
- Input data the app asks for

Note. in order to upload VCF files in the RShiny app, you will need to be using the server directly instead of connecting to the server via SSH. This is because the app searches for vcf data to upload from your local computer, but it can't because the files are too large.

This will output the file "pipeline_input.txt" that is called upon by WESdiag_pipeline_dev.sh

## Further information

The journal article associated with VCF-DART can be found [here](https://www.sciencedirect.com/science/article/abs/pii/S1525157819303538?via%3Dihub)
