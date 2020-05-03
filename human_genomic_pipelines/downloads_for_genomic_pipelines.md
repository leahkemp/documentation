# Database downloads for human_genomics_pipeline and vcf_annotation_pipeline

Created: 2020/03/11 11:25:43
Last modified: 2020/04/23 15:20:02

- **Aim:** Document where data/databases required for [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) and [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline) can be downloaded and what data processing they require before they can be used
- **Prerequisite software:**  [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html), [tabix](http://www.htslib.org/doc/tabix.html), [bgzip](http://www.htslib.org/doc/bgzip.html), [gunzip](https://linux.die.net/man/1/gunzip), [bwa](http://bio-bwa.sourceforge.net/), [samtools](http://www.htslib.org/), [gatk](https://anaconda.org/bioconda/gatk4)
- **OS:** Ubuntu 16.04 (Wintermute - research server)

## Table of contents

- [Database downloads for human_genomics_pipeline and vcf_annotation_pipeline](#database-downloads-for-humangenomicspipeline-and-vcfannotationpipeline)
  - [Table of contents](#table-of-contents)
    - [Reference human genome](#reference-human-genome)
      - [Option one: download from the GATK resource bundle](#option-one-download-from-the-gatk-resource-bundle)
      - [Option two: download from NCBI](#option-two-download-from-ncbi)
      - [Option three: download from UCSC](#option-three-download-from-ucsc)
    - [dbSNP database](#dbsnp-database)
      - [Option one: download from NCBI](#option-one-download-from-ncbi)
      - [Option two: download from the GATK resource bundle](#option-two-download-from-the-gatk-resource-bundle)
    - [VEP database](#vep-database)
    - [dbNSFP](#dbnsfp)
    - [Other vcf annotation databases](#other-vcf-annotation-databases)
      - [Option one: download from NCBI](#option-one-download-from-ncbi-1)
      - [Option two: download from the GATK resource bundle](#option-two-download-from-the-gatk-resource-bundle-1)
    - [CADD database](#cadd-database)

There are several large data file/databases that need to be downloaded to run the[human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) and [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline). There are also often many places from which you can download them from and many verisons or builds of the databases. Here I will describe how I downloaded and prepared the data/databases required to run the human_genomics_pipeline and vcf_annotation_pipeline against both GRCh37 and GRCh38. It won't cover all ways of obtaining this data or all versions/builds of the data, but it will cover a few sources/builds. The data I recommend downloading and using in these pipelines will be indicated by a :star:.

If possible, try and use databases downloaded from the same source, since issues typically arise when using reference human genomes and dbSNP databases from different sources. This is because different sources may label data differently (eg. [chromosome labeling and length](https://gatkforums.broadinstitute.org/gatk/discussion/11359/input-files-reference-and-features-have-incompatible-contigs)).

----

*For both the GRCh37 and GRCh38 runs of the human_genomics_pipeline and vcf_annotation_pipeline, I used the reference human genomes from the GATK resource bundle and the dbSNP databases from NCBI. I chose to use these reference genomes since it is in the same format as several vcf annotation databases that need to be used in the vcf_annotation_pipeline after a human_genomics_pipeline run. I chose to use the dbSNP databases from NCBI since they provide a much more recent build of the dbSNP database compared to that available in the GATK resource bundle. Additionally, they provide them in a GATK formatted version. However, I used the slightly older version of the dbSNP database (build 151) since the newest version (build153) format is incompatible with the reference human genome and downloaded from the GATK resource bundle that are required for the human_genomics_pipeline and vcf_annotation_pipeline (see [Compare dbSNP files/builds](https://github.com/leahkemp/documentation/blob/master/compare_dbSNP_databases.md) and [A cheeky peek at big genomic data files for human genomic pipelines](https://github.com/leahkemp/documentation/blob/master/cheeky_peek_big_genomic_data_files.md)).*

----

### Reference human genome

Along with the reference human genome sequence, you will also require it's associated fasta sequence dictionary file (.dict) and fasta index file (.fai) files. Some sources will have these files available, but if it isn't, these files will [need to be created]((https://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference)) using [GATK4](https://anaconda.org/bioconda/gatk4) and [samtools](http://www.htslib.org/).

Secondly, the index files (.amb, .ann, .bwt, .pac, .sa) for the reference human genome sequence need to be created using [bwa](http://bio-bwa.sourceforge.net/) with the ['bwtsw' argument](http://seqanswers.com/forums/showthread.php?t=3547) to ensure that the correct algorithm is used to handle a large whole genome sequence.

#### Option one: download from the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)

GATK resource bundle ftp site: ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle

:star:

```bash
# GRCh37
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.fasta.gz
gunzip ucsc.hg19.fasta.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.dict.gz
gunzip ucsc.hg19.dict.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.fasta.fai.gz
gunzip ucsc.hg19.fasta.fai.gz
bwa index -a bwtsw ucsc.hg19.fasta
```

```bash
# GRCh38
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.gz
gunzip Homo_sapiens_assembly38.fasta.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.dict
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.fai
bwa index -a bwtsw Homo_sapiens_assembly38.fasta
```

:star:

#### Option two: download from [NCBI](https://www.ncbi.nlm.nih.gov/genome/guide/human/?)

NCBI ftp site: ftp://ftp.ncbi.nlm.nih.gov:21/

```bash
# GRCh37
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.fna.gz
gunzip GRCh37_latest_genomic.fna.gz
bwa index -a bwtsw GRCh37_latest_genomic.fna
gatk CreateSequenceDictionary -R GRCh37_latest_genomic.fna
samtools faidx GRCh37_latest_genomic.fna
```

```bash
# GRCh38
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
gunzip GRCh38_latest_genomic.fna.gz
bwa index -a bwtsw GRCh38_latest_genomic.fna
gatk CreateSequenceDictionary -R GRCh38_latest_genomic.fna
samtools faidx GRCh38_latest_genomic.fna
```

#### Option three: download from [UCSC](https://hgdownload.soe.ucsc.edu/downloads.html)

UCSC ftp site: ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/

```bash
# GRCh37
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
bwa index -a bwtsw hg19.fa
gatk CreateSequenceDictionary -R hg19.fa
samtools faidx hg19.fa
```

```bash
# GRCh38
wget ftp://hgdownload.soe.ucsc.edu:21/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
bwa index -a bwtsw hg38.fa
gatk CreateSequenceDictionary -R hg38.fa
samtools faidx hg38.fa
```

### dbSNP database

Along with the dbSNP database, you will also need to download it's associated tabix index file (.tbi). Sometimes it is not available for download, so one can be created with [tabix](http://www.htslib.org/doc/tabix.html). In order to do this, our dbsnp vcf file [needs to be bgzf-compressed file](https://github.com/samtools/bcftools/issues/668) (also see [here](https://www.biostars.org/p/138514/)). Check the format of your dbsnp file with:

```bash
htsfile your_file.vcf.gz
```

Output should say something like:

```bash
your_file.vcf.gz:    VCF version 4.0 BGZf-compressed variant calling data
```

If the dbSNP databases is compressed in the wrong format, it will need to be unzipped with [gunzip](https://linux.die.net/man/1/gunzip), and rezipped into bgzip format with [bgzip](http://www.htslib.org/doc/bgzip.html) before [tabix](http://www.htslib.org/doc/tabix.html) can be used to create a tabix index file (.tbi)

#### Option one: download from [NCBI](https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/)

NCBI ftp site: ftp://ftp.ncbi.nlm.nih.gov:21/

```bash
# Build 153 - GRCh37
wget ftp://ftp.ncbi.nlm.nih.gov:21/snp/latest_release/VCF/GCF_000001405.25.gz
wget ftp://ftp.ncbi.nlm.nih.gov:21/snp/latest_release/VCF/GCF_000001405.25.gz.tbi
```

```bash
# Build 153 - GRCh38
wget ftp://ftp.ncbi.nlm.nih.gov:21/snp/latest_release/VCF/GCF_000001405.38.gz
wget ftp://ftp.ncbi.nlm.nih.gov:21/snp/latest_release/VCF/GCF_000001405.38.gz.tbi
```

```bash
# Build 151 - GRCh37
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz.tbi
```

```bash
# Build 151 - GRCh38
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz.tbi
```

:star:

```bash
# Build 151 - GRCh37 - GATK formatted
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/GATK/All_20180423.vcf.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/GATK/All_20180423.vcf.gz.tbi
```

```bash
# Build 151 - GRCh38 - GATK formatted
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/All_20180418.vcf.gz
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/All_20180418.vcf.gz.tbi
```

:star:

```bash
# Build 150 - GRCh37
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz
```

```bash
# Build 150 - GRCh38
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/All_20170710.vcf.gz
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/All_20170710.vcf.gz.tbi
```

```bash
# Build 150 - GRCh37 - GATK formatted
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/GATK/All_20170710.vcf.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/GATK/All_20170710.vcf.gz.tbi
```

```bash
# Build 150 - GRCh38 - GATK formatted
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/GATK/All_20170710.vcf.gz
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/GATK/All_20170710.vcf.gz.tbi
```

```bash
# Build 146 - GRCh38
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/dbsnp_146.hg38.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/dbsnp_146.hg38.vcf.gz.tbi
```

#### Option two: download from the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)

GATK resource bundle ftp site: ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle

```bash
# Build 138 - GRCh37
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/dbsnp_138.hg19.vcf.gz
gunzip dbsnp_138.hg19.vcf.gz
bgzip dbsnp_138.hg19.vcf
tabix dbsnp_138.hg19.vcf.gz
```

```bash
# Build 146 - GRCh38
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/dbsnp_146.hg38.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/dbsnp_146.hg38.vcf.gz.tbi
```

### VEP database

Create a conda environment

```bash
conda create --name download_data_env python=3.7
conda activate download_data_env
```

Install conda package of VEP

```bash
conda install -c bioconda ensembl-vep=99.2
```

Download VEP database with with the [ensembl-vep conda package](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/ensembl-vep/meta.yaml)

*These are very large files and will likely take some time to download*

```bash
# GRCh37/hg19
vep_install -a cf -s homo_sapiens -y GRCh37 -c /store/lkemp/publicData/vep/GRCh37 --CONVERT
# GRCh38/hg38
vep_install -a cf -s homo_sapiens -y GRCh38 -c /store/lkemp/publicData/vep/GRCh38 --CONVERT
```

### dbNSFP

We made a custom build database, but you can download a dbNSFP database from...

### Other vcf annotation databases

Along with these databases, you will also need to download it's associated tabix index file (.tbi). Sometimes it is not available for download, so one can be created with [tabix](http://www.htslib.org/doc/tabix.html). In order to do this, our dbsnp vcf file [needs to be bgzf-compressed file](https://github.com/samtools/bcftools/issues/668) (also see [here](https://www.biostars.org/p/138514/)). Check the format of your dbsnp file with:

```bash
htsfile your_file.vcf.gz
```

Output should say something like:

```bash
your_file.vcf.gz:    VCF version 4.0 BGZf-compressed variant calling data
```

If the dbSNP databases is compressed in the wrong format, it will need to be unzipped with [gunzip](https://linux.die.net/man/1/gunzip), and rezipped into bgzip format with [bgzip](http://www.htslib.org/doc/bgzip.html) before [tabix](http://www.htslib.org/doc/tabix.html) can be used to create a tabix index file (.tbi)

#### Option one: download from NCBI

NCBI ftp site: ftp://ftp.ncbi.nlm.nih.gov:21/

```bash
# Mills - GRCh37
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk:21/vol1/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz.tbi
```

```bash
# Mills - GRCh38

```

#### Option two: download from the GATK resource bundle

GATK resource bundle ftp site: ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle

:star:

```bash
# Mills - GRCh37
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
gunzip Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
bgzip Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
tabix Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz

# 1000G indel - GRCh37
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.gz
gunzip 1000G_phase1.indels.hg19.sites.vcf.gz
bgzip 1000G_phase1.indels.hg19.sites.vcf
tabix 1000G_phase1.indels.hg19.sites.vcf.gz

# 1000G snp - GRCh37
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
gunzip 1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
bgzip 1000G_phase1.snps.high_confidence.hg19.sites.vcf
tabix 1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz

# Omni - GRCh37
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/1000G_omni2.5.hg19.sites.vcf.gz
gunzip 1000G_omni2.5.hg19.sites.vcf.gz
bgzip 1000G_omni2.5.hg19.sites.vcf
tabix 1000G_omni2.5.hg19.sites.vcf.gz

# Hapmap3 - GRCh37
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/hapmap_3.3.hg19.sites.vcf.gz
gunzip hapmap_3.3.hg19.sites.vcf.gz
bgzip hapmap_3.3.hg19.sites.vcf
tabix hapmap_3.3.hg19.sites.vcf.gz
```

:star:

:star:

```bash
# Mills - GRCh38
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

# 1000G snp - GRCh38
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi

# Omni - GRCh38
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/1000G_omni2.5.hg38.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/1000G_omni2.5.hg38.vcf.gz.tbi

# Hapmap3 - GRCh38
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/hapmap_3.3.hg38.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/hapmap_3.3.hg38.vcf.gz.tbi

# Axiom exome plus - GRCh38
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi
```

:star:

### CADD database

This database is **massive** so it will take a long time to download (around 85G). aria2c can speed things up when downloading large files like these.

:star:

```bash
# GRCh37
aria2c -c -s 10 https://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/whole_genome_SNVs.tsv.gz
wget https://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/whole_genome_SNVs.tsv.gz.tbi
```

```bash
# GRCh38
aria2c -c -s 10 https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/whole_genome_SNVs.tsv.gz
wget https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/whole_genome_SNVs.tsv.gz.tbi
```

:star:
