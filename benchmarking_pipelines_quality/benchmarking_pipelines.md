# Benchmarking genomic pipelines

Created: 2020-04-22 13:37:04
Last modified: 2020/05/25 18:44:13

- **Aim:** Undertake benchmarking of genomics pipelines to test their quality for clinical use. 
- **Prerequisite software:** [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html), [bgzip](http://www.htslib.org/doc/bgzip.html), [tabix](http://www.htslib.org/doc/tabix.html)
- **OS:** Ubuntu 16.04 (Wintermute - research server)

The idea is to run these pipelines against the Genome In A Bottle (GIAB) sample [NA12878](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) and compare the quality of variant calls in their output vcf files. See the complementary docs for benchmarking the Nvidia Parabricks pipeline [here](https://github.com/ESR-NZ/ESR-Parabricks). Benchmarking will follow the best practices described in [this paper](https://www.nature.com/articles/s41587-019-0054-x).

## Table of contents

- [Benchmarking genomic pipelines](#benchmarking-genomic-pipelines)
  - [Table of contents](#table-of-contents)
  - [Setup](#setup)
    - [Install benchmarking software](#install-benchmarking-software)
    - [Download and prepare data](#download-and-prepare-data)
      - [Fastq data to run through pipelines](#fastq-data-to-run-through-pipelines)
      - [Truth vcf to compare output to](#truth-vcf-to-compare-output-to)
    - [Other setup](#other-setup)
      - [sdf files for hap.py + RTG vcfeval](#sdf-files-for-happy--rtg-vcfeval)
      - [bed regions files for hap.py + RTG vcfeval](#bed-regions-files-for-happy--rtg-vcfeval)
      - [Formatting vcf files for all vcf comparisons](#formatting-vcf-files-for-all-vcf-comparisons)
  - [Benchmarking](#benchmarking)
    - [intra_truth_comparison](#intratruthcomparison)
      - [Compared with bedtools intersect](#compared-with-bedtools-intersect)
        - [Compared with hap.py](#compared-with-happy)
        - [NIST7035 ('baseline') compared to NIST7086 ('truth')](#nist7035-baseline-compared-to-nist7086-truth)
        - [NIST7086 ('baseline') compared to NIST7035 ('truth')](#nist7086-baseline-compared-to-nist7035-truth)
    - [bench 1.0](#bench-10)
      - [human_genomics_pipeline + minimal vcf_annotation_pipeline](#humangenomicspipeline--minimal-vcfannotationpipeline)
        - [Compared with bedtools intersect](#compared-with-bedtools-intersect-1)
        - [Compared with hap.py + RTG tools](#compared-with-happy--rtg-tools)
      - [parabricks germline pipeline](#parabricks-germline-pipeline)
        - [Compared with bedtools intersect](#compared-with-bedtools-intersect-2)
        - [Compared with hap.py + RTG tools](#compared-with-happy--rtg-tools-1)
    - [bench 1.1](#bench-11)
      - [human_genomics_pipeline + minimal vcf_annotation_pipeline](#humangenomicspipeline--minimal-vcfannotationpipeline-1)
        - [Compared with bedtools intersect](#compared-with-bedtools-intersect-3)
        - [Compared with hap.py + RTG tools](#compared-with-happy--rtg-tools-2)
      - [parabricks germline pipeline](#parabricks-germline-pipeline-1)
        - [Compared with bedtools intersect](#compared-with-bedtools-intersect-4)
        - [Compared with hap.py + RTG tools](#compared-with-happy--rtg-tools-3)
  - [Results of benchmarking](#results-of-benchmarking)

## Setup

### Install benchmarking software

Create a conda environment for installing hap.py (we need python 2 in order to get htslib2 that is required for installation of hap.py)

```bash
conda create -n happy_install_env python=2.7.15
conda activate happy_install_env
```

Install dependencies for [hap.py](https://github.com/Illumina/hap.py) and [RTG Tools](https://github.com/RealTimeGenomics/rtg-tools/tree/eb13bbb82d2fbeab7d54a92e8493ddd2acf0d349)

```bash
conda install -c conda-forge cython=0.29.15
conda install -c conda-forge scipy=1.2.1 # Also installs numpy=1.16.5 dependency
conda install -c conda-forge pandas=0.24.2
conda install -c bioconda pybedtools=0.8.1 # Also installs pysam=0.15.3 dependency
conda install -c bioconda bx-python=0.8.8
```

Install hap.py and RTG Tools (using the python helper script)

```bash
cd /store/lkemp/exome_project/benchmarking/
git clone git@github.com:Illumina/hap.py.git
cd hap.py
python install.py /store/lkemp/exome_project/benchmarking/hap.py-install \
--with-rtgtools \
--no-tests
```

Install other benchmarking software in a new conda environment

```bash
conda deactivate
conda create -n benchmarking_env python=3.7
conda activate benchmarking_env
conda install -c bioconda bedtools=2.29.2
conda install -c bioconda bcftools=1.10.2
conda install -c bioconda bedops=2.4.39
```

### Download and prepare data

#### Fastq data to run through pipelines

Download

```bash
cd /store/lkemp/publicData/exomes/NA12878_exome/
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/Garvan_NA12878_HG001_HiSeq_Exome.README
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L002_R1_001.fastq.gz
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L002_R2_001.fastq.gz
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7086_CGTACTAG_L001_R1_001.fastq.gz
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7086_CGTACTAG_L001_R2_001.fastq.gz
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7086_CGTACTAG_L002_R1_001.fastq.gz
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7086_CGTACTAG_L002_R2_001.fastq.gz
```

Note: there can be some issues with correctly downloading the fastq files (likely due to the ESR proxy). The below can be run to check the integrity of the gunzip files, no error suggests the archive is OK.

```bash
for i in NIST*; do
  echo "$i" ;
  gunzip -t "$i"
  echo "...done..."
done
```

Collapse pooled runs

```bash
# NIST7035
cat NIST7035*_R1_001.fastq.gz > NIST7035_NIST_R1.fastq.gz
cat NIST7035*_R2_001.fastq.gz > NIST7035_NIST_R2.fastq.gz

# NIST7086
cat NIST7086*_R1_001.fastq.gz > NIST7086_NIST_R1.fastq.gz
cat NIST7086*_R2_001.fastq.gz > NIST7086_NIST_R2.fastq.gz
```

*Note. NIST7035 and NIST7086 represent 2 vials of the same sample (NA12878)*

#### Truth vcf to compare output to

```bash
# VCF
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST.hc.snps.indels.vcf
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST.hc.snps.indels.vcf.idx
# regions bed
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz
```

Create separate vcf files for columns in known vcf (NIST7035 and NIST7086)

```bash
bcftools view \
-s NIST7035 \
project.NIST.hc.snps.indels.vcf \
-o project.NIST.hc.snps.indels.NIST7035.vcf \
--exclude-private

bcftools view \
-s NIST7086 \
project.NIST.hc.snps.indels.vcf \
-o project.NIST.hc.snps.indels.NIST7086.vcf \
--exclude-private
```

See [here](https://github.com/ga4gh/benchmarking-tools/blob/master/resources/high-confidence-sets/giab.md) for more info on GIAB resources

### Other setup

#### sdf files for hap.py + RTG vcfeval

Create an sdf file for the reference human genome (that was used in the benchmarking pipeline run). Create this with the RTG-tools format function:

- ucsc.hg19.fasta

```bash
/store/lkemp/exome_project/benchmarking/hap.py-install/libexec/rtg-tools-install/rtg format \
--output /store/lkemp/exome_project/benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
/store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta
```

#### bed regions files for hap.py + RTG vcfeval

Create a bed regions file from the known vcfs (using bedops)

```bash
vcf2bed < project.NIST.hc.snps.indels.NIST7035.vcf > project.NIST.hc.snps.indels.NIST7035.bed
vcf2bed < project.NIST.hc.snps.indels.NIST7086.vcf > project.NIST.hc.snps.indels.NIST7086.bed
```

#### Formatting vcf files for all vcf comparisons

The vcf files for comparison need to be bgzipped and have a tabix index file (.tbi) (write to new files so as to not modify the original files):

- Known vcf

```bash
bgzip < project.NIST.hc.snps.indels.NIST7035.vcf > project.NIST.hc.snps.indels.NIST7035.vcf.gz
tabix project.NIST.hc.snps.indels.NIST7035.vcf.gz
bgzip < project.NIST.hc.snps.indels.NIST7086.vcf > project.NIST.hc.snps.indels.NIST7086.vcf.gz
tabix project.NIST.hc.snps.indels.NIST7086.vcf.gz
```

- bench1.0

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/bench1.0/

bgzip < ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf > ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz
bgzip < ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf > ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz
```

- bench1.1

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/bench1.1/
# Pipeline output vcf (NIST7035)
bgzip < ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf > ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz
# Pipeline output vcf (NIST7086)
bgzip < ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf > ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz
```

## Benchmarking

human_genomic_pipeline will undertake pre-processing and variant calling. Because variant filtering occurs with vcf annotation_pipeline, we will benchmark the vcf files that have gone through both pipelines. However, we will use a minimal version of the vcf_annotation_pipeline since the annotation rules are not required for benchmarking.

See the results and settings of the pipeline runs on the Genome In A Bottle (GIAB) sample [NA12878](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) (NIST7035 and NIST7086) for:

- human_genomics_pipeline at https://github.com/ESR-NZ/human_genomics_pipeline/tree/bench1.0
- vcf_annotation_pipeline at https://github.com/ESR-NZ/vcf_annotation_pipeline/tree/bench1.0
- parabricks germline pipeline at ...

### intra_truth_comparison

Compare the two truth sample vcfs (NIST7035 and NIST7086)

#### Compared with bedtools intersect

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/intra_truth_comparison/
mkdir intersect_project.NIST.hc.snps.indels.NIST7035_v_project.NIST.hc.snps.indels.NIST7086
cd intersect_project.NIST.hc.snps.indels.NIST7035_v_project.NIST.hc.snps.indels.NIST7086
```

```bash
# Common
bedtools intersect \
-a /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
-b /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
> common_project.NIST.hc.snps.indels.NIST7035_v_project.NIST.hc.snps.indels.NIST7086.vcf

# Unique NIST7035
bedtools intersect \
-a /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
-b /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
-v \
> unique_project.NIST.hc.snps.indels.NIST7035.vcf

# Unique NIST7086
bedtools intersect \
-a /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
-b /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
-v \
> unique_project.NIST.hc.snps.indels.NIST7086.vcf
```

```bash
grep -v "#" common_project.NIST.hc.snps.indels.NIST7035_v_project.NIST.hc.snps.indels.NIST7086.vcf | wc -l
grep -v "#" unique_project.NIST.hc.snps.indels.NIST7035.vcf | wc -l
grep -v "#" unique_project.NIST.hc.snps.indels.NIST7086.vcf | wc -l
```

##### Compared with hap.py

##### NIST7035 ('baseline') compared to NIST7086 ('truth')

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/intra_truth_comparison/
mkdir happy_project.NIST.hc.snps.indels.NIST7035_v_project.NIST.hc.snps.indels.NIST7086
cd happy_project.NIST.hc.snps.indels.NIST7035_v_project.NIST.hc.snps.indels.NIST7086
```

```bash
/store/lkemp/exome_project/benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_project.NIST.hc.snps.indels.NIST7035_v_project.NIST.hc.snps.indels.NIST7086 \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

##### NIST7086 ('baseline') compared to NIST7035 ('truth')

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/intra_truth_comparison/
mkdir happy_project.NIST.hc.snps.indels.NIST7086_v_project.NIST.hc.snps.indels.NIST7035
cd happy_project.NIST.hc.snps.indels.NIST7086_v_project.NIST.hc.snps.indels.NIST7035
```

```bash
/store/lkemp/exome_project/benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_project.NIST.hc.snps.indels.NIST7086_v_project.NIST.hc.snps.indels.NIST7035 \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

### bench 1.0

#### human_genomics_pipeline + minimal vcf_annotation_pipeline

##### Compared with bedtools intersect

- NIST7035

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/bench1.0/
mkdir intersect_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035
cd intersect_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035
```

```bash
# Common
bedtools intersect \
-a ../vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz \
-b /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
> common_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035.vcf

# Unique truth
bedtools intersect \
-a ../vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz \
-b /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
-v \
> unique_NIST7035_NIST_filtered.vcf

# Unique query
bedtools intersect \
-a /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
-b ../vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz \
-v \
> unique_project.NIST.hc.snps.indels.NIST7035.vcf
```

```bash
grep -v "#" common_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035.vcf | wc -l
grep -v "#" unique_NIST7035_NIST_filtered.vcf | wc -l
grep -v "#" unique_project.NIST.hc.snps.indels.NIST7035.vcf | wc -l
```

- NIST7086

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/bench1.0/
mkdir intersect_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086
cd intersect_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086
```

```bash
# Common
bedtools intersect \
-a ../vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz \
-b /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
> common_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086.vcf

# Unique truth
bedtools intersect \
-a ../vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz \
-b /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
-v \
> unique_NIST7086_NIST_filtered.vcf

# Unique query
bedtools intersect \
-a /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
-b ../vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz \
-v \
> unique_project.NIST.hc.snps.indels.NIST7086.vcf
```

```bash
grep -v "#" common_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086.vcf | wc -l
grep -v "#" unique_NIST7086_NIST_filtered.vcf | wc -l
grep -v "#" unique_project.NIST.hc.snps.indels.NIST7086.vcf | wc -l
```

##### Compared with hap.py + RTG tools

- NIST7035

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/bench1.0/
mkdir happy_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035
cd happy_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035
```

```bash
/store/lkemp/exome_project/benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

- NIST7086

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/bench1.0/
mkdir happy_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086
cd happy_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086
```

```bash
/store/lkemp/exome_project/benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

#### parabricks germline pipeline

##### Compared with bedtools intersect

- NIST7035

```bash
cd /home/lkemp/benchmarking_parabricks/bench1.0/ # ORAC
mkdir intersect_NIST7035_NIST_v_project.NIST.hc.snps.indels.NIST7035
cd intersect_NIST7035_NIST_v_project.NIST.hc.snps.indels.NIST7035
```

```bash
# Common
bedtools intersect \
-a ../NIST7035_NIST.vcf \
-b /NGS/scratch/KSCBIOM/HumanGenomics/publicData/human_refs/GIAB/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
> common_NIST7035_NIST_v_project.NIST.hc.snps.indels.NIST7035.vcf

# Unique truth
bedtools intersect \
-a ../NIST7035_NIST.vcf \
-b /NGS/scratch/KSCBIOM/HumanGenomics/publicData/human_refs/GIAB/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
-v \
> unique_NIST7035_NIST.vcf

# Unique query
bedtools intersect \
-a /NGS/scratch/KSCBIOM/HumanGenomics/publicData/human_refs/GIAB/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
-b ../NIST7035_NIST.vcf \
-v \
> unique_project.NIST.hc.snps.indels.NIST7035.vcf
```

```bash
grep -v "#" common_NIST7035_NIST_v_project.NIST.hc.snps.indels.NIST7035.vcf | wc -l
grep -v "#" unique_NIST7035_NIST.vcf | wc -l
grep -v "#" unique_project.NIST.hc.snps.indels.NIST7035.vcf | wc -l
```

- NIST7086

```bash
cd /home/lkemp/benchmarking_parabricks/bench1.0/ # ORAC
mkdir intersect_NIST7086_NIST_v_project.NIST.hc.snps.indels.NIST7086
cd intersect_NIST7086_NIST_v_project.NIST.hc.snps.indels.NIST7086
```

```bash
# Common
bedtools intersect \
-a ../NIST7086_NIST.vcf \
-b /NGS/scratch/KSCBIOM/HumanGenomics/publicData/human_refs/GIAB/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
> common_NIST7086_NIST_v_project.NIST.hc.snps.indels.NIST7086.vcf

# Unique truth
bedtools intersect \
-a ../NIST7086_NIST.vcf \
-b /NGS/scratch/KSCBIOM/HumanGenomics/publicData/human_refs/GIAB/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
-v \
> unique_NIST7086_NIST.vcf

# Unique query
bedtools intersect \
-a /NGS/scratch/KSCBIOM/HumanGenomics/publicData/human_refs/GIAB/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
-b ../NIST7086_NIST.vcf \
-v \
> unique_project.NIST.hc.snps.indels.NIST7086.vcf
```

```bash
grep -v "#" common_NIST7086_NIST_v_project.NIST.hc.snps.indels.NIST7086.vcf | wc -l
grep -v "#" unique_NIST7086_NIST.vcf | wc -l
grep -v "#" unique_project.NIST.hc.snps.indels.NIST7086.vcf | wc -l
```

##### Compared with hap.py + RTG tools


### bench 1.1

#### human_genomics_pipeline + minimal vcf_annotation_pipeline

##### Compared with bedtools intersect

- NIST7035

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/bench1.1/
mkdir intersect_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035
cd intersect_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035
```

```bash
# Common
bedtools intersect \
-a ../vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz \
-b /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
> common_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035.vcf

# Unique truth
bedtools intersect \
-a ../vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz \
-b /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
-v \
> unique_NIST7035_NIST_filtered.vcf

# Unique query
bedtools intersect \
-a /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
-b ../vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz \
-v \
> unique_project.NIST.hc.snps.indels.NIST7035.vcf
```

```bash
grep -v "#" common_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035.vcf | wc -l
grep -v "#" unique_NIST7035_NIST_filtered.vcf | wc -l
grep -v "#" unique_project.NIST.hc.snps.indels.NIST7035.vcf | wc -l
```

- NIST7086

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/bench1.1/
mkdir intersect_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086
cd intersect_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086
```

```bash
# Common
bedtools intersect \
-a ../vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz \
-b /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
> common_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086.vcf

# Unique truth
bedtools intersect \
-a ../vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz \
-b /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
-v \
> unique_NIST7086_NIST_filtered.vcf

# Unique query
bedtools intersect \
-a /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
-b ../vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz \
-v \
> unique_project.NIST.hc.snps.indels.NIST7086.vcf
```

```bash
grep -v "#" common_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086.vcf | wc -l
grep -v "#" unique_NIST7086_NIST_filtered.vcf | wc -l
grep -v "#" unique_project.NIST.hc.snps.indels.NIST7086.vcf | wc -l
```

##### Compared with hap.py + RTG tools

- NIST7035

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/bench1.1/
mkdir happy_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035
cd happy_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035
```

```bash
/store/lkemp/exome_project/benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
../vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

- NIST7086

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/bench1.1/
mkdir happy_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086
cd happy_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086
```

```bash
/store/lkemp/exome_project/benchmarking/hap.py-install/bin/hap.py \
/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
../vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz \
-f /store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/benchmarking/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
```

#### parabricks germline pipeline

##### Compared with bedtools intersect

##### Compared with hap.py + RTG tools

## Results of benchmarking

See [here](https://leahkemp.github.io/documentation/benchmarking_pipelines/benchmarking_pipeline_results.html)
