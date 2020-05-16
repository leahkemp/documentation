# Benchmarking genomic pipelines

Created: 2020-04-22 13:37:04
Last modified: 2020/05/16 20:43:23

- **Aim:** Undertake benchmarking of genomics pipelines to test their quality for clinical use. 
- **Prerequisite software:** [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html), [bgzip](http://www.htslib.org/doc/bgzip.html), [tabix](http://www.htslib.org/doc/tabix.html)
- **OS:** Ubuntu 16.04 (Wintermute - research server)

The idea is to run these pipelines against the Genome In A Bottle (GIAB) sample [NA12878](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) and compare the quality of variant calls in their output vcf files. See the complementary docs for benchmarking the Nvidia Parabricks pipeline [here](https://github.com/ESR-NZ/ESR-Parabricks). Benchmarking will follow the best practices described in [this paper](https://www.nature.com/articles/s41587-019-0054-x).

## Table of contents

- [Benchmarking genomic pipelines](#benchmarking-genomic-pipelines)
  - [Table of contents](#table-of-contents)
  - [Setup](#setup)
    - [Install benchmarking software](#install-benchmarking-software)
      - [Option one - vcftools](#option-one---vcftools)
      - [Option two - hap.py and [RTG Tools](https://github.com/RealTimeGenomics/rtg-tools/tree/eb13bbb82d2fbeab7d54a92e8493ddd2acf0d349)](#option-two---happy-and-rtg-tools)
    - [Download and prepare data](#download-and-prepare-data)
      - [Fastq data to run through pipelines](#fastq-data-to-run-through-pipelines)
      - [Truth vcf to compare output to](#truth-vcf-to-compare-output-to)
  - [Benchmarking](#benchmarking)
    - [human_genomics_pipeline and [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline)](#humangenomicspipeline-and-vcfannotationpipeline)
      - [Compare the truth and query vcf](#compare-the-truth-and-query-vcf)
        - [Option one - bcftools](#option-one---bcftools)
        - [Option two - vcftools](#option-two---vcftools)
        - [Option three - rtg vcfeval](#option-three---rtg-vcfeval)
    - [Nvidia parabricks germline](#nvidia-parabricks-germline)
  - [Results of benchmarking](#results-of-benchmarking)
  - [Notes](#notes)

## Setup

### Install benchmarking software

#### Option one - [vcftools](http://vcftools.sourceforge.net/)

Create a conda environment with vcftools installed

```bash
conda create -n vcftools_env python=3.7
conda activate vcftools_env
conda install -c bioconda=0.1.16
```

#### Option two - [hap.py](https://github.com/Illumina/hap.py) and [RTG Tools](https://github.com/RealTimeGenomics/rtg-tools/tree/eb13bbb82d2fbeab7d54a92e8493ddd2acf0d349)

Create a conda environment (we need python 2 in order to get htslib2 that is required for installation)

```bash
conda create -n benchmarking_env python=2.7.5
conda activate benchmarking_env
```

Install dependencies

```bash
conda install -c conda-forge cython=0.29.15
conda install -c conda-forge scipy=1.2.1 # Also installs numpy=1.16.5 dependency
conda install -c conda-forge pandas=0.24.2
conda install -c bioconda pybedtools=0.8.1 # Also installs pysam=0.15.3 dependency
conda install -c bioconda bx-python=0.8.8
```

Install hap.py and RTG Tools (using the python and the helper script)

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/
git clone git@github.com:Illumina/hap.py.git
cd hap.py
python install.py /store/lkemp/exome_project/benchmarking/NA12878_exome/hap.py-install \
--with-rtgtools \
--no-tests
```

### Download and prepare data

#### Fastq data to run through pipelines

Download

```bash
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L002_R1_001.fastq.gz
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L002_R2_001.fastq.gz
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7086_CGTACTAG_L001_R1_001.fastq.gz
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7086_CGTACTAG_L001_R2_001.fastq.gz
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7086_CGTACTAG_L002_R1_001.fastq.gz
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7086_CGTACTAG_L002_R2_001.fastq.gz
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

See [here](https://github.com/ga4gh/benchmarking-tools/blob/master/resources/high-confidence-sets/giab.md) for more info on GIAB resources

## Benchmarking

human_genomic_pipeline will undertake pre-processing and variant calling. Because variant filtering occurs with vcf annotation_pipeline, we will benchmark the vcf files that have gone through both pipelines. However, we will use a minimal version of the vcf_annotation_pipeline since the annotation rules are not required for benchmarking.

### [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) and [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline)

See the results and settings of the pipeline runs on the Genome In A Bottle (GIAB) sample [NA12878](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) (NIST7035 and NIST7086) for:

- human_genomics_pipeline at https://github.com/ESR-NZ/human_genomics_pipeline/tree/bench1.0
- vcf_annotation_pipeline at https://github.com/ESR-NZ/vcf_annotation_pipeline/tree/bench1.0

#### Compare the truth and query vcf

##### Option one - bcftools

```bash
conda install -c bioconda bcftools=1.10.2
```

```bash
bcftools isec \
./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz \
../known/project.NIST.hc.snps.indels.vcf.gz \
-p NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels
```

```bash
bcftools isec \
./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz \
../known/project.NIST.hc.snps.indels.vcf.gz \
-p NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels
```

##### Option two - vcftools

```bash
vcftools \
--gzvcf ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz \
--gzdiff ../known/project.NIST.hc.snps.indels.vcf.gz \
--diff-site \
--out ./NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels
```

##### Option three - rtg vcfeval

The known sites file and query vcf file needs to be bgzipped and have a tabix index file (.tbi) (write to new files so as to not modify the original files):

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/bench1.0/
# Query vcf (NIST7035)
bgzip < ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf > ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz
# Query vcf (NIST7086)
bgzip < ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf > ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz
# Known vcf
bgzip < ../known/project.NIST.hc.snps.indels.vcf> ../known/project.NIST.hc.snps.indels.vcf.gz
tabix ../known/project.NIST.hc.snps.indels.vcf.gz
```

We also need to create an sdf file for the reference human genome (that was used in the benchmarking pipeline run). Create this with the rgt-tools format function:
FIX
```bash
../hap.py-install/libexec/rtg-tools-install/rtg format \
--output ../hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
/store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta
```

Run vcfeval

```bash
# NIST7035
../hap.py-install/libexec/rtg-tools-install/rtg vcfeval \
--baseline ../known/project.NIST.hc.snps.indels.vcf.gz \
--calls ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz \
--sample NIST7035,NIST7035_NIST \
--output vcfeval_NIST7035 \
--template ../hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--evaluation-regions ../known/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
--output-mode split \
--threads 12

# NIST7086
../hap.py-install/libexec/rtg-tools-install/rtg vcfeval \
--baseline ../known/project.NIST.hc.snps.indels.vcf.gz \
--calls ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz \
--sample NIST7086,NIST7086_NIST \
--output vcfeval_NIST7086 \
--template ../hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf \
--evaluation-regions ../known/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
--output-mode split \
--threads 12
```

This will generate VCF files containing called variants that were in the truth VCF (tp), called variants that were not in the truth VCF (fp) and truth variants that were not in the called variants (fn) (for a more in depth explanation see [here](https://github.com/ga4gh/benchmarking-tools/blob/master/doc/ref-impl/intermediate.md)). It will also output a count of variants in each of these categories.

View ROC plots using `rtg rocplot`

```bash

```

### [Nvidia parabricks germline](https://github.com/ESR-NZ/ESR-Parabricks)

See the results and settings of the pipeline runs on the Genome In A Bottle (GIAB) sample [NA12878](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) (NIST7035 and NIST7086) for:

- Nvidia parabricks germline at ...

Located at `/usr/local/bin/pbrun`

Run:

```bash

```

## Results of benchmarking

See [this doc](https://leahkemp.github.io/documentation/benchmarking_pipelines/benchmarking_pipeline_results.html)

## Notes

- I found it difficult to use the hap.py wrapper for vcfeval since there were parameters I wasn't able to pass to vcfeval (such as specifying a sample in a multi sample vcf file and specifying the output mode).
