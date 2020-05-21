# Benchmarking genomic pipelines

Created: 2020-04-22 13:37:04
Last modified: 2020/05/21 16:38:59

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
      - [Option two - bcftools](#option-two---bcftools)
      - [Option three and four - rtg vcfeval and hap.py](#option-three-and-four---rtg-vcfeval-and-happy)
    - [Download and prepare data](#download-and-prepare-data)
      - [Fastq data to run through pipelines](#fastq-data-to-run-through-pipelines)
      - [Truth vcf to compare output to](#truth-vcf-to-compare-output-to)
  - [Benchmarking](#benchmarking)
    - [human_genomics_pipeline and vcf_annotation_pipeline](#humangenomicspipeline-and-vcfannotationpipeline)
      - [Compare the truth and query vcf](#compare-the-truth-and-query-vcf)
        - [Option one - bcftools](#option-one---bcftools)
        - [Option two - vcftools](#option-two---vcftools)
        - [Option three - rtg vcfeval](#option-three---rtg-vcfeval)
        - [Option four - hap.py](#option-four---happy)
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

#### Option two - [bcftools](http://samtools.github.io/bcftools/bcftools.html)

```bash
conda install -c bioconda bcftools=1.10.2
conda install -c bioconda vcftools=1.1.16
```

#### Option three and four - rtg vcfeval and hap.py

[hap.py](https://github.com/Illumina/hap.py) and [RTG Tools](https://github.com/RealTimeGenomics/rtg-tools/tree/eb13bbb82d2fbeab7d54a92e8493ddd2acf0d349)

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

Install hap.py and RTG Tools (using the python helper script)

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

### human_genomics_pipeline and vcf_annotation_pipeline

See the results and settings of the pipeline runs on the Genome In A Bottle (GIAB) sample [NA12878](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) (NIST7035 and NIST7086) for:

- human_genomics_pipeline at https://github.com/ESR-NZ/human_genomics_pipeline/tree/bench1.0
- vcf_annotation_pipeline at https://github.com/ESR-NZ/vcf_annotation_pipeline/tree/bench1.0

#### Compare the truth and query vcf

Some chromosomes/chromosome patches are unique to the pipeline output vcf or the truth vcf. Variant calls on these regions could inflate the number of 'false-positive' and 'false-negative' calls. Therefore we will create a vcf file for both the pipeline output vcf and the truth vcf that only contain variants for chromosomes/chromosome patches that are common between them.

First, create separate vcf files for columns in known vcf (NIST7035 and NIST7086)

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/bench1.0/

bcftools view \
-s NIST7035 \
../known/project.NIST.hc.snps.indels.vcf \
-o ../known/project.NIST.hc.snps.indels.NIST7035.vcf

bcftools view \
-s NIST7086 \
../known/project.NIST.hc.snps.indels.vcf \
-o ../known/project.NIST.hc.snps.indels.NIST7086.vcf
```

Create regions file for all vcf files for comparison

```bash
# Pipeline output vcf (NIST7035)
grep -v "#" ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf | awk '{print $1}' | uniq > NIST7035_NIST_filtered.regions

# Known vcf (NIST7035)
grep -v "#" ../known/project.NIST.hc.snps.indels.NIST7035.vcf | awk '{print $1}' | uniq > ../known/project.NIST.hc.snps.indels.NIST7035.regions

# Pipeline output vcf (NIST7086)
grep -v "#" ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf | awk '{print $1}' | uniq > NIST7086_NIST_filtered.regions

# Known vcf (NIST7086)
grep -v "#" ../known/project.NIST.hc.snps.indels.NIST7086.vcf | awk '{print $1}' | uniq > ../known/project.NIST.hc.snps.indels.NIST7086.regions
```

Create a regions file that intersects both the pipeline output and truth vcf (NIST7035 and NIST7086)

```bash
# NIST7035
sort \
NIST7035_NIST_filtered.regions \
../known/project.NIST.hc.snps.indels.NIST7035.regions |
uniq -d > common_chroms_NIST7035_regions

# NIST7086
sort \
NIST7086_NIST_filtered.regions \
../known/project.NIST.hc.snps.indels.NIST7086.regions |
uniq -d > common_chroms_NIST7086_regions
```

The known sites file and query vcf file needs to be bgzipped and have a tabix index file (.tbi) (write to new files so as to not modify the original files):

```bash
# Pipeline output vcf (NIST7035)
bgzip < ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf > ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz
# Pipeline output vcf NIST7086)
bgzip < ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf > ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz
# Known vcf (NIST7035)
bgzip < ../known/project.NIST.hc.snps.indels.NIST7035.vcf > ../known/project.NIST.hc.snps.indels.NIST7035.vcf.gz
tabix ../known/project.NIST.hc.snps.indels.NIST7035.vcf.gz
# Known vcf (NIST7086)
bgzip < ../known/project.NIST.hc.snps.indels.NIST7086.vcf > ../known/project.NIST.hc.snps.indels.NIST7086.vcf.gz
tabix ../known/project.NIST.hc.snps.indels.NIST7086.vcf.gz
```

##### Option one - bcftools

1. Create intersection and complements of two sets saving the output in dir/*

- All chromosomes

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/bench1.0/
bcftools isec \
./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz \
../known/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
-p NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels \
--threads 8
```

```bash
bcftools isec \
./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz \
../known/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
-p NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels \
--threads 8
```

- Chromosome overlap adjusted

```bash
bcftools isec \
./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz \
../known/project.NIST.hc.snps.indels.NIST7035.vcf.gz \
--regions chr1,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr17_ctg5_hap1,chr17_gl000203_random,chr17_gl000204_random,chr17_gl000205_random,chr18,chr18_gl000207_random,chr19,chr19_gl000208_random,chr19_gl000209_random,chr1_gl000191_random,chr1_gl000192_random,chr2,chr20,chr21,chr22,chr3,chr4,chr4_ctg9_hap1,chr4_gl000193_random,chr4_gl000194_random,chr5,chr6,chr6_cox_hap2,chr6_dbb_hap3,chr6_mann_hap4,chr6_mcf_hap5,chr6_qbl_hap6,chr6_ssto_hap7,chr7,chr7_gl000195_random,chr8,chr9,chr9_gl000198_random,chr9_gl000199_random,chrM,chrUn_gl000211,chrUn_gl000212,chrUn_gl000213,chrUn_gl000214,chrUn_gl000215,chrUn_gl000216,chrUn_gl000217,chrUn_gl000218,chrUn_gl000219,chrUn_gl000220,chrUn_gl000221,chrUn_gl000222,chrUn_gl000223,chrUn_gl000224,chrUn_gl000225,chrUn_gl000226,chrUn_gl000228,chrUn_gl000229,chrUn_gl000230,chrUn_gl000231,chrUn_gl000232,chrUn_gl000233,chrUn_gl000234,chrUn_gl000235,chrUn_gl000237,chrUn_gl000238,chrUn_gl000239,chrUn_gl000240,chrUn_gl000241,chrUn_gl000242,chrUn_gl000243,chrUn_gl000247,chrX,chrY \
-p NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.chrom.adjusted \
--threads 8
```

```bash
bcftools isec \
./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz \
../known/project.NIST.hc.snps.indels.NIST7086.vcf.gz \
--regions chr1,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr17_ctg5_hap1,chr17_gl000203_random,chr17_gl000204_random,chr17_gl000205_random,chr18,chr18_gl000207_random,chr19,chr19_gl000208_random,chr19_gl000209_random,chr1_gl000191_random,chr1_gl000192_random,chr2,chr20,chr21,chr22,chr3,chr4,chr4_ctg9_hap1,chr4_gl000193_random,chr4_gl000194_random,chr5,chr6,chr6_cox_hap2,chr6_dbb_hap3,chr6_mann_hap4,chr6_mcf_hap5,chr6_qbl_hap6,chr6_ssto_hap7,chr7,chr7_gl000195_random,chr8,chr9,chr9_gl000198_random,chr9_gl000199_random,chrM,chrUn_gl000211,chrUn_gl000212,chrUn_gl000213,chrUn_gl000214,chrUn_gl000215,chrUn_gl000216,chrUn_gl000217,chrUn_gl000218,chrUn_gl000219,chrUn_gl000220,chrUn_gl000221,chrUn_gl000222,chrUn_gl000224,chrUn_gl000225,chrUn_gl000226,chrUn_gl000228,chrUn_gl000229,chrUn_gl000230,chrUn_gl000231,chrUn_gl000232,chrUn_gl000233,chrUn_gl000234,chrUn_gl000235,chrUn_gl000237,chrUn_gl000238,chrUn_gl000239,chrUn_gl000240,chrUn_gl000241,chrUn_gl000243,chrUn_gl000247,chrX,chrY \
-p NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.chrom.adjusted \
--threads 8
```

ERROR: We should be able to pass a regions file to the `--regions-file` flag in bcftols isec, however the current regions file is faulty, it either needs to be a tab-delimited file with CHROM, POS, columns. Otherwise it can be a VCF or BED file (see [here](http://samtools.github.io/bcftools/bcftools.html#common_options)). For now I'll state the chromosomes manually.

To get counts of variants/lines in each vcf output file, eg:

```bash
zgrep -v "#" ./NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels/0000.vcf | wc -l
zgrep -v "#" ./NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.chrom.adjusted/0000.vcf | wc -l
```

To get genotype quality scores (GQ)

```bash
vcftools \
--vcf ./NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.chrom.adjusted/0001.vcf \
--extract-FORMAT-info GQ \
--out test.vcf
```

Produce summary of vcf file

```bash
bcftools stats ./NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.chrom.adjusted/0001.vcf > ./NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.chrom.adjusted/0001.vcf.summary
```

##### Option two - vcftools

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/bench1.0/
# Extract all the GQ values
vcftools \
--gzvcf ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz \
--gzdiff ../known/project.NIST.hc.snps.indels.vcf.gz \
--diff-site \
--out ./NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels

# Summarise the GQ values

```

##### Option three - rtg vcfeval

We need to create an sdf file for the reference human genome (that was used in the benchmarking pipeline run). Create this with the rgt-tools format function:

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/bench1.0/
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

##### Option four - hap.py

Create a bed regions file from the known vcf (using bedops)

```bash
conda install -c bioconda bedops=2.4.39
```

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/bench1.0/
# NIST7035
vcf2bed < /store/mbenton/benchmarking/NIST7035_NIST.vcf > /store/lkemp/exome_project/benchmarking/NA12878_exome/known/NIST7035.bed
```

```bash
# NIST7035
mkdir happy_NIST7035
cd happy_NIST7035/

# CPU pipelines (h_g_p and v_a_p)
/store/lkemp/exome_project/benchmarking/NA12878_exome/hap.py-install/bin/hap.py \
/store/mbenton/benchmarking/project.NIST.hc.snps.indels.NIST7035.vcf \
/store/mbenton/benchmarking/NIST7035_NIST_filtered.vcf \
-f /store/lkemp/exome_project/benchmarking/NA12878_exome/known/NIST7035.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7035_cpu \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/benchmarking/NA12878_exome/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf

# GPU pipeline (parabricks)
/store/lkemp/exome_project/benchmarking/NA12878_exome/hap.py-install/bin/hap.py \
/store/mbenton/benchmarking/project.NIST.hc.snps.indels.NIST7035.vcf \
/store/mbenton/benchmarking/NIST7035_NIST.vcf \
-f /store/lkemp/exome_project/benchmarking/NA12878_exome/known/NIST7035.bed \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-o happy_NIST7035_gpu \
--engine=vcfeval \
--engine-vcfeval-template /store/lkemp/exome_project/benchmarking/NA12878_exome/hap.py-install/libexec/rtg-tools-install/ucsc.hg19.fasta.sdf
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
