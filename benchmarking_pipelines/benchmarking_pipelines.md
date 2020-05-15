# Benchmarking genomic pipelines

Created: 2020-04-22 13:37:04
Last modified: 2020/05/15 12:15:59

- **Aim:** Undertake benchmarking of genomics pipelines to test their quality for clinical use. 
- **Prerequisite software:** [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html), [bgzip](http://www.htslib.org/doc/bgzip.html), [tabix](http://www.htslib.org/doc/tabix.html)
- **OS:** Ubuntu 16.04 (Wintermute - research server)

The ides is to run these pipelines against the Genome In A Bottle (GIAB) sample [NA12878](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) and compare the quality of varaiant calls in their output vcf files. See the complementary docs for benchmarking the Nvidia Parabricks pipeline [here](https://github.com/ESR-NZ/ESR-Parabricks). Benchmakring will follow the best practices described in [this paper](https://www.nature.com/articles/s41587-019-0054-x).

## Table of contents

- [Benchmarking genomic pipelines](#benchmarking-genomic-pipelines)
  - [Table of contents](#table-of-contents)
  - [Setup](#setup)
    - [Install benchmarking software](#install-benchmarking-software)
      - [hap.py and [RTG Tools](https://github.com/RealTimeGenomics/rtg-tools/tree/eb13bbb82d2fbeab7d54a92e8493ddd2acf0d349)](#happy-and-rtg-tools)
    - [Download and prepare data](#download-and-prepare-data)
      - [Fastq data to run through pipelines](#fastq-data-to-run-through-pipelines)
      - [Truth vcf to compare output to](#truth-vcf-to-compare-output-to)
  - [Benchmarking](#benchmarking)
    - [human_genomics_pipeline and [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline)](#humangenomicspipeline-and-vcfannotationpipeline)
      - [Compare the truth and query vcf with vcfeval](#compare-the-truth-and-query-vcf-with-vcfeval)
    - [Nvidia parabricks germline](#nvidia-parabricks-germline)
  - [Notes](#notes)

## Setup

### Install benchmarking software

#### [hap.py](https://github.com/Illumina/hap.py) and [RTG Tools](https://github.com/RealTimeGenomics/rtg-tools/tree/eb13bbb82d2fbeab7d54a92e8493ddd2acf0d349)

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

#### Compare the truth and query vcf with vcfeval

The known sites file and query vcf file needs to be bgzipped and have a tabix index file (.tbi) (write to new files so as to not modify the original files):

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/
# Query vcf (NIST7035)
bgzip < ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf > ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz
# Query vcf (NIST7086)
bgzip < ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf > ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz
tabix ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz
# Known vcf
bgzip < ./known/project.NIST.hc.snps.indels.vcf> ./known/project.NIST.hc.snps.indels.vcf.gz
tabix ./known/project.NIST.hc.snps.indels.vcf.gz
```

We also need to create an sdf file for the reference human genome (that was used in the benchmarking pipeline run). Create this with the rgt-tools format function:

```bash
./hap.py-install/libexec/rtg-tools-install/rtg format \
--output ./hap.py-install/libexec/rtg-tools-install/Homo_sapiens_assembly38.fasta.sdf \
/store/lkemp/publicData/referenceGenome/gatkBundle/GRCh38/Homo_sapiens_assembly38.fasta
```

Run `rtg vcfeval`

```bash
# NIST7035
./hap.py-install/libexec/rtg-tools-install/rtg vcfeval \
--baseline ./known/project.NIST.hc.snps.indels.vcf.gz \
--calls ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz \
--sample NIST7035,NIST7035_NIST \
--output vcfeval_NIST7035 \
--template ./hap.py-install/libexec/rtg-tools-install/Homo_sapiens_assembly38.fasta.sdf \
--evaluation-regions ./known/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
--output-mode split \
--threads 12

# NIST7086
./hap.py-install/libexec/rtg-tools-install/rtg vcfeval \
--baseline ./known/project.NIST.hc.snps.indels.vcf.gz \
--calls ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz \
--sample NIST7086,NIST7086_NIST \
--output vcfeval_NIST7086 \
--template ./hap.py-install/libexec/rtg-tools-install/Homo_sapiens_assembly38.fasta.sdf \
--evaluation-regions ./known/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
--output-mode split \
--threads 12
```

This will generate VCF files containing called variants that were in the truth VCF (tp), called variants that were not in the truth VCF (fp) and truth variants that were not in the called variants (fn) (for a more in depth explanation see [here](https://github.com/ga4gh/benchmarking-tools/blob/master/doc/ref-impl/intermediate.md)). It will also output a count of variants in each of these categories.

View ROC plots using `rtg rocplot`

```bash

```

Run `vcftools`

- Chromosomes in known

```bash
zgrep -v "#" ./known/project.NIST.hc.snps.indels.vcf.gz | awk '{print $1}' | uniq
```

```bash
chr21
chr22
chrX
chrY
chr1_gl000191_random
chr1_gl000192_random
chr4_ctg9_hap1
chr4_gl000193_random
chr4_gl000194_random
chr6_apd_hap1
chr6_cox_hap2
chr6_dbb_hap3
chr6_mann_hap4
chr6_mcf_hap5
chr6_qbl_hap6
chr6_ssto_hap7
chr7_gl000195_random
chr9_gl000198_random
chr9_gl000199_random
chr17_ctg5_hap1
chr17_gl000203_random
chr17_gl000204_random
chr17_gl000205_random
chr18_gl000207_random
chr19_gl000208_random
chr19_gl000209_random
chrUn_gl000211
chrUn_gl000212
chrUn_gl000213
chrUn_gl000214
chrUn_gl000215
chrUn_gl000216
chrUn_gl000217
chrUn_gl000218
chrUn_gl000219
chrUn_gl000220
chrUn_gl000221
chrUn_gl000222
chrUn_gl000223
chrUn_gl000224
chrUn_gl000225
chrUn_gl000226
chrUn_gl000228
chrUn_gl000229
chrUn_gl000230
chrUn_gl000231
chrUn_gl000232
chrUn_gl000233
chrUn_gl000234
chrUn_gl000235
chrUn_gl000237
chrUn_gl000238
chrUn_gl000239
chrUn_gl000240
chrUn_gl000241
chrUn_gl000242
chrUn_gl000243
chrUn_gl000247
chrUn_gl000248
```

- Chromosomes in pipeline output

```bash
zgrep -v "#" ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf.gz | awk '{print $1}' | uniq
```

```bash
chr1
chr2
chr3
chr4
chr5
chr6
chr7
chr8
chr9
chr10
chr11
chr12
chr13
chr14
chr15
chr16
chr17
chr18
chr19
chr20
chr21
chr22
chrX
chrY
chrM
chr1_KI270706v1_random
chr1_KI270707v1_random
chr1_KI270708v1_random
chr1_KI270709v1_random
chr1_KI270710v1_random
chr1_KI270711v1_random
chr1_KI270712v1_random
chr1_KI270713v1_random
chr1_KI270714v1_random
chr2_KI270715v1_random
chr2_KI270716v1_random
chr3_GL000221v1_random
chr4_GL000008v2_random
chr5_GL000208v1_random
chr9_KI270718v1_random
chr9_KI270719v1_random
chr9_KI270720v1_random
chr11_KI270721v1_random
chr14_GL000009v2_random
chr14_GL000225v1_random
chr14_KI270722v1_random
chr14_GL000194v1_random
chr14_KI270723v1_random
chr14_KI270724v1_random
chr14_KI270725v1_random
chr14_KI270726v1_random
chr16_KI270728v1_random
chr17_GL000205v2_random
chr17_KI270729v1_random
chr17_KI270730v1_random
chr22_KI270731v1_random
chr22_KI270732v1_random
chr22_KI270733v1_random
chr22_KI270734v1_random
chr22_KI270735v1_random
chr22_KI270736v1_random
chr22_KI270737v1_random
chr22_KI270738v1_random
chr22_KI270739v1_random
chrUn_KI270322v1
chrUn_KI270316v1
chrUn_KI270315v1
chrUn_KI270317v1
chrUn_KI270442v1
chrUn_KI270466v1
chrUn_KI270465v1
chrUn_KI270467v1
chrUn_KI270435v1
chrUn_KI270438v1
chrUn_KI270468v1
chrUn_KI270510v1
chrUn_KI270509v1
chrUn_KI270518v1
chrUn_KI270508v1
chrUn_KI270516v1
chrUn_KI270512v1
chrUn_KI270519v1
chrUn_KI270522v1
chrUn_KI270511v1
chrUn_KI270515v1
chrUn_KI270507v1
chrUn_KI270517v1
chrUn_KI270529v1
chrUn_KI270538v1
chrUn_KI270583v1
chrUn_KI270587v1
chrUn_KI270580v1
chrUn_KI270589v1
chrUn_KI270590v1
chrUn_KI270584v1
chrUn_KI270582v1
chrUn_KI270588v1
chrUn_KI270591v1
chrUn_KI270330v1
chrUn_KI270333v1
chrUn_KI270337v1
chrUn_KI270362v1
chrUn_KI270448v1
chrUn_KI270521v1
chrUn_GL000195v1
chrUn_GL000219v1
chrUn_GL000220v1
chrUn_GL000224v1
chrUn_GL000213v1
chrUn_KI270743v1
chrUn_KI270744v1
chrUn_KI270745v1
chrUn_KI270746v1
chrUn_KI270748v1
chrUn_KI270749v1
chrUn_KI270750v1
chrUn_KI270751v1
chrUn_KI270753v1
chrUn_KI270754v1
chrUn_KI270756v1
chrUn_KI270757v1
chrUn_GL000214v1
chrUn_KI270742v1
chrUn_GL000216v2
chrUn_GL000218v1
chr1_KI270762v1_alt
chr1_KI270766v1_alt
chr1_KI270760v1_alt
chr1_KI270765v1_alt
chr1_KI270764v1_alt
chr2_KI270773v1_alt
chr2_KI270774v1_alt
chr2_KI270772v1_alt
chr2_KI270775v1_alt
chr2_KI270768v1_alt
chr2_GL383522v1_alt
chr3_KI270780v1_alt
chr3_GL383526v1_alt
chr3_KI270779v1_alt
chr3_KI270784v1_alt
chr4_GL000257v2_alt
chr4_KI270785v1_alt
chr4_KI270789v1_alt
chr5_KI270792v1_alt
chr5_GL383532v1_alt
chr5_GL949742v1_alt
chr5_GL339449v2_alt
chr5_GL383531v1_alt
chr5_KI270795v1_alt
chr7_KI270809v1_alt
chr7_KI270803v1_alt
chr7_KI270805v1_alt
chr8_KI270821v1_alt
chr8_KI270813v1_alt
chr8_KI270822v1_alt
chr8_KI270810v1_alt
chr8_KI270819v1_alt
chr8_KI270820v1_alt
chr8_KI270817v1_alt
chr9_GL383542v1_alt
chr10_GL383545v1_alt
chr10_GL383546v1_alt
chr11_KI270832v1_alt
chr11_KI270830v1_alt
chr11_KI270831v1_alt
chr11_JH159137v1_alt
chr11_KI270827v1_alt
chr12_GL877875v1_alt
chr12_KI270835v1_alt
chr12_GL383550v2_alt
chr13_KI270840v1_alt
chr13_KI270843v1_alt
chr13_KI270841v1_alt
chr13_KI270838v1_alt
chr14_KI270844v1_alt
chr14_KI270847v1_alt
chr14_KI270845v1_alt
chr15_KI270851v1_alt
chr15_KI270848v1_alt
chr15_KI270849v1_alt
chr15_GL383555v2_alt
chr15_KI270850v1_alt
chr16_KI270854v1_alt
chr16_KI270856v1_alt
chr16_KI270853v1_alt
chr16_GL383556v1_alt
chr17_GL383563v3_alt
chr17_KI270862v1_alt
chr17_KI270857v1_alt
chr17_JH159146v1_alt
chr17_JH159147v1_alt
chr17_GL000258v2_alt
chr17_KI270858v1_alt
chr17_KI270859v1_alt
chr17_KI270860v1_alt
chr18_KI270864v1_alt
chr18_GL383568v1_alt
chr18_KI270863v1_alt
chr19_KI270865v1_alt
chr19_GL383573v1_alt
chr19_GL383575v2_alt
chr19_GL383574v1_alt
chr19_KI270866v1_alt
chr20_KI270869v1_alt
chr22_KI270878v1_alt
chr22_KI270879v1_alt
chr22_KI270877v1_alt
chr22_GL383582v2_alt
chrX_KI270880v1_alt
chr2_KI270894v1_alt
chr3_KI270895v1_alt
chr4_KI270896v1_alt
chr5_KI270897v1_alt
chr5_KI270898v1_alt
chr6_GL000251v2_alt
chr7_KI270899v1_alt
chr8_KI270901v1_alt
chr8_KI270900v1_alt
chr11_KI270902v1_alt
chr11_KI270903v1_alt
chr15_KI270905v1_alt
chr17_KI270909v1_alt
chr17_KI270908v1_alt
chrX_KI270913v1_alt
chr3_KI270924v1_alt
chr4_KI270925v1_alt
chr6_GL000252v2_alt
chr8_KI270926v1_alt
chr11_KI270927v1_alt
chr22_KI270928v1_alt
chr6_GL000253v2_alt
chr3_KI270935v1_alt
chr6_GL000254v2_alt
chr6_GL000255v2_alt
chr6_GL000256v2_alt
chr19_GL949752v1_alt
chr6_KI270758v1_alt
chr19_KI270938v1_alt
chrEBV
chrUn_KN707607v1_decoy
chrUn_KN707608v1_decoy
chrUn_KN707609v1_decoy
chrUn_KN707610v1_decoy
chrUn_KN707612v1_decoy
chrUn_KN707613v1_decoy
chrUn_KN707626v1_decoy
chrUn_KN707631v1_decoy
chrUn_KN707636v1_decoy
chrUn_KN707638v1_decoy
chrUn_KN707639v1_decoy
chrUn_KN707643v1_decoy
chrUn_KN707645v1_decoy
chrUn_KN707646v1_decoy
chrUn_KN707647v1_decoy
chrUn_KN707649v1_decoy
chrUn_KN707650v1_decoy
chrUn_KN707651v1_decoy
chrUn_KN707652v1_decoy
chrUn_KN707653v1_decoy
chrUn_KN707658v1_decoy
chrUn_KN707659v1_decoy
chrUn_KN707660v1_decoy
chrUn_KN707661v1_decoy
chrUn_KN707662v1_decoy
chrUn_KN707668v1_decoy
chrUn_KN707686v1_decoy
chrUn_KN707687v1_decoy
chrUn_KN707688v1_decoy
chrUn_KN707691v1_decoy
chrUn_KN707693v1_decoy
chrUn_KN707696v1_decoy
chrUn_KN707702v1_decoy
chrUn_KN707715v1_decoy
chrUn_KN707720v1_decoy
chrUn_KN707728v1_decoy
chrUn_KN707733v1_decoy
chrUn_KN707734v1_decoy
chrUn_KN707736v1_decoy
chrUn_KN707740v1_decoy
chrUn_KN707748v1_decoy
chrUn_KN707749v1_decoy
chrUn_KN707750v1_decoy
chrUn_KN707765v1_decoy
chrUn_KN707775v1_decoy
chrUn_KN707779v1_decoy
chrUn_KN707780v1_decoy
chrUn_KN707782v1_decoy
chrUn_KN707783v1_decoy
chrUn_KN707784v1_decoy
chrUn_KN707798v1_decoy
chrUn_KN707807v1_decoy
chrUn_KN707811v1_decoy
chrUn_KN707816v1_decoy
chrUn_KN707822v1_decoy
chrUn_KN707826v1_decoy
chrUn_KN707828v1_decoy
chrUn_KN707834v1_decoy
chrUn_KN707857v1_decoy
chrUn_KN707860v1_decoy
chrUn_KN707862v1_decoy
chrUn_KN707863v1_decoy
chrUn_KN707864v1_decoy
chrUn_KN707865v1_decoy
chrUn_KN707866v1_decoy
chrUn_KN707867v1_decoy
chrUn_KN707872v1_decoy
chrUn_KN707873v1_decoy
chrUn_KN707874v1_decoy
chrUn_KN707875v1_decoy
chrUn_KN707878v1_decoy
chrUn_KN707879v1_decoy
chrUn_KN707881v1_decoy
chrUn_KN707882v1_decoy
chrUn_KN707883v1_decoy
chrUn_KN707884v1_decoy
chrUn_KN707885v1_decoy
chrUn_KN707886v1_decoy
chrUn_KN707887v1_decoy
chrUn_KN707888v1_decoy
chrUn_KN707889v1_decoy
chrUn_KN707890v1_decoy
chrUn_KN707891v1_decoy
chrUn_KN707892v1_decoy
chrUn_KN707896v1_decoy
chrUn_KN707901v1_decoy
chrUn_KN707902v1_decoy
chrUn_KN707903v1_decoy
chrUn_KN707904v1_decoy
chrUn_KN707905v1_decoy
chrUn_KN707906v1_decoy
chrUn_KN707911v1_decoy
chrUn_KN707915v1_decoy
chrUn_KN707919v1_decoy
chrUn_KN707920v1_decoy
chrUn_KN707947v1_decoy
chrUn_KN707953v1_decoy
chrUn_KN707954v1_decoy
chrUn_KN707957v1_decoy
chrUn_KN707959v1_decoy
chrUn_KN707961v1_decoy
chrUn_KN707962v1_decoy
chrUn_KN707963v1_decoy
chrUn_KN707964v1_decoy
chrUn_KN707965v1_decoy
chrUn_KN707966v1_decoy
chrUn_KN707967v1_decoy
chrUn_KN707968v1_decoy
chrUn_KN707969v1_decoy
chrUn_KN707970v1_decoy
chrUn_KN707971v1_decoy
chrUn_KN707972v1_decoy
chrUn_KN707973v1_decoy
chrUn_KN707974v1_decoy
chrUn_KN707975v1_decoy
chrUn_KN707979v1_decoy
chrUn_KN707980v1_decoy
chrUn_KN707981v1_decoy
chrUn_KN707982v1_decoy
chrUn_KN707985v1_decoy
chrUn_KN707986v1_decoy
chrUn_KN707987v1_decoy
chrUn_KN707988v1_decoy
chrUn_KN707990v1_decoy
chrUn_JTFH01000001v1_decoy
chrUn_JTFH01000002v1_decoy
chrUn_JTFH01000009v1_decoy
chrUn_JTFH01000012v1_decoy
chrUn_JTFH01000013v1_decoy
chrUn_JTFH01000017v1_decoy
chrUn_JTFH01000018v1_decoy
chrUn_JTFH01000020v1_decoy
chrUn_JTFH01000022v1_decoy
chrUn_JTFH01000023v1_decoy
chrUn_JTFH01000029v1_decoy
chrUn_JTFH01000030v1_decoy
chrUn_JTFH01000031v1_decoy
chrUn_JTFH01000040v1_decoy
chrUn_JTFH01000041v1_decoy
chrUn_JTFH01000050v1_decoy
chrUn_JTFH01000054v1_decoy
chrUn_JTFH01000064v1_decoy
chrUn_JTFH01000068v1_decoy
chrUn_JTFH01000070v1_decoy
chrUn_JTFH01000071v1_decoy
chrUn_JTFH01000072v1_decoy
chrUn_JTFH01000074v1_decoy
chrUn_JTFH01000077v1_decoy
chrUn_JTFH01000082v1_decoy
chrUn_JTFH01000089v1_decoy
chrUn_JTFH01000090v1_decoy
chrUn_JTFH01000093v1_decoy
chrUn_JTFH01000096v1_decoy
chrUn_JTFH01000098v1_decoy
chrUn_JTFH01000099v1_decoy
chrUn_JTFH01000101v1_decoy
chrUn_JTFH01000106v1_decoy
chrUn_JTFH01000107v1_decoy
chrUn_JTFH01000111v1_decoy
chrUn_JTFH01000112v1_decoy
chrUn_JTFH01000115v1_decoy
chrUn_JTFH01000116v1_decoy
chrUn_JTFH01000117v1_decoy
chrUn_JTFH01000119v1_decoy
chrUn_JTFH01000123v1_decoy
chrUn_JTFH01000124v1_decoy
chrUn_JTFH01000126v1_decoy
chrUn_JTFH01000127v1_decoy
chrUn_JTFH01000129v1_decoy
chrUn_JTFH01000132v1_decoy
chrUn_JTFH01000133v1_decoy
chrUn_JTFH01000134v1_decoy
chrUn_JTFH01000135v1_decoy
chrUn_JTFH01000136v1_decoy
chrUn_JTFH01000140v1_decoy
chrUn_JTFH01000144v1_decoy
chrUn_JTFH01000145v1_decoy
chrUn_JTFH01000148v1_decoy
chrUn_JTFH01000150v1_decoy
chrUn_JTFH01000151v1_decoy
chrUn_JTFH01000153v1_decoy
chrUn_JTFH01000154v1_decoy
chrUn_JTFH01000155v1_decoy
chrUn_JTFH01000158v1_decoy
chrUn_JTFH01000159v1_decoy
chrUn_JTFH01000164v1_decoy
chrUn_JTFH01000167v1_decoy
chrUn_JTFH01000178v1_decoy
chrUn_JTFH01000179v1_decoy
chrUn_JTFH01000180v1_decoy
chrUn_JTFH01000181v1_decoy
chrUn_JTFH01000182v1_decoy
chrUn_JTFH01000185v1_decoy
chrUn_JTFH01000190v1_decoy
chrUn_JTFH01000191v1_decoy
chrUn_JTFH01000192v1_decoy
chrUn_JTFH01000194v1_decoy
chrUn_JTFH01000196v1_decoy
chrUn_JTFH01000200v1_decoy
chrUn_JTFH01000202v1_decoy
chrUn_JTFH01000205v1_decoy
chrUn_JTFH01000206v1_decoy
chrUn_JTFH01000207v1_decoy
chrUn_JTFH01000212v1_decoy
chrUn_JTFH01000214v1_decoy
chrUn_JTFH01000216v1_decoy
chrUn_JTFH01000217v1_decoy
chrUn_JTFH01000218v1_decoy
chrUn_JTFH01000223v1_decoy
chrUn_JTFH01000226v1_decoy
chrUn_JTFH01000230v1_decoy
chrUn_JTFH01000241v1_decoy
chrUn_JTFH01000242v1_decoy
chrUn_JTFH01000244v1_decoy
chrUn_JTFH01000246v1_decoy
chrUn_JTFH01000249v1_decoy
chrUn_JTFH01000257v1_decoy
chrUn_JTFH01000258v1_decoy
chrUn_JTFH01000263v1_decoy
chrUn_JTFH01000264v1_decoy
chrUn_JTFH01000266v1_decoy
chrUn_JTFH01000267v1_decoy
chrUn_JTFH01000268v1_decoy
chrUn_JTFH01000269v1_decoy
chrUn_JTFH01000272v1_decoy
chrUn_JTFH01000273v1_decoy
chrUn_JTFH01000274v1_decoy
chrUn_JTFH01000277v1_decoy
chrUn_JTFH01000280v1_decoy
chrUn_JTFH01000281v1_decoy
chrUn_JTFH01000285v1_decoy
chrUn_JTFH01000292v1_decoy
chrUn_JTFH01000295v1_decoy
chrUn_JTFH01000297v1_decoy
chrUn_JTFH01000302v1_decoy
chrUn_JTFH01000303v1_decoy
chrUn_JTFH01000306v1_decoy
chrUn_JTFH01000307v1_decoy
chrUn_JTFH01000315v1_decoy
chrUn_JTFH01000319v1_decoy
chrUn_JTFH01000323v1_decoy
chrUn_JTFH01000324v1_decoy
chrUn_JTFH01000326v1_decoy
chrUn_JTFH01000329v1_decoy
chrUn_JTFH01000340v1_decoy
chrUn_JTFH01000341v1_decoy
chrUn_JTFH01000342v1_decoy
chrUn_JTFH01000343v1_decoy
chrUn_JTFH01000346v1_decoy
chrUn_JTFH01000347v1_decoy
chrUn_JTFH01000348v1_decoy
chrUn_JTFH01000349v1_decoy
chrUn_JTFH01000351v1_decoy
chrUn_JTFH01000361v1_decoy
chrUn_JTFH01000366v1_decoy
chrUn_JTFH01000378v1_decoy
chrUn_JTFH01000380v1_decoy
chrUn_JTFH01000383v1_decoy
chrUn_JTFH01000384v1_decoy
chrUn_JTFH01000392v1_decoy
chrUn_JTFH01000395v1_decoy
chrUn_JTFH01000396v1_decoy
chrUn_JTFH01000397v1_decoy
chrUn_JTFH01000402v1_decoy
chrUn_JTFH01000403v1_decoy
chrUn_JTFH01000407v1_decoy
chrUn_JTFH01000410v1_decoy
chrUn_JTFH01000412v1_decoy
chrUn_JTFH01000415v1_decoy
chrUn_JTFH01000417v1_decoy
chrUn_JTFH01000418v1_decoy
chrUn_JTFH01000419v1_decoy
chrUn_JTFH01000420v1_decoy
chrUn_JTFH01000423v1_decoy
chrUn_JTFH01000426v1_decoy
chrUn_JTFH01000429v1_decoy
chrUn_JTFH01000430v1_decoy
chrUn_JTFH01000433v1_decoy
chrUn_JTFH01000434v1_decoy
chrUn_JTFH01000436v1_decoy
chrUn_JTFH01000443v1_decoy
chrUn_JTFH01000447v1_decoy
chrUn_JTFH01000450v1_decoy
chrUn_JTFH01000451v1_decoy
chrUn_JTFH01000456v1_decoy
chrUn_JTFH01000458v1_decoy
chrUn_JTFH01000459v1_decoy
chrUn_JTFH01000466v1_decoy
chrUn_JTFH01000480v1_decoy
chrUn_JTFH01000493v1_decoy
chrUn_JTFH01000495v1_decoy
chrUn_JTFH01000510v1_decoy
chrUn_JTFH01000513v1_decoy
chrUn_JTFH01000515v1_decoy
chrUn_JTFH01000517v1_decoy
chrUn_JTFH01000519v1_decoy
chrUn_JTFH01000526v1_decoy
chrUn_JTFH01000528v1_decoy
chrUn_JTFH01000530v1_decoy
chrUn_JTFH01000537v1_decoy
chrUn_JTFH01000544v1_decoy
chrUn_JTFH01000546v1_decoy
chrUn_JTFH01000550v1_decoy
chrUn_JTFH01000553v1_decoy
chrUn_JTFH01000554v1_decoy
chrUn_JTFH01000557v1_decoy
chrUn_JTFH01000564v1_decoy
chrUn_JTFH01000566v1_decoy
chrUn_JTFH01000567v1_decoy
chrUn_JTFH01000578v1_decoy
chrUn_JTFH01000589v1_decoy
chrUn_JTFH01000594v1_decoy
chrUn_JTFH01000600v1_decoy
chrUn_JTFH01000611v1_decoy
chrUn_JTFH01000616v1_decoy
chrUn_JTFH01000618v1_decoy
chrUn_JTFH01000619v1_decoy
chrUn_JTFH01000621v1_decoy
chrUn_JTFH01000627v1_decoy
chrUn_JTFH01000628v1_decoy
chrUn_JTFH01000631v1_decoy
chrUn_JTFH01000635v1_decoy
chrUn_JTFH01000638v1_decoy
chrUn_JTFH01000641v1_decoy
chrUn_JTFH01000645v1_decoy
chrUn_JTFH01000650v1_decoy
chrUn_JTFH01000651v1_decoy
chrUn_JTFH01000653v1_decoy
chrUn_JTFH01000657v1_decoy
chrUn_JTFH01000661v1_decoy
chrUn_JTFH01000663v1_decoy
chrUn_JTFH01000666v1_decoy
chrUn_JTFH01000667v1_decoy
chrUn_JTFH01000672v1_decoy
chrUn_JTFH01000675v1_decoy
chrUn_JTFH01000686v1_decoy
chrUn_JTFH01000700v1_decoy
chrUn_JTFH01000704v1_decoy
chrUn_JTFH01000707v1_decoy
chrUn_JTFH01000714v1_decoy
chrUn_JTFH01000715v1_decoy
chrUn_JTFH01000717v1_decoy
chrUn_JTFH01000730v1_decoy
chrUn_JTFH01000732v1_decoy
chrUn_JTFH01000740v1_decoy
chrUn_JTFH01000746v1_decoy
chrUn_JTFH01000749v1_decoy
chrUn_JTFH01000750v1_decoy
chrUn_JTFH01000759v1_decoy
chrUn_JTFH01000762v1_decoy
chrUn_JTFH01000784v1_decoy
chrUn_JTFH01000791v1_decoy
chrUn_JTFH01000795v1_decoy
chrUn_JTFH01000796v1_decoy
chrUn_JTFH01000797v1_decoy
chrUn_JTFH01000799v1_decoy
chrUn_JTFH01000800v1_decoy
chrUn_JTFH01000801v1_decoy
chrUn_JTFH01000802v1_decoy
chrUn_JTFH01000806v1_decoy
chrUn_JTFH01000808v1_decoy
chrUn_JTFH01000810v1_decoy
chrUn_JTFH01000820v1_decoy
chrUn_JTFH01000821v1_decoy
chrUn_JTFH01000824v1_decoy
chrUn_JTFH01000826v1_decoy
chrUn_JTFH01000829v1_decoy
chrUn_JTFH01000843v1_decoy
chrUn_JTFH01000844v1_decoy
chrUn_JTFH01000845v1_decoy
chrUn_JTFH01000846v1_decoy
chrUn_JTFH01000850v1_decoy
chrUn_JTFH01000852v1_decoy
chrUn_JTFH01000865v1_decoy
chrUn_JTFH01000866v1_decoy
chrUn_JTFH01000870v1_decoy
chrUn_JTFH01000872v1_decoy
chrUn_JTFH01000875v1_decoy
chrUn_JTFH01000876v1_decoy
chrUn_JTFH01000879v1_decoy
chrUn_JTFH01000896v1_decoy
chrUn_JTFH01000899v1_decoy
chrUn_JTFH01000900v1_decoy
chrUn_JTFH01000902v1_decoy
chrUn_JTFH01000906v1_decoy
chrUn_JTFH01000910v1_decoy
chrUn_JTFH01000920v1_decoy
chrUn_JTFH01000923v1_decoy
chrUn_JTFH01000927v1_decoy
chrUn_JTFH01000928v1_decoy
chrUn_JTFH01000937v1_decoy
chrUn_JTFH01000946v1_decoy
chrUn_JTFH01000948v1_decoy
chrUn_JTFH01000952v1_decoy
chrUn_JTFH01000961v1_decoy
chrUn_JTFH01000963v1_decoy
chrUn_JTFH01000968v1_decoy
chrUn_JTFH01000971v1_decoy
chrUn_JTFH01000972v1_decoy
chrUn_JTFH01000973v1_decoy
chrUn_JTFH01000975v1_decoy
chrUn_JTFH01000976v1_decoy
chrUn_JTFH01000977v1_decoy
chrUn_JTFH01000978v1_decoy
chrUn_JTFH01000979v1_decoy
chrUn_JTFH01000982v1_decoy
chrUn_JTFH01000983v1_decoy
chrUn_JTFH01000984v1_decoy
chrUn_JTFH01000985v1_decoy
chrUn_JTFH01000986v1_decoy
chrUn_JTFH01000987v1_decoy
chrUn_JTFH01000989v1_decoy
chrUn_JTFH01000992v1_decoy
chrUn_JTFH01000995v1_decoy
chrUn_JTFH01000996v1_decoy
chrUn_JTFH01000998v1_decoy
chrUn_JTFH01000999v1_decoy
chrUn_JTFH01001000v1_decoy
chrUn_JTFH01001002v1_decoy
chrUn_JTFH01001003v1_decoy
chrUn_JTFH01001006v1_decoy
chrUn_JTFH01001007v1_decoy
chrUn_JTFH01001008v1_decoy
chrUn_JTFH01001009v1_decoy
chrUn_JTFH01001011v1_decoy
chrUn_JTFH01001013v1_decoy
chrUn_JTFH01001017v1_decoy
chrUn_JTFH01001018v1_decoy
chrUn_JTFH01001020v1_decoy
chrUn_JTFH01001032v1_decoy
chrUn_JTFH01001034v1_decoy
chrUn_JTFH01001036v1_decoy
chrUn_JTFH01001037v1_decoy
chrUn_JTFH01001039v1_decoy
chrUn_JTFH01001043v1_decoy
chrUn_JTFH01001044v1_decoy
chrUn_JTFH01001045v1_decoy
chrUn_JTFH01001046v1_decoy
chrUn_JTFH01001047v1_decoy
chrUn_JTFH01001049v1_decoy
chrUn_JTFH01001054v1_decoy
chrUn_JTFH01001056v1_decoy
chrUn_JTFH01001058v1_decoy
chrUn_JTFH01001064v1_decoy
chrUn_JTFH01001067v1_decoy
chrUn_JTFH01001068v1_decoy
chrUn_JTFH01001070v1_decoy
chrUn_JTFH01001072v1_decoy
chrUn_JTFH01001073v1_decoy
chrUn_JTFH01001074v1_decoy
chrUn_JTFH01001075v1_decoy
chrUn_JTFH01001087v1_decoy
chrUn_JTFH01001089v1_decoy
chrUn_JTFH01001097v1_decoy
chrUn_JTFH01001098v1_decoy
chrUn_JTFH01001101v1_decoy
chrUn_JTFH01001108v1_decoy
chrUn_JTFH01001109v1_decoy
chrUn_JTFH01001113v1_decoy
chrUn_JTFH01001114v1_decoy
chrUn_JTFH01001126v1_decoy
chrUn_JTFH01001142v1_decoy
chrUn_JTFH01001144v1_decoy
chrUn_JTFH01001145v1_decoy
chrUn_JTFH01001147v1_decoy
chrUn_JTFH01001148v1_decoy
chrUn_JTFH01001153v1_decoy
chrUn_JTFH01001155v1_decoy
chrUn_JTFH01001156v1_decoy
chrUn_JTFH01001161v1_decoy
chrUn_JTFH01001164v1_decoy
chrUn_JTFH01001167v1_decoy
chrUn_JTFH01001168v1_decoy
chrUn_JTFH01001172v1_decoy
chrUn_JTFH01001176v1_decoy
chrUn_JTFH01001180v1_decoy
chrUn_JTFH01001184v1_decoy
chrUn_JTFH01001188v1_decoy
chrUn_JTFH01001189v1_decoy
chrUn_JTFH01001193v1_decoy
chrUn_JTFH01001194v1_decoy
chrUn_JTFH01001204v1_decoy
chrUn_JTFH01001206v1_decoy
chrUn_JTFH01001207v1_decoy
chrUn_JTFH01001208v1_decoy
chrUn_JTFH01001212v1_decoy
chrUn_JTFH01001217v1_decoy
chrUn_JTFH01001223v1_decoy
chrUn_JTFH01001226v1_decoy
chrUn_JTFH01001228v1_decoy
chrUn_JTFH01001235v1_decoy
chrUn_JTFH01001239v1_decoy
chrUn_JTFH01001241v1_decoy
chrUn_JTFH01001242v1_decoy
chrUn_JTFH01001243v1_decoy
chrUn_JTFH01001246v1_decoy
chrUn_JTFH01001252v1_decoy
chrUn_JTFH01001255v1_decoy
chrUn_JTFH01001257v1_decoy
chrUn_JTFH01001269v1_decoy
chrUn_JTFH01001271v1_decoy
chrUn_JTFH01001282v1_decoy
chrUn_JTFH01001284v1_decoy
chrUn_JTFH01001297v1_decoy
chrUn_JTFH01001305v1_decoy
chrUn_JTFH01001317v1_decoy
chrUn_JTFH01001319v1_decoy
chrUn_JTFH01001323v1_decoy
chrUn_JTFH01001328v1_decoy
chrUn_JTFH01001329v1_decoy
chrUn_JTFH01001332v1_decoy
chrUn_JTFH01001336v1_decoy
chrUn_JTFH01001337v1_decoy
chrUn_JTFH01001345v1_decoy
chrUn_JTFH01001351v1_decoy
chrUn_JTFH01001353v1_decoy
chrUn_JTFH01001357v1_decoy
chrUn_JTFH01001365v1_decoy
chrUn_JTFH01001368v1_decoy
chrUn_JTFH01001372v1_decoy
chrUn_JTFH01001374v1_decoy
chrUn_JTFH01001375v1_decoy
chrUn_JTFH01001377v1_decoy
chrUn_JTFH01001386v1_decoy
chrUn_JTFH01001387v1_decoy
chrUn_JTFH01001388v1_decoy
chrUn_JTFH01001390v1_decoy
chrUn_JTFH01001391v1_decoy
chrUn_JTFH01001394v1_decoy
chrUn_JTFH01001395v1_decoy
chrUn_JTFH01001398v1_decoy
chrUn_JTFH01001405v1_decoy
chrUn_JTFH01001416v1_decoy
chrUn_JTFH01001424v1_decoy
chrUn_JTFH01001425v1_decoy
chrUn_JTFH01001430v1_decoy
chrUn_JTFH01001432v1_decoy
chrUn_JTFH01001433v1_decoy
chrUn_JTFH01001435v1_decoy
chrUn_JTFH01001446v1_decoy
chrUn_JTFH01001454v1_decoy
chrUn_JTFH01001465v1_decoy
chrUn_JTFH01001478v1_decoy
chrUn_JTFH01001506v1_decoy
chrUn_JTFH01001512v1_decoy
chrUn_JTFH01001513v1_decoy
chrUn_JTFH01001514v1_decoy
chrUn_JTFH01001520v1_decoy
chrUn_JTFH01001539v1_decoy
chrUn_JTFH01001541v1_decoy
chrUn_JTFH01001549v1_decoy
chrUn_JTFH01001553v1_decoy
chrUn_JTFH01001556v1_decoy
chrUn_JTFH01001557v1_decoy
chrUn_JTFH01001565v1_decoy
chrUn_JTFH01001574v1_decoy
chrUn_JTFH01001590v1_decoy
chrUn_JTFH01001596v1_decoy
chrUn_JTFH01001608v1_decoy
chrUn_JTFH01001613v1_decoy
chrUn_JTFH01001615v1_decoy
chrUn_JTFH01001632v1_decoy
chrUn_JTFH01001664v1_decoy
chrUn_JTFH01001669v1_decoy
chrUn_JTFH01001677v1_decoy
chrUn_JTFH01001680v1_decoy
chrUn_JTFH01001701v1_decoy
chrUn_JTFH01001724v1_decoy
chrUn_JTFH01001738v1_decoy
chrUn_JTFH01001739v1_decoy
chrUn_JTFH01001740v1_decoy
chrUn_JTFH01001746v1_decoy
chrUn_JTFH01001748v1_decoy
chrUn_JTFH01001750v1_decoy
chrUn_JTFH01001752v1_decoy
chrUn_JTFH01001753v1_decoy
chrUn_JTFH01001754v1_decoy
chrUn_JTFH01001755v1_decoy
chrUn_JTFH01001757v1_decoy
chrUn_JTFH01001758v1_decoy
chrUn_JTFH01001760v1_decoy
chrUn_JTFH01001764v1_decoy
chrUn_JTFH01001768v1_decoy
chrUn_JTFH01001776v1_decoy
chrUn_JTFH01001780v1_decoy
chrUn_JTFH01001781v1_decoy
chrUn_JTFH01001783v1_decoy
chrUn_JTFH01001785v1_decoy
chrUn_JTFH01001793v1_decoy
chrUn_JTFH01001794v1_decoy
chrUn_JTFH01001797v1_decoy
chrUn_JTFH01001800v1_decoy
chrUn_JTFH01001805v1_decoy
chrUn_JTFH01001807v1_decoy
chrUn_JTFH01001820v1_decoy
chrUn_JTFH01001822v1_decoy
chrUn_JTFH01001823v1_decoy
chrUn_JTFH01001824v1_decoy
chrUn_JTFH01001829v1_decoy
chrUn_JTFH01001831v1_decoy
chrUn_JTFH01001835v1_decoy
chrUn_JTFH01001838v1_decoy
chrUn_JTFH01001844v1_decoy
chrUn_JTFH01001847v1_decoy
chrUn_JTFH01001849v1_decoy
chrUn_JTFH01001850v1_decoy
chrUn_JTFH01001858v1_decoy
chrUn_JTFH01001862v1_decoy
chrUn_JTFH01001866v1_decoy
chrUn_JTFH01001867v1_decoy
chrUn_JTFH01001868v1_decoy
chrUn_JTFH01001869v1_decoy
chrUn_JTFH01001870v1_decoy
chrUn_JTFH01001875v1_decoy
chrUn_JTFH01001876v1_decoy
chrUn_JTFH01001877v1_decoy
chrUn_JTFH01001878v1_decoy
chrUn_JTFH01001884v1_decoy
chrUn_JTFH01001885v1_decoy
chrUn_JTFH01001888v1_decoy
chrUn_JTFH01001889v1_decoy
chrUn_JTFH01001890v1_decoy
chrUn_JTFH01001893v1_decoy
chrUn_JTFH01001894v1_decoy
chrUn_JTFH01001896v1_decoy
chrUn_JTFH01001898v1_decoy
chrUn_JTFH01001899v1_decoy
chrUn_JTFH01001906v1_decoy
chrUn_JTFH01001907v1_decoy
chrUn_JTFH01001911v1_decoy
chrUn_JTFH01001914v1_decoy
chrUn_JTFH01001915v1_decoy
chrUn_JTFH01001916v1_decoy
chrUn_JTFH01001919v1_decoy
chrUn_JTFH01001920v1_decoy
chrUn_JTFH01001923v1_decoy
chrUn_JTFH01001925v1_decoy
chrUn_JTFH01001927v1_decoy
chrUn_JTFH01001929v1_decoy
chrUn_JTFH01001931v1_decoy
chrUn_JTFH01001934v1_decoy
chrUn_JTFH01001937v1_decoy
chrUn_JTFH01001939v1_decoy
chrUn_JTFH01001940v1_decoy
chrUn_JTFH01001941v1_decoy
chrUn_JTFH01001943v1_decoy
chrUn_JTFH01001944v1_decoy
chrUn_JTFH01001946v1_decoy
chrUn_JTFH01001947v1_decoy
chrUn_JTFH01001948v1_decoy
chrUn_JTFH01001949v1_decoy
chrUn_JTFH01001953v1_decoy
chrUn_JTFH01001956v1_decoy
chrUn_JTFH01001957v1_decoy
chrUn_JTFH01001959v1_decoy
chrUn_JTFH01001960v1_decoy
chrUn_JTFH01001961v1_decoy
chrUn_JTFH01001963v1_decoy
chrUn_JTFH01001966v1_decoy
chrUn_JTFH01001970v1_decoy
chrUn_JTFH01001972v1_decoy
chrUn_JTFH01001973v1_decoy
chrUn_JTFH01001974v1_decoy
chrUn_JTFH01001976v1_decoy
chrUn_JTFH01001978v1_decoy
chrUn_JTFH01001979v1_decoy
chrUn_JTFH01001980v1_decoy
chrUn_JTFH01001981v1_decoy
chrUn_JTFH01001984v1_decoy
chrUn_JTFH01001985v1_decoy
chrUn_JTFH01001986v1_decoy
chrUn_JTFH01001991v1_decoy
chrUn_JTFH01001992v1_decoy
chrUn_JTFH01001995v1_decoy
chrUn_JTFH01001997v1_decoy
chrUn_JTFH01001998v1_decoy
```


```bash
vcftools \
--gzvcf ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf.gz \
--gzdiff ./known/project.NIST.hc.snps.indels.vcf.gz \
--diff-site \
--out ./NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels \
--chr chr[1-22] \
--chr chrX \
--chr chrY \
--chr chrM
```

### [Nvidia parabricks germline](https://github.com/ESR-NZ/ESR-Parabricks)

See the results and settings of the pipeline runs on the Genome In A Bottle (GIAB) sample [NA12878](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) (NIST7035 and NIST7086) for:

- Nvidia parabricks germline at ...

Located at `/usr/local/bin/pbrun`

Run:

```bash

```

## Notes

- I found it difficult to use the hap.py wrapper for vcfeval since there were parameters I wasn't able to pass to vcfeval (such as specifying a sample in a multi sample vcf file and specifying the output mode).
