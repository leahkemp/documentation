# Benchmarking genomic pipelines

Created: 2020-04-22 13:37:04
Last modified: 2020/04/23 17:32:42

- **Aim:** Undertake benchmarking of genomics pipelines to test their quality for clinical use
- **Prerequisite software:** [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html), [bgzip](http://www.htslib.org/doc/bgzip.html), [tabix](http://www.htslib.org/doc/tabix.html)
- **OS:** Ubuntu 16.04 (Wintermute - research server)

## Table of contents

- [Benchmarking genomic pipelines](#benchmarking-genomic-pipelines)
  - [Table of contents](#table-of-contents)
  - [Setup](#setup)
    - [Install software](#install-software)
      - [RTG Tools (for vcfeval)](#rtg-tools-for-vcfeval)
      - [hap.py](#happy)
    - [Download data](#download-data)
  - [Benchmarking](#benchmarking)
    - [human_genomics_pipeline](#humangenomicspipeline)
      - [Compare the truth and query vcf](#compare-the-truth-and-query-vcf)
      - [Compare ...](#compare)
    - [Nvidia parabricks](#nvidia-parabricks)

## Setup

### Install software

#### [RTG Tools (for vcfeval)](https://github.com/RealTimeGenomics/rtg-tools/tree/eb13bbb82d2fbeab7d54a92e8493ddd2acf0d349)

This installation requires ant and Java 1.8 or later

```bash
git clone https://github.com/RealTimeGenomics/rtg-tools.git
cd rtg-tools
```

Compile/run tests

```bash
ant runalltests
```

Build package

```bash
ant zip-nojre
```

Install

```bash
cd /store/lkemp/exome_project/benchmarking/
unzip /store/lkemp/exome_project/benchmarking/rtg-tools/build/rtg-tools-3.11-39691f9-base.zip
```

#### hap.py

```bash
conda create -n benchmarking_env python=2.7.5
conda activate benchmarking_env
```

(Need python 2 in order to get htslib2 that is required for installation)
Install [hap.py](https://github.com/Illumina/hap.py) dependencies

```bash
conda install -c conda-forge cython=0.29.15
conda install -c conda-forge scipy=1.2.1 # Also installs numpy=1.16.5 dependency
conda install -c conda-forge pandas=0.24.2
conda install -c bioconda pybedtools=0.8.1 # Also installs pysam=0.15.3 dependency
conda install -c bioconda bx-python=0.8.8
```

Install [hap.py](https://github.com/Illumina/hap.py) (using the python and the helper script)

```bash
git clone git@github.com:Illumina/hap.py.git
cd hap.py
python install.py ~/hap.py-install --no-tests
```

### Download data

## Benchmarking

See [this paper](https://www.nature.com/articles/s41587-019-0054-x) for best practices for benchmarking germline small-variant calls in human genomes

### [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline)

#### Compare the truth and query vcf

Create copies of all files into one location since their structure will be manipulated

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/
# Query vcf
cp /store/lkemp/exome_project/benchmarking/NA12878_exome/human_genomics_pipeline/vcf/NA12878_NIST.raw.snps.indels.AS.g.vcf .
# Known vcf
cp /store/lkemp/exome_project/benchmarking/known/project.NIST.hc.snps.indels.vcf .
# Reference human genome and associated files
cp /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh38/* .
```

The known sites file and query vcf file needs to be bgzipped and have a tabix index file (.tbi)

```bash
bgzip project.NIST.hc.snps.indels.vcf
tabix project.NIST.hc.snps.indels.vcf.gz
bgzip NA12878_NIST.raw.snps.indels.AS.g.vcf
tabix NA12878_NIST.raw.snps.indels.AS.g.vcf.gz
```

NEEDS AN SDF FILE?

Get SDF by formatting the reference genome

```bash
/store/lkemp/exome_project/benchmarking/rtg-tools-3.11-39691f9/rtg format \
-I Homo_sapiens_assembly38.fasta \
-o Homo_sapiens_assembly38.sdf.fasta
```

With vcfeval

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/
/store/lkemp/exome_project/benchmarking/rtg-tools-3.11-39691f9/rtg vcfeval \
-b project.NIST.hc.snps.indels.vcf.gz \
-c NA12878_NIST.raw.snps.indels.AS.g.vcf.gz \
-t reference.txt \
-o vcfeval_human_genomics_pipeline
```

#### Compare ...

With [hap.py](https://github.com/Illumina/hap.py)

```bash
python /home/lkemp/hap.py-install/bin/hap.py \
/store/lkemp/exome_project/benchmarking/known/project.NIST.hc.snps.indels.vcf \
/store/lkemp/exome_project/benchmarking/NA12878_exome/human_genomics_pipeline/vcf/NA12878_NIST.raw.snps.indels.AS.g.vcf \
-f /store/lkemp/exome_project/benchmarking/known/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-o happyOutput \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh38/Homo_sapiens_assembly38.fasta
```

Error: "too many AD fields"

See [this thread](https://github.com/Illumina/hap.py/issues/86)

- Look at the vcf file output by human_genomics_pipeline

```bash
cat NA12878_NIST.raw.snps.indels.AS.g.vcf | grep -v '^##' | head -n 50
```

```bash
#CHROM  POS     ID      REF     ALT             QUAL    FILTER  INFO                                                                   FORMAT                  20
chr1    1       .       N       <NON_REF>       .       .       END=10353                                                              GT:DP:GQ:MIN_DP:PL      0/0:0:0:0:0,0,0
chr1    10354   .       C       <NON_REF>       .       .       END=10357                                                              GT:DP:GQ:MIN_DP:PL      0/0:1:3:1:0,3,30
chr1    10358   .       A       <NON_REF>       .       .       END=10378                                                              GT:DP:GQ:MIN_DP:PL      0/0:2:6:2:0,6,38
chr1    10379   .       C       <NON_REF>       .       .       END=10379                                                              GT:DP:GQ:MIN_DP:PL      0/0:0:0:0:0,0,0
chr1    10380   .       C       <NON_REF>       .       .       END=10383                                                              GT:DP:GQ:MIN_DP:PL      0/0:2:6:2:0,6,47
chr1    10384   .       C       <NON_REF>       .       .       END=10431                                                              GT:DP:GQ:MIN_DP:PL      0/0:1:3:1:0,3,14
chr1    10432   .       A       <NON_REF>       .       .       END=12915                                                              GT:DP:GQ:MIN_DP:PL      0/0:0:0:0:0,0,0
chr1    12916   .       T       <NON_REF>       .       .       END=13005                                                              GT:DP:GQ:MIN_DP:PL      0/0:1:3:1:0,3,21
chr1    13006   .       G       <NON_REF>       .       .       END=13057                                                              GT:DP:GQ:MIN_DP:PL      0/0:0:0:0:0,0,0
chr1    13058   .       C       <NON_REF>       .       .       END=13130                                                              GT:DP:GQ:MIN_DP:PL      0/0:1:3:1:0,3,20
chr1    13131   .       T       <NON_REF>       .       .       END=13131                                                              GT:DP:GQ:MIN_DP:PL      0/0:1:0:1:0,0,0
chr1    13132   .       G       <NON_REF>       .       .       END=13148                                                              GT:DP:GQ:MIN_DP:PL      0/0:1:3:1:0,3,29
chr1    13149   .       G       <NON_REF>       .       .       END=13585                                                              GT:DP:GQ:MIN_DP:PL      0/0:0:0:0:0,0,0
chr1    13586   .       A       <NON_REF>       .       .       END=13595                                                              GT:DP:GQ:MIN_DP:PL      0/0:2:6:2:0,6,65
chr1    13596   .       G       <NON_REF>       .       .       END=13612                                                              GT:DP:GQ:MIN_DP:PL      0/0:4:12:4:0,12,126
chr1    13613   .       T       <NON_REF>       .       .       END=13613                                                              GT:DP:GQ:MIN_DP:PL      0/0:4:0:4:0,0,64
chr1    13614   .       T       <NON_REF>       .       .       END=13615                                                              GT:DP:GQ:MIN_DP:PL      0/0:4:12:4:0,12,134
chr1    57158   .       G       C,<NON_REF>     37.31   .       DP=2;ExcessHet=3.0103;MLEAC=1,0;MLEAF=0.500,0.00;RAW_MQandDP=1250,2    GT:AD:DP:GQ:PL:SB       1/1:0,2,0:2:6:49,6,0,49,6,49:0,0,1,1
```

Field descriptions

```bash
cat NA12878_NIST.raw.snps.indels.AS.g.vcf | grep '^##'
```

```bash
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele not already represented at this location by REF and ALT">
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another; will always be heterozygous and is not intended to describe called alleles">
##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phasing set (typically the position of the first variant in the set)">
##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.">

##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=ExcessHet,Number=1,Type=Float,Description="Phred-scaled p-value for exact test of excess heterozygosity">
##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
##INFO=<ID=RAW_MQandDP,Number=2,Type=Integer,Description="Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
```

- Look at the vcf file of the known variants

```bash
cat project.NIST.hc.snps.indels.vcf | grep -v '^##' | head -n 50
```

```bash
#CHROM  POS     ID              REF     ALT     QUAL            FILTER  INFO                                                                                                                                                               FORMAT          NIST7035                        NIST7086
chrM    73      .               G       A       7451.20         .       AC=4;AF=1.00;AN=4;DP=250;FS=0.000;MLEAC=4;MLEAF=1.00;MQ=53.56;MQ0=0;QD=29.80                                                                                       GT:AD:DP:GQ:PL  1/1:0,122:122:99:3604,369,0     1/1:0,128:128:99:3874,393,0
chrM    150     .               T       C       14828.20        .       AC=4;AF=1.00;AN=4;DP=325;FS=0.000;MLEAC=4;MLEAF=1.00;MQ=57.01;MQ0=0;QD=27.70                                                                                       GT:AD:DP:GQ:PL  1/1:0,166:166:99:7622,515,0     1/1:0,159:159:99:7233,488,0
chrM    152     rs117135796     T       C       14828.20        .       AC=4;AF=1.00;AN=4;DB;DP=329;FS=0.000;MLEAC=4;MLEAF=1.00;MQ=57.10;MQ0=0;QD=24.92                                                                                    GT:AD:DP:GQ:PL  1/1:0,169:169:99:7622,515,0     1/1:0,160:160:99:7233,488,0
chrM    195     .               C       T       12097.20        .       AC=4;AF=1.00;AN=4;DP=369;FS=0.000;MLEAC=4;MLEAF=1.00;MQ=58.11;MQ0=0;QD=32.78                                                                                       GT:AD:DP:GQ:PL  1/1:0,187:187:99:6223,570,0     1/1:0,182:182:99:5901,547,0
chrM    302     .               A       AC      355.40          .       AC=2;AF=0.500;AN=4;BaseQRankSum=-1.005;ClippingRankSum=-1.368;DP=273;FS=14.323;MLEAC=2;MLEAF=0.500;MQ=58.90;MQ0=0;MQRankSum=1.775;QD=1.30;ReadPosRankSum=-0.051    GT:AD:DP:GQ:PL  0/1:89,23:112:99:142,0,2214     0/1:96,24:120:99:251,0,2232
chrM    410     .               A       T       10654.20        .       AC=4;AF=1.00;AN=4;DP=322;FS=0.000;MLEAC=4;MLEAF=1.00;MQ=58.30;MQ0=0;QD=33.09                                                                                       GT:AD:DP:GQ:PL  1/1:0,163:163:99:5355,489,0     1/1:0,159:159:99:5326,478,0
chrM    2261    .               C       T       10148.20        .       AC=4;AF=1.00;AN=4;BaseQRankSum=0.177;ClippingRankSum=0.108;DP=345;FS=2.374;MLEAC=4;MLEAF=1.00;MQ=56.67;MQ0=0;MQRankSum=1.695;QD=29.42;ReadPosRankSum=-0.532        GT:AD:DP:GQ:PL  1/1:3,170:173:99:5002,434,0     1/1:0,172:172:99:5173,517,0
chrM    2354    .               C       T       7843.20         .       AC=4;AF=1.00;AN=4;BaseQRankSum=0.823;ClippingRankSum=1.008;DP=244;FS=0.000;MLEAC=4;MLEAF=1.00;MQ=57.37;MQ0=0;MQRankSum=0.170;QD=32.14;ReadPosRankSum=-0.568        GT:AD:DP:GQ:PL  1/1:1,107:108:99:3433,297,0     1/1:0,136:136:99:4437,409,0
chrM    2485    .               C       T       7774.20         .       AC=4;AF=1.00;AN=4;DP=231;FS=0.000;MLEAC=4;MLEAF=1.00;MQ=52.91;MQ0=0;QD=33.65                                                                                       GT:AD:DP:GQ:PL  1/1:0,102:102:99:3485,315,0     1/1:0,127:127:99:4316,393,0
chrM    2708    .               G       A       7088.20         .       AC=4;AF=1.00;AN=4;DP=217;FS=0.000;MLEAC=4;MLEAF=1.00;MQ=56.31;MQ0=0;QD=32.66                                                                                       GT:AD:DP:GQ:PL  1/1:0,106:106:99:3468,319,0     1/1:0,110:110:99:3647,331,0
```

Field decriptions

```bash
cat project.NIST.hc.snps.indels.vcf | grep '^##'
```

```bash
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">

##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
##INFO=<ID=ClippingRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP Membership">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description="Consistency of the site with at most two segregating haplotypes">
##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##INFO=<ID=MQ0,Number=1,Type=Integer,Description="Total Mapping Quality Zero Reads">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
```

### [Nvidia parabricks](https://github.com/ESR-NZ/ESR-Parabricks)