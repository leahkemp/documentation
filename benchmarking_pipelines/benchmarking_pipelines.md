# Benchmarking genomic pipelines

Created: 2020-04-22 13:37:04
Last modified: 2020/05/12 18:16:04

- **Aim:** Undertake benchmarking of genomics pipelines to test their quality for clinical use. 
- **Prerequisite software:** [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html), [bgzip](http://www.htslib.org/doc/bgzip.html), [tabix](http://www.htslib.org/doc/tabix.html)
- **OS:** Ubuntu 16.04 (Wintermute - research server)

The pipelines were run against the Genome In A Bottle (GIAB) sample [NA12878](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) and the quality of their outputs evaluated. See the complementary docs for benchmarking the Nvidia Parabricks pipeline [here](https://github.com/ESR-NZ/ESR-Parabricks).

## Table of contents

- [Benchmarking genomic pipelines](#benchmarking-genomic-pipelines)
  - [Table of contents](#table-of-contents)
  - [Setup](#setup)
    - [Install benchmarking software](#install-benchmarking-software)
      - [hap.py and [RTG Tools](https://github.com/RealTimeGenomics/rtg-tools/tree/eb13bbb82d2fbeab7d54a92e8493ddd2acf0d349)](#happy-and-rtg-tools)
    - [Download and prepare data](#download-and-prepare-data)
  - [Benchmarking](#benchmarking)
    - [human_genomics_pipeline](#humangenomicspipeline)
      - [Compare the truth and query vcf](#compare-the-truth-and-query-vcf)
    - [Nvidia parabricks germline](#nvidia-parabricks-germline)

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

Truth vcf for [NA12878](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) sample was accessed and prepared how described [here]([here](https://github.com/ESR-NZ/ESR-Parabricks)).

## Benchmarking

See [this paper](https://www.nature.com/articles/s41587-019-0054-x) for best practices for benchmarking germline small-variant calls in human genomes

### [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline)

See the results of the human_genomics_pipeline run on the Genome In A Bottle (GIAB) sample [NA12878](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) at ....

#### Compare the truth and query vcf

"Comparing variants at the level of the genomic haplotypes that the variants represent as a way to overcome the problems associated with comparing complex variants, in which alternative yet equivalent variant representations can confound direct comparison methods"

The known sites file and query vcf file needs to be bgzipped and have a tabix index file (.tbi) (write to new files so as to not modify the original files)

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/
# Query vcf
bgzip < ./human_genomics_pipeline/vcf/NA12878_NIST_raw_snps_indels_AS_g.vcf > ./human_genomics_pipeline/vcf/NA12878_NIST_raw_snps_indels_AS_g.vcf.gz
tabix ./human_genomics_pipeline/vcf/NA12878_NIST_raw_snps_indels_AS_g.vcf.gz
# Known vcf
bgzip < ./known/project.NIST.hc.snps.indels.vcf > ./known/project.NIST.hc.snps.indels.vcf.gz
tabix ./known/project.NIST.hc.snps.indels.vcf.gz
```

We also need to create an sdf file for the reference human genome (that was used in the benchmarking pipeline run). Create this with the rgt-tools format function

```bash
./hap.py-install/libexec/rtg-tools-install/rtg format \
--output ./hap.py-install/libexec/rtg-tools-install/Homo_sapiens_assembly38.fasta.sdf \
/store/lkemp/publicData/referenceGenome/gatkBundle/GRCh38/Homo_sapiens_assembly38.fasta
```

*Note. the two columns in the vcf file (NIST7035 and NIST7086) represent 2 vials of the same sample (NA12878)*

```bash
./hap.py-install/bin/hap.py \
./known/project.NIST.hc.snps.indels.vcf \
./human_genomics_pipeline/vcf/NA12878_NIST_raw_snps_indels_AS_g.vcf.gz \
-f ./known/nexterarapidcapture_expandedexome_targetedregions.bed.gz \
-o happy_human_genomics_pipeline \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh38/Homo_sapiens_assembly38.fasta \
--engine=vcfeval \
--engine-vcfeval-template ./hap.py-install/libexec/rtg-tools-install/Homo_sapiens_assembly38.fasta.sdf \
--threads 12
```

STUCK FROM HERE: need to tell it which columns to use, ie. pass to vcfeval the `--sample NIST7035,20` argument

This will generate VCF files containing called variants that were in the truth VCF (tp), called variants that were not in the truth VCF (fp) and truth variants that were not in the called variants (fn) (for a more in depth explanation see [here](https://github.com/ga4gh/benchmarking-tools/blob/master/doc/ref-impl/intermediate.md)).

### [Nvidia parabricks germline](https://github.com/ESR-NZ/ESR-Parabricks)

See the results of the Nvidia parabricks run on the Genome In A Bottle (GIAB) sample [NA12878](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) at ....
