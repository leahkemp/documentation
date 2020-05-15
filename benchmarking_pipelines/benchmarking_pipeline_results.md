# Benchmarking results

## Table of contents

- [Benchmarking results](#benchmarking-results)
  - [Table of contents](#table-of-contents)
  - [Known vcf](#known-vcf)
  - [bench1.0](#bench10)
    - [Run parameters/settings](#run-parameterssettings)
    - [Results](#results)
  - [bench1.1](#bench11)

## Known vcf

Parameters used to create the Genome In A Bottle (GIAB) sample [NA12878](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) (both NIST7035 and NIST7086) known vcf file:

- **Inputs:**
  - Reference genome: ucsc.hg19.fasta
  - dbSNP database: dbsnp_137.hg19.vcf

## bench1.0

### Run parameters/settings

- **Aim:** Benchmarking against against the Genome In A Bottle (GIAB) sample [NA12878](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) (both NIST7035 and NIST7086)
- **Pipeline:** [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) and minimal [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline) (no annotation)
- **Inputs human_genomics_pipeline:**
  - Reference genome: ucsc.hg19.fasta
  - dbSNP database: dbsnp_138.hg19.vcf 
  - (see [run settings/output](https://github.com/ESR-NZ/human_genomics_pipeline/tree/bench1.0) for more detail)

- **Inputs vcf_annotation_pipeline:**
  - Reference genome: ucsc.hg19.fasta
  - dbSNP database: dbsnp_138.hg19.vcf
  - Hapmap: hapmap_3.3.hg19.sites.vcf.gz
  - Mills: Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
  - (see the [run settings/output](https://github.com/ESR-NZ/vcf_annotation_pipeline/tree/bench1.0) for more detail)

- **Other settings:**
  - No padding
  - No intervals

*Note. dbsnp_137.hg19.vcf doesn't appear to be available anymore, therefore I used the closest version available*

### Results

## bench1.1
