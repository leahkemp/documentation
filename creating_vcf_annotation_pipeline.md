# Create snakemake wrapper for vcf_annotation_pipeline

Created: 2020/03/13 11:00:43
Last modified: 2020/03/16 14:05:23

- **Aim:** convert a bash script pipeline written by Miles Benton into a snakemake pipeline: [vcf_annotation_pipeline](https://github.com/leahkemp/vcf_annotation_pipeline.git)
- **Prerequisite software:**
- **Prerequisite data:**
- **OS:** Ubuntu 16.04 (Wintermute - research server)

## VariantRecalibrator arguments

[See here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890851-Variant-annotations)

The annotations called by -an/--use-annotation are produced (in this case) by gatk4_HaplotypeCaller in human_genomics_pipeline

The


To set the parameters for the annotation argument (-n/--use-annotation) for the Variant Recalibrator step, see the input VCF file's INFO field for a list of all available annotations. Do this by printing the first section of one of the final vcf file output by human_genomics_pipeline:

```bash
head ../vcf_annnotation_pipeline/vcf/LJ_17BL1030_S8.raw.snps.indels.AS.g.vcf -n 200
```

Partial output:

```bash
#CHROM  POS     ID      REF     ALT             QUAL    FILTER  INFO                                                                 FORMAT                  20
chrM    1       .       G       <NON_REF>       .       .       END=43                                                               GT:DP:GQ:MIN_DP:PL      0/0:0:0:0:0,0,0
chrM    44      .       C       <NON_REF>       .       .       END=70                                                               GT:DP:GQ:MIN_DP:PL      0/0:1:3:1:0,3,16
chrM    71      .       G       <NON_REF>       .       .       END=72                                                               GT:DP:GQ:MIN_DP:PL      0/0:2:6:2:0,6,45
chrM    73      .       G       A,<NON_REF>     37.31   .       DP=2;ExcessHet=3.0103;MLEAC=1,0;MLEAF=0.500,0.00;RAW_MQandDP=7200,2  GT:AD:DP:GQ:PL:SB       1/1:0,2,0:2:6:49,6,0,49,6,49:0,0,1,1
chrM    74      .       T       <NON_REF>       .       .       END=84                                                               GT:DP:GQ:MIN_DP:PL      0/0:2:6:2:0,6,40
chrM    85      .       G       <NON_REF>       .       .       END=133                                                              GT:DP:GQ:MIN_DP:PL      0/0:3:9:3:0,9,56
chrM    134     .       T       <NON_REF>       .       .       END=149                                                              GT:DP:GQ:MIN_DP:PL      0/0:3:6:2:0,6,56
```

```bash
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
```

## Access to other databases:

```bash
# GRCh37
wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh37/vcf/GRCh37.variant_call.all.vcf.gz

#GRCh38
wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh38/vcf/GRCh38.variant_call.all.vcf.gz
```
