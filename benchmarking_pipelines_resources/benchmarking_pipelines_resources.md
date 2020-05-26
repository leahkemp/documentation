# Benchmarking genomic pipelines - resources

Created: 2020-04-22 13:37:04
Last modified: 2020/05/26 11:55:17

- **Aim:** Undertake benchmarking of genomics pipelines to optimise their resource use.
- **Prerequisite software:** [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html), [samtools 1.9](http://www.htslib.org/), [bedtools 2.25](https://bedtools.readthedocs.io/en/latest/)
- **OS:** Ubuntu 16.04 (Wintermute - research server) (*add info about production OS*)

The idea is to run these pipelines ([human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) and [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline)) against a reduced Genome In A Bottle (GIAB) sample [NIST7035](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) dataset and evaluate the real time and user time for each step in the pipelines. Once the resources are optimised for running these pipelines on Wintermute (research server), resource benchmarking will be run again on production to fine tune the optimal resource use for these pipelines running on this machine.

## Table of contents

- [Benchmarking genomic pipelines - resources](#benchmarking-genomic-pipelines---resources)
  - [Table of contents](#table-of-contents)
  - [Create a test dataset](#create-a-test-dataset)
    - [Download and prepare the test WES data (NIST7035)](#download-and-prepare-the-test-wes-data-nist7035)
    - [Reduce the test dataset size for benchmarking](#reduce-the-test-dataset-size-for-benchmarking)
      - [Map NIST7035 to the reference genome (ORAC)](#map-nist7035-to-the-reference-genome-orac)
      - [Extract chr1 from the mapped bam files (retain paired reads) (wintermute)](#extract-chr1-from-the-mapped-bam-files-retain-paired-reads-wintermute)
      - [Extract the paired fastq reads (R1 and R2) from the raw fastq reads that map to chr1 (wintermute)](#extract-the-paired-fastq-reads-r1-and-r2-from-the-raw-fastq-reads-that-map-to-chr1-wintermute)
  - [Testing](#testing)
    - [Setup](#setup)
    - [Wintermute](#wintermute)
      - [human_genomics_pipeline](#humangenomicspipeline)
        - [Single sample](#single-sample)
        - [Cohort sample](#cohort-sample)
      - [vcf_annotation_pipeline](#vcfannotationpipeline)
        - [Single sample](#single-sample-1)
        - [Cohort sample](#cohort-sample-1)
        - [Both single and cohort samples](#both-single-and-cohort-samples)
    - [Production](#production)
      - [human_genomics_pipeline](#humangenomicspipeline-1)
      - [vcf_annotation_pipeline](#vcfannotationpipeline-1)

## Create a test dataset

### Download and prepare the test WES data (NIST7035)

Download

```bash
cd /store/lkemp/publicData/exomes/NA12878_exome/
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/Garvan_NA12878_HG001_HiSeq_Exome.README
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L002_R1_001.fastq.gz
wget https://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L002_R2_001.fastq.gz
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
cat NIST7035*_R1_001.fastq.gz > NIST7035_NIST_R1.fastq.gz
cat NIST7035*_R2_001.fastq.gz > NIST7035_NIST_R2.fastq.gz
```

### Reduce the test dataset size for benchmarking

I will extract the fastq reads that map to chr1 so as to reduce the size of the dataset (width) without reducing the depth of reads that could influence the performance of the downstream rules (which could skew the resource benchmarking results for these steps)

#### Map NIST7035 to the reference genome (ORAC)

To increase speed, I will use the bwa-mem step of parabricks on orac

```bash
pbrun fq2bam \
--ref /NGS/scratch/KSCBIOM/HumanGenomics/publicData/human_refs/GRCh37/ucsc.hg19.fasta \
--in-fq /NGS/scratch/KSCBIOM/HumanGenomics/publicData/human_refs/GIAB/NA12878_exome/NIST7035_NIST_R1.fastq.gz /NGS/scratch/KSCBIOM/HumanGenomics/publicData/human_refs/GIAB/NA12878_exome/NIST7035_NIST_R2.fastq.gz \
--out-bam /home/lkemp/resource_benchmarking/NIST7035_NIST.bam
```

#### Extract chr1 from the mapped bam files (retain paired reads) (wintermute)

See [this documentation](https://gist.github.com/darencard/72ddd9e6c08aaff5ff64ca512a04a6dd))

```bash
cd /home/lkemp/resource_benchmarking/

samtools view -u -f 1 -F 12 NIST7035_NIST.bam chr1 -b > NIST7035_NIST_chr1.bam
```

#### Extract the paired fastq reads (R1 and R2) from the raw fastq reads that map to chr1 (wintermute)

Extract:

- reads that mapped properly as pairs
- reads that didn’t map properly as pairs (both didn’t map, or one didn’t map)

```bash
# R1 unmapped, R2 mapped
samtools view -u -f 4 -F 264 NIST7035_NIST_chr1.bam > NIST7035_NIST_chr1_unmap_map.bam
# R1 mapped, R2 unmapped
samtools view -u -f 8 -F 260 NIST7035_NIST_chr1.bam > NIST7035_NIST_chr1_map_unmap.bam
# R1 & R2 unmapped
samtools view -u -f 12 -F 256 NIST7035_NIST_chr1.bam > NIST7035_NIST_chr1_unmap_unmap.bam
```

Merge the three files that contain at least one unmapped pair

```bash
samtools merge -u NIST7035_NIST_chr1_unmapped.bam NIST7035_NIST_chr1_unmap_map.bam NIST7035_NIST_chr1_map_unmap.bam NIST7035_NIST_chr1_unmap_unmap.bam
```

Sort

```bash
samtools sort -n NIST7035_NIST_chr1_map_map.bam NIST7035_NIST_chr1_mapped.sort
samtools sort -n NIST7035_NIST_chr1_unmapped.bam NIST7035_NIST_chr1_unmapped.sort
```

Check the files

```bash
samtools flagstat NIST7035_NIST_chr1.sorted.md.bam
```

Extract the FASTQ reads into two paired read files

```bash
bamToFastq -i NIST7035_NIST_chr1_mapped.sort.bam -fq NIST7035_NIST_chr1_mapped.1.fastq -fq2 NIST7035_NIST_chr1_mapped.2.fastq
bamToFastq -i NIST7035_NIST_chr1_unmapped.sort.bam -fq NIST7035_NIST_chr1_unmapped.1.fastq -fq2 NIST7035_NIST_chr1_unmapped.2.fastq
```

Combine both the first and paired reads together from the mapped and unmapped files

```bash
cat NIST7035_NIST_chr1_mapped.1.fastq NIST7035_NIST_chr1_unmapped.1.fastq > NIST7035_NIST_chr1.1.fastq
cat NIST7035_NIST_chr1_mapped.2.fastq NIST7035_NIST_chr1_unmapped.2.fastq > NIST7035_NIST_chr1.2.fastq
```

## Testing

Approach:

- Run each rule separately
- Run with doubling threads: 1, 2, 4, 8, 16 etc. (until runtime plateaus)
- Compare real time vs. user time (minimise any divergence between them)
- Re-run benchmarking on new machines (eg. production) to fine-tune resource allocation

### Setup

Replace the full fastq input file with the reduced ones

```bash
cd /store/lkemp/exome_project/resource_benchmarking/fastq/
rm -r NIST70*
cp ../human_genomics_pipeline/mapped/NIST7035_NIST_chr1* .
cp ../human_genomics_pipeline/mapped/NIST7086_NIST_chr1* .
```

Clone a fresh pipeline and checkout the branch for resource benchmarking on Wintermute

```bash
cd ..
rm -r human_genomic_pipeline
git clone git@github.com:ESR-NZ/human_genomics_pipeline.git
cd human_genomics_pipeline
git checkout resource_benchmarking_wintermute
```

### Wintermute

#### human_genomics_pipeline

Set the number of threads in each rule to the maximum number that we will test (32), therefore we can control the number of threads used for each test with the `j` flag passed to snakemake on the command line (if the number of threads for a given rule are larger that the threads passed to this flag, they will be scaled down)

Run rules in order shown

*note. I removed all files created by a rule before re-running the same rule with a new number of threads. I also made sure a rules conda environment was created prior to running benhcmarking for a rule to ensure the time represents on;y activatea conda environment and not creating one from scratch. One sample was used NIST7035_NIST_chr1*

- fastqc

```bash
time snakemake -j 1 --use-conda --configfile config/resource_benchmarking.yml --until fastqc
time snakemake -j 2 --use-conda --configfile config/resource_benchmarking.yml --until fastqc
time snakemake -j 4 --use-conda --configfile config/resource_benchmarking.yml --until fastqc
time snakemake -j 8 --use-conda --configfile config/resource_benchmarking.yml --until fastqc
time snakemake -j 16 --use-conda --configfile config/resource_benchmarking.yml --until fastqc
time snakemake -j 32 --use-conda --configfile config/resource_benchmarking.yml --until fastqc
```

- multiqc

```bash
time snakemake -j 1 --use-conda --configfile config/resource_benchmarking.yml --until multiqc_pre_trim
time snakemake -j 2 --use-conda --configfile config/resource_benchmarking.yml --until multiqc_pre_trim
time snakemake -j 4 --use-conda --configfile config/resource_benchmarking.yml --until multiqc_pre_trim
time snakemake -j 8 --use-conda --configfile config/resource_benchmarking.yml --until multiqc_pre_trim
time snakemake -j 16 --use-conda --configfile config/resource_benchmarking.yml --until multiqc_pre_trim
time snakemake -j 32 --use-conda --configfile config/resource_benchmarking.yml --until multiqc_pre_trim
```

- trim_galore

```bash
time snakemake -j 1 --use-conda --configfile config/resource_benchmarking.yml --until trim_galore_pe
time snakemake -j 2 --use-conda --configfile config/resource_benchmarking.yml --until trim_galore_pe
time snakemake -j 4 --use-conda --configfile config/resource_benchmarking.yml --until trim_galore_pe
time snakemake -j 8 --use-conda --configfile config/resource_benchmarking.yml --until trim_galore_pe
time snakemake -j 16 --use-conda --configfile config/resource_benchmarking.yml --until trim_galore_pe
time snakemake -j 32 --use-conda --configfile config/resource_benchmarking.yml --until trim_galore_pe
```

- multiqc_post_trim

```bash
time snakemake -j 1 --use-conda --configfile config/resource_benchmarking.yml --until multiqc_post_trim
time snakemake -j 2 --use-conda --configfile config/resource_benchmarking.yml --until multiqc_post_trim
time snakemake -j 4 --use-conda --configfile config/resource_benchmarking.yml --until multiqc_post_trim
time snakemake -j 8 --use-conda --configfile config/resource_benchmarking.yml --until multiqc_post_trim
time snakemake -j 16 --use-conda --configfile config/resource_benchmarking.yml --until multiqc_post_trim
time snakemake -j 32 --use-conda --configfile config/resource_benchmarking.yml --until multiqc_post_trim
```

- bwa_map

```bash
time snakemake -j 1 --use-conda --configfile config/resource_benchmarking.yml --until bwa_map
time snakemake -j 2 --use-conda --configfile config/resource_benchmarking.yml --until bwa_map
time snakemake -j 4 --use-conda --configfile config/resource_benchmarking.yml --until bwa_map
time snakemake -j 8 --use-conda --configfile config/resource_benchmarking.yml --until bwa_map
time snakemake -j 16 --use-conda --configfile config/resource_benchmarking.yml --until bwa_map
time snakemake -j 32 --use-conda --configfile config/resource_benchmarking.yml --until bwa_map
```

- sambamba_sort

```bash
time snakemake -j 1 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_sort
time snakemake -j 2 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_sort
time snakemake -j 4 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_sort
time snakemake -j 8 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_sort
time snakemake -j 16 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_sort
time snakemake -j 32 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_sort
```

- sambamba_mkdups

```bash
time snakemake -j 1 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_mkdups
time snakemake -j 2 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_mkdups
time snakemake -j 4 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_mkdups
time snakemake -j 8 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_mkdups
time snakemake -j 16 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_mkdups
time snakemake -j 32 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_mkdups
```

- sambamba_index

```bash
time snakemake -j 1 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_index
time snakemake -j 2 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_index
time snakemake -j 4 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_index
time snakemake -j 8 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_index
time snakemake -j 16 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_index
time snakemake -j 32 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_index
```

- gatk_add_replace_read_groups

```bash
time snakemake -j 1 --use-conda --configfile config/resource_benchmarking.yml --until gatk_add_replace_read_groups
time snakemake -j 2 --use-conda --configfile config/resource_benchmarking.yml --until gatk_add_replace_read_groups
time snakemake -j 4 --use-conda --configfile config/resource_benchmarking.yml --until gatk_add_replace_read_groups
time snakemake -j 8 --use-conda --configfile config/resource_benchmarking.yml --until gatk_add_replace_read_groups
time snakemake -j 16 --use-conda --configfile config/resource_benchmarking.yml --until gatk_add_replace_read_groups
time snakemake -j 32 --use-conda --configfile config/resource_benchmarking.yml --until gatk_add_replace_read_groups
```

- sambamba_index_rgadd

```bash
time snakemake -j 1 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_index_rgadd
time snakemake -j 2 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_index_rgadd
time snakemake -j 4 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_index_rgadd
time snakemake -j 8 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_index_rgadd
time snakemake -j 16 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_index_rgadd
time snakemake -j 32 --use-conda --configfile config/resource_benchmarking.yml --until sambamba_index_rgadd
```

- gatk_base_recalibrator

```bash
time snakemake -j 1 --use-conda --configfile config/resource_benchmarking.yml --until gatk_base_recalibrator
time snakemake -j 2 --use-conda --configfile config/resource_benchmarking.yml --until gatk_base_recalibrator
time snakemake -j 4 --use-conda --configfile config/resource_benchmarking.yml --until gatk_base_recalibrator
time snakemake -j 8 --use-conda --configfile config/resource_benchmarking.yml --until gatk_base_recalibrator
time snakemake -j 16 --use-conda --configfile config/resource_benchmarking.yml --until gatk_base_recalibrator
time snakemake -j 32 --use-conda --configfile config/resource_benchmarking.yml --until gatk_base_recalibrator
```

- gatk_apply_bqsr

```bash
time snakemake -j 1 --use-conda --configfile config/resource_benchmarking.yml --until gatk_apply_bqsr
time snakemake -j 2 --use-conda --configfile config/resource_benchmarking.yml --until gatk_apply_bqsr
time snakemake -j 4 --use-conda --configfile config/resource_benchmarking.yml --until gatk_apply_bqsr
time snakemake -j 8 --use-conda --configfile config/resource_benchmarking.yml --until gatk_apply_bqsr
time snakemake -j 16 --use-conda --configfile config/resource_benchmarking.yml --until gatk_apply_bqsr
time snakemake -j 32 --use-conda --configfile config/resource_benchmarking.yml --until gatk_apply_bqsr
```

##### Single sample

- gatk_haplotype_caller_single

```bash
time snakemake -j 1 --use-conda --configfile config/resource_benchmarking.yml --until gatk_haplotype_caller_single
time snakemake -j 2 --use-conda --configfile config/resource_benchmarking.yml --until gatk_haplotype_caller_single
time snakemake -j 4 --use-conda --configfile config/resource_benchmarking.yml --until gatk_haplotype_caller_single
time snakemake -j 8 --use-conda --configfile config/resource_benchmarking.yml --until gatk_haplotype_caller_single
time snakemake -j 16 --use-conda --configfile config/resource_benchmarking.yml --until gatk_haplotype_caller_single
time snakemake -j 32 --use-conda --configfile config/resource_benchmarking.yml --until gatk_haplotype_caller_single
```

##### Cohort sample

- gatk_haplotype_caller_gvcf

- gatk_combine_gvcf

- gatk_genotype_gvcf

#### vcf_annotation_pipeline

resource_benchmarking_wintermute branch?

##### Single sample

- gatk_cnn_score_variants

- gatk_filter_variant_tranches

##### Cohort sample

- gatk_variant_recalibrator_indel

- gatk_variant_recalibrator_snp

- gatk_vqsr_indel

- gatk_vqsr_snp

##### Both single and cohort samples

- snpsift_dbnsfp

- vep

- genmod_cadd

### Production

#### human_genomics_pipeline

#### vcf_annotation_pipeline