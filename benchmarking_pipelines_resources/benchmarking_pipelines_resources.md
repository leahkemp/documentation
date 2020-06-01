# Benchmarking genomic pipelines - resources

Created: 2020-04-22 13:37:04
Last modified: 2020/05/28 18:51:21

- **Aim:** Undertake benchmarking of genomics pipelines to optimise their resource use.
- **Prerequisite software:** [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html), [samtools 1.9](http://www.htslib.org/), [bedtools 2.25](https://bedtools.readthedocs.io/en/latest/), [bgzip 1.2.1](http://www.htslib.org/doc/bgzip.html)
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
    - [Production](#production)
      - [human_genomics_pipeline](#human_genomics_pipeline)
      - [vcf_annotation_pipeline](#vcf_annotation_pipeline)

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

Extract chromosome one

```bash
cd /store/lkemp/exome_project/resource_benchmarking/
samtools view NIST7035_NIST.bam chr1 -b > NIST7035_NIST_chr1.bam
```

Extract the reads for which both paired reads were mapped

See [this documentation](https://gist.github.com/darencard/72ddd9e6c08aaff5ff64ca512a04a6dd))

```bash
samtools view -u -f 1 -F 12 NIST7035_NIST_chr1.bam > NIST7035_NIST_chr1_map_map.bam
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
samtools sort -n NIST7035_NIST_chr1_map_map.bam -o NIST7035_NIST_chr1_mapped.sort
samtools sort -n NIST7035_NIST_chr1_unmapped.bam -o NIST7035_NIST_chr1_unmapped.sort
```

Check that the number of unmapped and mapped reads total the number of reads in the original bam file (for chr1)

```bash
# Original file
samtools flagstat NIST7035_NIST_chr1.bam
```

My output:

```bash
8223836 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
2650 + 0 supplementary
933687 + 0 duplicates
8214293 + 0 mapped (99.88% : N/A)
8221186 + 0 paired in sequencing
4111231 + 0 read1
4109955 + 0 read2
8152956 + 0 properly paired (99.17% : N/A)
8202100 + 0 with itself and mate mapped
9543 + 0 singletons (0.12% : N/A)
8456 + 0 with mate mapped to a different chr
5219 + 0 with mate mapped to a different chr (mapQ>=5)
```

```bash
# Mapped
samtools view -c NIST7035_NIST_chr1_mapped.sort
# Unmapped
samtools view -c NIST7035_NIST_chr1_unmapped.sort
```

My output:

```bash
8204729
19107
```

Extract the mapped FASTQ reads into two paired read files

```bash
bamToFastq -i NIST7035_NIST_chr1_mapped.sort -fq NIST7035_NIST_chr1_mapped.1.fastq -fq2 NIST7035_NIST_chr1_mapped.2.fastq
bamToFastq -i NIST7035_NIST_chr1_unmapped.sort -fq NIST7035_NIST_chr1_unmapped.1.fastq -fq2 NIST7035_NIST_chr1_unmapped.2.fastq
```

Combine both the first and paired reads together from the mapped and unmapped files

```bash
cat NIST7035_NIST_chr1_mapped.1.fastq NIST7035_NIST_chr1_unmapped.1.fastq > NIST7035_NIST_chr1_R1.fastq
cat NIST7035_NIST_chr1_mapped.2.fastq NIST7035_NIST_chr1_unmapped.2.fastq > NIST7035_NIST_chr1_R2.fastq
```

Bgzip

```bash
bgzip NIST7035_NIST_chr1_R1.fastq
bgzip NIST7035_NIST_chr1_R2.fastq
```

## Testing

Approach:

- Run each rule separately
- Run with doubling threads: 1, 2, 4, 8, 16 etc. (until runtime plateaus)
- Compare real time vs. user time (minimise any divergence between them)
- Re-run benchmarking on new machines (eg. production) to fine-tune resource allocation

### Setup

Copy the reduced dataset (fastq files for chr1) to where we will do the resource benchmarking

```bash
mkdir fastq
cd fastq/

cp ../NIST7035_NIST_chr1_R1.fastq.gz .
cp ../NIST7035_NIST_chr1_R2.fastq.gz .
```

Clone pipeline and checkout the branch for resource benchmarking on Wintermute

```bash
cd /store/lkemp/exome_project/resource_benchmarking/
git clone git@github.com:ESR-NZ/human_genomics_pipeline.git
cd human_genomics_pipeline
git checkout resource_benchmarking_wintermute
```

### Wintermute

- Create a resource_benchmarking branch for the pipeline

- Set the number of threads in each rule to the maximum number that we will test (32), therefore we can control the number of threads used for each test with the `j` flag passed to snakemake on the command line (if the number of threads for a given rule are larger that the threads passed to this flag, they will be scaled down)

- Create a 'times' dir

- Wrap each rule script with ( time rule_script 2> rule.stderr ) 2> times/rule_time.txt

- This will write the output of the time command to a file (within times/) for each rule (and extract the output messages to the .stderr files)

- Create shell scripts that will prompt the pipeline to run

- For each run (shell script), double the number of threads used each run, as set by the `-j` parameter: 1, 2, 4, 8, 16, 32

- Repeat for runs against single samples and cohort samples to benchmark all the rules in the pipelines

Run through each pipeline run

```bash
# human_genomics_pipeline
bash run_1_threads.sh
bash run_2_threads.sh
bash run_4_threads.sh
bash run_8_threads.sh
bash run_16_threads.sh
bash run_32_threads.sh

# vcf_annotation_pipeline
bash run_1_thread.sh
bash run_2_thread.sh
bash run_4_thread.sh
bash run_8_thread.sh
bash run_16_thread.sh
bash run_32_thread.sh
```

Extract the times from the output times files and plot

### Production

#### human_genomics_pipeline

#### vcf_annotation_pipeline