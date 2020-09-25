# Benchmarking genomic pipelines - resources - wintermute - per rule runtime

Created: 2020-09-11 13:37:04
Last modified: 2020/09/25 16:42:07

- **Aim:** Undertake benchmarking of genomics pipelines to optimise the threading of each rule in the pipelines.
- **Prerequisite software:** [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html), [wget](https://www.gnu.org/software/wget/)
- **OS:** Ubuntu 16.04 (Wintermute - research server)

The idea is to run these pipelines ([human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) and [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline)) against the Genome In A Bottle (GIAB) sample [NIST7035](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) exome and evaluate the clock time for each step in the pipelines with increased threading to evaluate the point of diminishing return. This will also provide overall pipeline run times. These tests will be undertaken on Wintermute.

## Table of contents

- [Benchmarking genomic pipelines - resources - wintermute - per rule runtime](#benchmarking-genomic-pipelines---resources---wintermute---per-rule-runtime)
  - [Table of contents](#table-of-contents)
  - [Create a test dataset](#create-a-test-dataset)
    - [Setup](#setup)
    - [Benchmarking](#benchmarking)
    - [Results](#results)
  - [Previous resource benchmarking](#previous-resource-benchmarking)

## Create a test dataset

Download

```bash
cd /store/lkemp/publicData/NA12878_exome/
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
cat NIST7035*_R1_001.fastq.gz > NIST7035_NIST_1.fastq.gz
cat NIST7035*_R2_001.fastq.gz > NIST7035_NIST_2.fastq.gz
```

Create dummy vcf files for cohort runs

```bash
cp NIST7035_NIST_1.fastq.gz dummy1_NIST_1.fastq.gz
cp NIST7035_NIST_2.fastq.gz dummy1_NIST_2.fastq.gz
cp NIST7035_NIST_1.fastq.gz dummy2_NIST_1.fastq.gz
cp NIST7035_NIST_2.fastq.gz dummy2_NIST_2.fastq.gz
```
Create a dummy pedigree file for cohort runs

```txt
NIST7035_NIST	dummy1_NIST	0	0	1	1
NIST7035_NIST	dummy2_NIST	0	0	2	1
NIST7035_NIST	NIST7035_NIST	dummy1_NIST	dummy2_NIST	1	2

```

### Setup

Setup folders

```bash
cd store/lkemp/
mkdir resource_benchmarking
cd resource_benchmarking
mkdir resource_bench_01_thread
mkdir resource_bench_02_thread
mkdir resource_bench_04_thread
mkdir resource_bench_08_thread
mkdir resource_bench_16_thread
mkdir resource_bench_32_thread
cd resource_bench_01_thread
mkdir single # repeat for all directories
mkdir cohort # repeat for all directories
```

Copy the reduced dataset to where we will do the resource benchmarking

```bash
cd /store/lkemp/resource_benchmarking/resource_bench_01_thread/single/
mkdir fastq
cd fastq/ # repeat for all directories
cp /store/lkemp/publicData/NA12878_exome/NIST7035_NIST_1.fastq.gz . # repeat for all directories
cp /store/lkemp/publicData/NA12878_exome/NIST7035_NIST_2.fastq.gz . # repeat for all directories
cp /store/lkemp/publicData/NA12878_exome/dummy1_NIST_1.fastq.gz . # for cohort runs, repeat for all directories
cp /store/lkemp/publicData/NA12878_exome/dummy1_NIST_2.fastq.gz . # for cohort runs, repeat for all directories
cp /store/lkemp/publicData/NA12878_exome/dummy2_NIST_1.fastq.gz . # for cohort runs, repeat for all directories
cp /store/lkemp/publicData/NA12878_exome/dummy2_NIST_2.fastq.gz . # for cohort runs, repeat for all directories
```

Clone forked pipeline and create/checkout the branch for resource benchmarking

```bash
cd /store/lkemp/resource_benchmarking/resource_bench_01_thread/single/
git clone git@github.com:leahkemp/human_genomics_pipeline.git
cd human_genomics_pipeline
git branch resource_benchmarking
```

```bash
cd /store/lkemp/exome_project/resource_benchmarking/
git clone git@github.com:ESR-NZ/human_genomics_pipeline.git
cd human_genomics_pipeline
git checkout resource_benchmarking
```

### Benchmarking

Workflow:

- Create a resource_benchmarking branch for the pipeline

- Set the number of threads *within each rule* to the maximum number that we will test (32), therefore we can control the number of threads used for each test with the `-j` flag passed to Snakemake on the command line (if the number of threads for a given rule are larger that the threads passed to this flag, they will be scaled down)

- Wrap each benchmarking file in each rule with `repeat("benchmarks/rule/{sample}.tsv, 3)` to get the rule to run 3 times to get a measure of the variability of measurements

- Double the number of threads used each run (in run.sh script), as set by the `-j` parameter: 1, 2, 4, 8, 16, 32

- Repeat for runs against single samples and cohort samples to benchmark all the rules in the pipelines

- Collate the benchmarking output with a python script (merge_resource_benchmarking.py)

- Plot the times from the benchmarking files

- Repeat for vcf_annotation_pipeline
  
### Results

Full run settings and outputs for each threading level for can be found at:

- [resource_bench_1_thread](https://github.com/leahkemp/human_genomics_pipeline/tree/resource_bench_1_thread)
- [resource_bench_2_thread](https://github.com/leahkemp/human_genomics_pipeline/tree/resource_bench_2_threads)
- [resource_bench_4_thread](https://github.com/leahkemp/human_genomics_pipeline/tree/resource_bench_4_threads)
- [resource_bench_8_thread](https://github.com/leahkemp/human_genomics_pipeline/tree/resource_bench_8_threads)
- [resource_bench_16_thread](https://github.com/leahkemp/human_genomics_pipeline/tree/resource_bench_16_threads)
- [resource_bench_32_thread](https://github.com/leahkemp/human_genomics_pipeline/tree/resource_bench_32_threads)

- [resource_bench_1_thread_cohort](https://github.com/leahkemp/human_genomics_pipeline/tree/resource_bench_1_thread_cohort)
- [resource_bench_2_thread_cohort](https://github.com/leahkemp/human_genomics_pipeline/tree/resource_bench_2_threads_cohort)
- [resource_bench_4_thread_cohort](https://github.com/leahkemp/human_genomics_pipeline/tree/resource_bench_4_threads_cohort)
- [resource_bench_8_thread_cohort](https://github.com/leahkemp/human_genomics_pipeline/tree/resource_bench_8_threads_cohort)
- [resource_bench_16_thread_cohort](https://github.com/leahkemp/human_genomics_pipeline/tree/resource_bench_16_threads_cohort)
- [resource_bench_32_thread_cohort](https://github.com/leahkemp/human_genomics_pipeline/tree/resource_bench_32_threads_cohort)
  
Plotting/report can be found [here](results/)

## Previous resource benchmarking

Previous resource benchmarking on an older version of the pipeline

Run settings and results: 

- [resource_bench2.0](https://github.com/ESR-NZ/human_genomics_pipeline/tree/resource_bench2.0)
- [resource_bench2.1](https://github.com/ESR-NZ/human_genomics_pipeline/tree/resource_bench2.1)
- [resource_bench2.2](https://github.com/ESR-NZ/human_genomics_pipeline/tree/resource_bench2.2)
- [resource_bench2.3](https://github.com/ESR-NZ/human_genomics_pipeline/tree/resource_bench2.3)
- [resource_bench2.4](https://github.com/ESR-NZ/human_genomics_pipeline/tree/resource_bench2.4)
- [resource_bench2.5](https://github.com/ESR-NZ/human_genomics_pipeline/tree/resource_bench2.5)

Extended data, methods and results: https://github.com/leahkemp/documentation/tree/resource_benchmarking/benchmarking_pipelines_resources