# Benchmarking pipeline resources - per rule - production cluster

Created: 2020-09-11 13:37:04
Last modified: 2020/10/07 14:51:05

- **Aim:** Undertake benchmarking of genomics pipelines to optimise the threading of each rule in the pipelines.
- **Prerequisite software:** [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html), [wget](https://www.gnu.org/software/wget/)
- **OS:** ESR production cluster

The idea is to run these pipelines ([human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) and [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline)) and evaluate the clock time for each step in the pipelines with increased threading to evaluate the point of diminishing return. These tests will extend those undertaken on Wintermute (see [these docs](../research_server/resource_benchmarking_per_rule_research_server.md)) by accounting for the number of samples analysed in a pipeline run at one time. These tests will be undertaken on the ESR production cluster in order to scale up the sample number in the pipeline runs as well as get more accurate pipeline run times since we will be analysing our clinical exomes on the ESR production cluster. These tests will be undertaken on exomes sequenced at ESR (as opposed to the public exome used in [these tests](../research_server/resource_benchmarking_per_rule_research_server.md)) because there are multiple exomes to use for testing pipeline runs with multiple samples and the size of the input fastq files will be more comparable to the size of the clinical exomes we will be analysing.

## Table of contents

- [Benchmarking pipeline resources - per rule - production cluster](#benchmarking-pipeline-resources---per-rule---production-cluster)
  - [Table of contents](#table-of-contents)
    - [Setup](#setup)
    - [Results](#results)
  - [Previous resource benchmarking](#previous-resource-benchmarking)

### Setup

Variables tested:

- cpu runs
- 1/8 samples/families
- single/cohort runs (to get runtimes for all rules)
- 1-32 threads per rule

*Note. one 'sample' for a cohort run will mean one family (and one sample is comprised of three individuals/exomes)

Create folder structure

```bash
cd /NGS/scratch/KSCBIOM/HumanGenomics/resource_benchmarking/

# Create initial folder structure
mkdir -p exome/cpu_run/{01_sample/{single_run,cohort_run},08_sample/{single_run,cohort_run}}

# Create threading directories for cpu_runs dirs
for i in exome/*/*/*; do
  mkdir $i/01_threads
  mkdir $i/02_threads
  mkdir $i/04_threads
  mkdir $i/08_threads
  mkdir $i/16_threads
  mkdir $i/32_threads
  done

# Create pipeline input folders for single runs
for i in exome/*/*/single_run/*; do
  mkdir $i/fastq
  mkdir $i/bams
  mkdir $i/vcf
  done

# Create pipeline input folders for cohort runs
for i in exome/*/*/cohort_run/*; do
  mkdir $i/fastq
  mkdir $i/pedigrees
  mkdir $i/vcf
  done
```

Create dummy vcf files for cohort runs and runs with multiple samples (first set the directory to the downloaded data in `make_dummy_data.sh`)

```bash
bash make_dummy_data.sh
```

Create dummy pedigree files for cohort runs (families 1-8) in the same location as the dummy data. For example:

```txt
dummy_proband01	dummy_father01	0	0	1	1
dummy_proband01	dummy_mother01	0	0	2	1
dummy_proband01	dummy_proband01	dummy_father01	dummy_mother01	1	2

```

Copy the fastq data and pedigree files to where we will do the resource benchmarking (first set the directory to the downloaded data in `populate_fastq_dirs.sh`)

```bash
bash populate_fastq_pedigree_dirs.sh
```

Clone forked pipeline and create/checkout the branch for resource benchmarking which includes modifications for resource benchmarking (each rule wrapped with `repeat("benchmarks/rule/{sample}.tsv, 3)`)

```bash
# human_genomics_pipeline
for filedir in /NGS/scratch/KSCBIOM/HumanGenomics/resource_benchmarking/exome/*/*/*/*; do
  cd $filedir
  git clone https://github.com/leahkemp/human_genomics_pipeline.git
  cd human_genomics_pipeline
  git checkout bc2f1ae711a58267708b2e0f6c5b9c9453ef3b77
  done

# vcf_annotation_pipeline
for filedir in /NGS/scratch/KSCBIOM/HumanGenomics/resource_benchmarking/exome/*/*/*/*; do
  cd $filedir
  git clone https://github.com/leahkemp/vcf_annotation_pipeline.git
  cd vcf_annotation_pipeline
  git checkout 2747b87237a08a3ab70765bca94a9a6c7386b1db
  done
```

Setup directories that contain all the variations of rule files, config files, run scripts and cluster configs that we will need. And add the appropriate modified files to these directories

```bash
cd /NGS/scratch/KSCBIOM/HumanGenomics/resource_benchmarking

# Set the threads in all rule files from 1-32 - for human_genomics_pipeline (hgp) and vcf_annotation_pipeline (vap)
mkdir rules_01_thread_hgp
mkdir rules_02_thread_hgp
mkdir rules_04_thread_hgp
mkdir rules_08_thread_hgp
mkdir rules_16_thread_hgp
mkdir rules_32_thread_hgp

mkdir rules_01_thread_vap
mkdir rules_02_thread_vap
mkdir rules_04_thread_vap
mkdir rules_08_thread_vap
mkdir rules_16_thread_vap
mkdir rules_32_thread_vap

# Set the 'DATA', 'GPU_ACCELERATED' and 'RECALIBRATION RESOURCES' parameters - for human_genomics_pipeline (hgp) and vcf_annotation_pipeline (vap)
mkdir config_cohort_cpu_hgp
mkdir config_single_cpu_hgp
mkdir config_cohort_gpu_hgp
mkdir config_single_gpu_hgp

mkdir config_cohort_cpu_vap
mkdir config_single_cpu_vap
mkdir config_cohort_gpu_vap
mkdir config_single_gpu_vap

# Set the overall threading to 32 (-j 32) and set any other parameters - for human_genomics_pipeline (hgp) and vcf_annotation_pipeline (vap)
mkdir runscript_hgp

mkdir runscript_vap

mkdir cluster_config
```

Copy the config files, run scripts and cluster configs to where we will do the resource benchmarking (first set the working directory in `populate_configs_hgp.sh`)

```bash
bash populate_configs.sh
```

Run all human_genomics_pipeline runs. Before running vcf_annotation_pipeline runs, move the outputs of human_genomics_pipeline to the directories needed by vcf_annotation_pipeline

```bash
bash move_input_files_for_vap.sh
```

Run all vcf_annotation_pipeline runs

Once the all pipeline runs have been carried out, merge the outputs of all benchmarking output files into one csv file

```bash
python merge_resource_benchmarking_production.py
```

Get the node information where available from the log files

```bash
python get_nodes.py
```

Merge the two csv files

```bash
python merge_csvs.py
```

### Results

## Previous resource benchmarking

Previous resource benchmarking on an older version of the pipeline

Full run settings and outputs for each threading level for can be found at:

- [resource_bench2.0](https://github.com/ESR-NZ/human_genomics_pipeline/tree/resource_bench2.0)
- [resource_bench2.1](https://github.com/ESR-NZ/human_genomics_pipeline/tree/resource_bench2.1)
- [resource_bench2.2](https://github.com/ESR-NZ/human_genomics_pipeline/tree/resource_bench2.2)
- [resource_bench2.3](https://github.com/ESR-NZ/human_genomics_pipeline/tree/resource_bench2.3)
- [resource_bench2.4](https://github.com/ESR-NZ/human_genomics_pipeline/tree/resource_bench2.4)
- [resource_bench2.5](https://github.com/ESR-NZ/human_genomics_pipeline/tree/resource_bench2.5)

Extended data, methods and results: https://github.com/leahkemp/documentation/tree/resource_benchmarking/benchmarking_pipelines_resources
