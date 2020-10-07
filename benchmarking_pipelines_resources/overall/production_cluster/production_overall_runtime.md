# Benchmarking genomic pipelines - resources - production - overall runtime

Created: 2020/09/25 16:29:41
Last modified: 2020/09/25 18:59:07

- **Aim:** Undertake benchmarking of genomics pipelines to optimise their overall runtimes and establish the number of samples at which running the pipeline **not** gpu accelerated is faster in overall runtime. These tests will include the optimal per-rule threading that was established in the related tests, [benchmarking genomic pipelines - resources - production - per rule runtime](../benchmarking_pipelines_resources/benchmarking_pipelines_resources_production_per_rule_runtime.md) and [benchmarking genomic pipelines - resources - wintermute - per rule runtime](../benchmarking_pipelines_resources/benchmarking_pipeline_resources_wintermute_per_rule_runtime.md).
- **Prerequisite software:** [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html), [wget](https://www.gnu.org/software/wget/)
- **OS:** ESR production cluster

## Table of contents

- [Benchmarking genomic pipelines - resources - production - overall runtime](#benchmarking-genomic-pipelines---resources---production---overall-runtime)
  - [Table of contents](#table-of-contents)
    - [Setup](#setup)

### Setup

Variables tested:

- cpu/gpu runs
- 1, 8, 16, 32 samples
- single runs
- 32, 64, 96 max threads per pipeline run

*Note. one 'sample' for a cohort run will mean one family (and one sample is comprised of three individuals/exomes)

Create folder structure

```bash
cd /NGS/scratch/KSCBIOM/HumanGenomics/resource_benchmarking/resource_benchmarking_total_runtime/

# Create initial folder structure
mkdir -p exome/{cpu_run/{01_sample/single_run/,08_sample/single_run,16_sample/single_run,32_sample/single_run},gpu_run/{01_sample/single_run/,08_sample/single_run,16_sample/single_run,32_sample/single_run}}

# Create threading directories for cpu_runs dirs
for i in exome/*/*/*; do
  mkdir $i/32_threads
  mkdir $i/64_threads
  mkdir $i/96_threads
  done

# Create pipeline input folders
for i in exome/*/*/single_run/*; do
  mkdir $i/fastq
  mkdir $i/bams
  mkdir $i/vcf
  done
```