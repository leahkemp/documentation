#!/bin/bash

workingdir="/NGS/scratch/KSCBIOM/HumanGenomics/resource_benchmarking/"

##### Populate configuration files ##### 

# Single run - cpu run - human_genomics_pipeline
for config_dir in exome/cpu_run/*/single_run/*/human_genomics_pipeline/config; do
    cp $workingdir/config_single_cpu_hgp/config.yaml $config_dir
done

# Single run - gpu run - human_genomics_pipeline
for config_dir in exome/gpu_run/*/single_run/*/human_genomics_pipeline/config; do
    cp $workingdir/config_single_gpu_hgp/config.yaml $config_dir
done

# Cohort run - cpu run - human_genomics_pipeline
for config_dir in exome/cpu_run/*/cohort_run/*/human_genomics_pipeline/config; do
    cp $workingdir/config_cohort_cpu_hgp/config.yaml $config_dir
done

# Cohort run - gpu run - human_genomics_pipeline
for config_dir in exome/gpu_run/*/cohort_run/*/human_genomics_pipeline/config; do
    cp $workingdir/config_cohort_gpu_vap/config.yaml $config_dir
done

# Single run - cpu run - vcf_annotation_pipeline
for config_dir in exome/cpu_run/*/single_run/*/vcf_annotation_pipeline/config; do
    cp $workingdir/config_single_cpu_vap/config.yaml $config_dir
done

# Single run - gpu run - vcf_annotation_pipeline
for config_dir in exome/gpu_run/*/single_run/*/vcf_annotation_pipeline/config; do
    cp $workingdir/config_single_gpu_vap/config.yaml $config_dir
done

# Cohort run - cpu run - vcf_annotation_pipeline
for config_dir in exome/cpu_run/*/cohort_run/*/vcf_annotation_pipeline/config; do
    cp $workingdir/config_cohort_cpu_vap/config.yaml $config_dir
done

# Cohort run - gpu run - vcf_annotation_pipeline
for config_dir in exome/gpu_run/*/cohort_run/*/vcf_annotation_pipeline/config; do
    cp $workingdir/config_cohort_gpu_vap/config.yaml $config_dir
done

##### Populate cluster configuration files #####

# All runs - human_genomics_pipeline
for config_dir in exome/*/*/*/*/human_genomics_pipeline/config; do
    cp $workingdir/cluster_config/cluster.json $config_dir
done

# All runs - vcf_annotation_pipeline
for config_dir in exome/*/*/*/*/vcf_annotation_pipeline/config; do
    cp $workingdir/cluster_config/cluster.json $config_dir
done

##### Populate runscripts #####

# All runs - human_genomics_pipeline
for config_dir in exome/*/*/*/*/human_genomics_pipeline/workflow; do
    cp $workingdir/runscript_hgp/run_hpc.sh $config_dir
done

# All runs - vcf_annotation_pipeline
for config_dir in exome/*/*/*/*/vcf_annotation_pipeline/workflow; do
    cp $workingdir/runscript_vap/run_hpc.sh $config_dir
done

##### Populate rule files #####

# 01 thread - human_genomics_pipeline
for config_dir in exome/*/*/*/01_threads/human_genomics_pipeline/workflow/rules; do
    cp $workingdir/rules_01_thread_hgp/* $config_dir
done

# 02 threads - human_genomics_pipeline
for config_dir in exome/*/*/*/02_threads/human_genomics_pipeline/workflow/rules; do
    cp $workingdir/rules_02_thread_hgp/* $config_dir
done

# 04 threads - human_genomics_pipeline
for config_dir in exome/*/*/*/04_threads/human_genomics_pipeline/workflow/rules; do
    cp $workingdir/rules_04_thread_hgp/* $config_dir
done

# 08 threads - human_genomics_pipeline
for config_dir in exome/*/*/*/08_threads/human_genomics_pipeline/workflow/rules; do
    cp $workingdir/rules_08_thread_hgp/* $config_dir
done

# 16 threads - human_genomics_pipeline
for config_dir in exome/*/*/*/16_threads/human_genomics_pipeline/workflow/rules; do
    cp $workingdir/rules_16_thread_hgp/* $config_dir
done

# 32 threads - human_genomics_pipeline
for config_dir in exome/*/*/*/32_threads/human_genomics_pipeline/workflow/rules; do
    cp $workingdir/rules_32_thread_hgp/* $config_dir
done

# 01 thread - vcf_annotation_pipeline
for config_dir in exome/*/*/*/01_threads/vcf_annotation_pipeline/workflow/rules; do
    cp $workingdir/rules_01_thread_vap/* $config_dir
done

# 02 threads - vcf_annotation_pipeline
for config_dir in exome/*/*/*/02_threads/vcf_annotation_pipeline/workflow/rules; do
    cp $workingdir/rules_02_thread_vap/* $config_dir
done

# 04 threads - vcf_annotation_pipeline
for config_dir in exome/*/*/*/04_threads/vcf_annotation_pipeline/workflow/rules; do
    cp $workingdir/rules_04_thread_vap/* $config_dir
done

# 08 threads - vcf_annotation_pipeline
for config_dir in exome/*/*/*/08_threads/vcf_annotation_pipeline/workflow/rules; do
    cp $workingdir/rules_08_thread_vap/* $config_dir
done

# 16 threads - vcf_annotation_pipeline
for config_dir in exome/*/*/*/16_threads/vcf_annotation_pipeline/workflow/rules; do
    cp $workingdir/rules_16_thread_vap/* $config_dir
done

# 32 threads - vcf_annotation_pipeline
for config_dir in exome/*/*/*/32_threads/vcf_annotation_pipeline/workflow/rules; do
    cp $workingdir/rules_32_thread_vap/* $config_dir
done