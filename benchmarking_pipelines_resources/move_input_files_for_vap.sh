#!/bin/bash

# Move vcf files output human_genomics_pipeline
for vcfdir in exome/*/*/*/*/vcf; do
    cp $vcfdir/../human_genomics_pipeline/results/called/* $vcfdir
done

# Move bam files output human_genomics_pipeline
for bamdir in exome/*/*/*/*/bams; do
    cp $bamdir/../human_genomics_pipeline/results/mapped/*_recalibrated.bai $bamdir
    cp $bamdir/../human_genomics_pipeline/results/mapped/*_recalibrated.bam $bamdir
done
