#!/bin/bash

# Move vcf files output human_genomics_pipeline
for vcfdir in exome/*/*/*/*/vcf; do
    cp $vcfdir/../human_genomics_pipeline/results/called/* $vcfdir
done

# Move bam files output human_genomics_pipeline
for bamdir in exome/*/*/*/*/bam; do
    cp $bamdir/../human_genomics_pipeline/results/mapped/* $bamdir
done
