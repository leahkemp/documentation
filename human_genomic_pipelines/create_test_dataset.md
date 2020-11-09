# Create test dataset for pipelines

Created: 2020/10/29 10:07:39
Last modified: 2020/11/10 12:13:05

- **Aim:** Create a test dataset for [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) and [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline)
- **Prerequisite software:**  [Conda 4.8.5](https://docs.conda.io/projects/conda/en/latest/index.html)
- **OS:** CentOS-7 (ORAC - ESR cluster)

## Table of contents

- [Create test dataset for pipelines](#create-test-dataset-for-pipelines)
  - [Table of contents](#table-of-contents)
  - [Setup](#setup)
  - [Get bams and vcf for publicly available trio](#get-bams-and-vcf-for-publicly-available-trio)
  - [Manually create pedigree file](#manually-create-pedigree-file)
  - [Reduce vcf](#reduce-vcf)
    - [Subset by exome capture regions](#subset-by-exome-capture-regions)
    - [Filter for variants only found in the proband](#filter-for-variants-only-found-in-the-proband)
    - [Randomly sub-sample variants](#randomly-sub-sample-variants)
  - [Create a bed file from the reduced vcf](#create-a-bed-file-from-the-reduced-vcf)
  - [Pull out fastq reads from bams](#pull-out-fastq-reads-from-bams)
  - [Run test data through pipelines](#run-test-data-through-pipelines)
    - [Run cohort test data through pipelines](#run-cohort-test-data-through-pipelines)
    - [Run single sample test data through pipelines](#run-single-sample-test-data-through-pipelines)


## Setup

```bash
mkdir create_test_dataset
cd create_test_dataset
```

## Get bams and vcf for publicly available trio

```bash
# vcf
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/analysis/RTG_RTGJointTrio_06062019/GRCh37/family.merged.vcf.gz

# bams
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG005_NA24631_son/NIST_Stanford_Illumina_6kb_matepair/bams/HG005.mate_pair.sorted.bam
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG005_NA24631_son/NIST_Stanford_Illumina_6kb_matepair/bams/HG005.mate_pair.sorted.bam.bai
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG006_NA24694-huCA017E_father/NIST_Stanford_Illumina_6kb_matepair/bams/HG006.mate_pair.sorted.bam
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG006_NA24694-huCA017E_father/NIST_Stanford_Illumina_6kb_matepair/bams/HG006.mate_pair.sorted.bam.bai
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG007_NA24695-hu38168_mother/NIST_Stanford_Illumina_6kb_matepair/bams/HG007.mate_pair.sorted.bam
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG007_NA24695-hu38168_mother/NIST_Stanford_Illumina_6kb_matepair/bams/HG007.mate_pair.sorted.bam.bai
```

See [here](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/analysis/RTG_RTGJointTrio_06062019/GRCh37/README.md) for information on this vcf data, as well as [HG005_NA24631_son](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG005_NA24631_son/NIST_Stanford_Illumina_6kb_matepair/README.NIST_Stanford_Illumina_6kb_matepair), [HG006_NA24694-huCA017E_father](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG006_NA24694-huCA017E_father/NIST_Stanford_Illumina_6kb_matepair/README.NIST_Stanford_Illumina_6kb_matepair) and [HG007_NA24695-hu38168_mother](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG007_NA24695-hu38168_mother/NIST_Stanford_Illumina_6kb_matepair/README.NIST_Stanford_Illumina_6kb_matepair) for information on this bam data.

## Manually create pedigree file

Label manually created pedigree file as `NA24631_pedigree.ped` and put in the `.pedigrees/` directory

## Reduce vcf

### Subset by exome capture regions

Create conda env with GATK4 installed

```bash
conda create -n gatk4 python=3.7.6
conda install -c bioconda gatk4=4.1.9.0
conda activate gatk4
```

Subset vcf by custom exome capture regions (file can be found [here](./create_test_dataset_files/custom.bed.gz))

```bash
gatk SelectVariants \
-V family.merged.vcf.gz \
-O family.merged.subset.vcf \
-L custom.bed.gz
```

This takes the vcf file from 6,226,947 variants to 329,619 variants

### Filter for variants only found in the proband

Code based on [this version of the pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline/blob/c82f0749a705ebc85d9589cd2b84f014cc8b735f/workflow/rules/SnpSift_filter_proband.smk)

```bash
conda create -n snpsift python=3.7.6
conda activate snpsift
conda install -c bioconda snpsift=4.3.1t

cat family.merged.subset.vcf | SnpSift filter '( isVariant ( GEN[NA24631] ) ) ' > family.merged.subset.probandonly.vcf
```

### Randomly sub-sample variants

```bash
conda activate gatk4

gatk SelectVariants \
-V family.merged.subset.probandonly.vcf \
-O family.merged.subset.probandonly.randomsubset.vcf \
--select-random-fraction 0.5
```

## Create a bed file from the reduced vcf

Get the positions of the variants - create a bed file

```bash
conda create -n bedops python=3.7.6
conda activate bedops
conda install -c bioconda bedops=2.4.39

vcf2bed < family.merged.subset.probandonly.randomsubset.vcf > family.merged.subset.probandonly.randomsubset.bed

# Re-sort bed
sort -V -k1,1 -k2,2 family.merged.subset.probandonly.randomsubset.bed > family.merged.subset.probandonly.randomsubset.sorted.bed

# Add chr prefix (so it is compatible with bams downstream, note. chrM/MT doesn't need to be fixed, not present in bed)
awk 'OFS="\t" {if (NR > 5) $1="chr"$1; print}' family.merged.subset.probandonly.randomsubset.sorted.bed > family.merged.subset.probandonly.randomsubset.sorted.chr.bed
```

*Note. the resulting `family.merged.subset.probandonly.randomsubset.sorted.chr.bed` file can be found [here](./create_test_dataset_files/family.merged.subset.probandonly.randomsubset.sorted.chr.bed)*

## Pull out fastq reads from bams

Use these variant positions (defined in the bed file) to extract the fastq reads from the bam files (following [this guide](https://gist.github.com/darencard/72ddd9e6c08aaff5ff64ca512a04a6dd))

```bash
conda create -n bedtools python=3.7.6
conda activate bedtools
conda install -c bioconda bedtools=2.29.2
conda install -c bioconda samtools=1.9

bams_to_process=("HG005.mate_pair.sorted" "HG006.mate_pair.sorted" "HG007.mate_pair.sorted")

for bam in "${bams_to_process[@]}"; do
# Subset bam with bed file
samtools view $bam.bam -L family.merged.subset.probandonly.randomsubset.sorted.chr.bed -@ 16 -b > $bam.reduced.bam
done
```

Start extracting FASTQ reads from bam

```bash
for bam in "${bams_to_process[@]}"; do
# R1 mapped, R2 mapped
samtools view -u -f 1 -F 12 $bam.reduced.bam > $bam.reduced_map_map.bam
# R1 unmapped, R2 mapped
samtools view -u -f 4 -F 264 $bam.reduced.bam > $bam.reduced_unmap_map.bam
# R1 mapped, R2 unmapped
samtools view -u -f 8 -F 260 $bam.reduced.bam > $bam.reduced_map_unmap.bam
# R1 unmapped, R2 unmapped
samtools view -u -f 12 -F 256 $bam.reduced.bam > $bam.reduced_unmap_unmap.bam
# Merge the three files that contain at least one unmapped pair
samtools merge -u $bam.reduced_unmapped.bam $bam.reduced_unmap_map.bam $bam.reduced_map_unmap.bam $bam.reduced_unmap_unmap.bam
# Sort
samtools sort -n $bam.reduced_map_map.bam -o $bam.reduced_mapped.sort
samtools sort -n $bam.reduced_unmapped.bam -o $bam.reduced_unmapped.sort
done
```

Create FASTQ read files

```bash
for bam in "${bams_to_process[@]}"; do
# Extract the mapped FASTQ reads into two paired read files
bamToFastq -i $bam.reduced_mapped.sort -fq $bam.reduced_mapped.1.fastq -fq2 $bam.reduced_mapped.2.fastq
bamToFastq -i $bam.reduced_unmapped.sort -fq $bam.reduced_unmapped.1.fastq -fq2 $bam.reduced_unmapped.2.fastq
# Combine both the first and paired reads together from the mapped and unmapped files
cat $bam.reduced_mapped.1.fastq $bam.reduced_unmapped.1.fastq > $bam.reduced_R1.fastq
cat $bam.reduced_mapped.2.fastq $bam.reduced_unmapped.2.fastq > $bam.reduced_R2.fastq
# Bgzip
bgzip $bam.reduced_R1.fastq
bgzip $bam.reduced_R2.fastq
done
```

Create file directory structure and rename files to work as inputs for pipelines

```bash
mkdir fastq

# Son
mv HG005.mate_pair.sorted.reduced_R1.fastq.gz fastq/NA24631_1.fastq.gz
mv HG005.mate_pair.sorted.reduced_R2.fastq.gz fastq/NA24631_2.fastq.gz

# Father
mv HG006.mate_pair.sorted.reduced_R1.fastq.gz fastq/NA24694_1.fastq.gz
mv HG006.mate_pair.sorted.reduced_R2.fastq.gz fastq/NA24694_2.fastq.gz

# Mother
mv HG007.mate_pair.sorted.reduced_R1.fastq.gz fastq/NA24695_1.fastq.gz
mv HG007.mate_pair.sorted.reduced_R2.fastq.gz fastq/NA24695_2.fastq.gz
```

Check the number of reads in each fastq file is consistent between reads 1 and 2 for all samples

```bash
fastqs_to_process=("NA24631_1.fastq.gz" "NA24631_2.fastq.gz" "NA24694_1.fastq.gz" "NA24694_2.fastq.gz" "NA24695_1.fastq.gz" "NA24695_2.fastq.gz")

for fastq in "${fastqs_to_process[@]}"; do
echo "Reads in $fastq:"
echo $(zcat fastq/$fastq | wc -l)/4| bc
done
```

My output:

```bash
Reads in NA24631_1.fastq.gz:
27730
Reads in NA24631_2.fastq.gz:
27730
Reads in NA24694_1.fastq.gz:
27758
Reads in NA24694_2.fastq.gz:
27758
Reads in NA24695_1.fastq.gz:
23746
Reads in NA24695_2.fastq.gz:
23746
```

## Run test data through pipelines

Setup

```bash
# Clone pipelines
git clone https://github.com/ESR-NZ/human_genomics_pipeline.git
cd human_genomics_pipeline
git checkout v1.0.0
cd ..
git clone https://github.com/ESR-NZ/vcf_annotation_pipeline.git
cd vcf_annotation_pipeline
git checkout v1.0.0
```

Manually configure pipelines

### Run cohort test data through pipelines

Manually set `DATA:` in `human_genomics_pipeline/config/config.yaml` and `vcf_annotation_pipeline/config/config.yaml` to `DATA: "Cohort"`

Run pipelines

```bash
screen -S test_data_test

conda activate pipeline_run_env

cd ../human_genomics_pipeline/workflow/
bash dryrun_hpc.sh
bash run_hpc.sh
```

Integrate cohort test data with github repos

### Run single sample test data through pipelines

Manually set `DATA:` in `human_genomics_pipeline/config/config.yaml` and `vcf_annotation_pipeline/config/config.yaml` to `DATA: "Single"`

Remove previous run results, metadata and unnecessary input data

```bash
rm ../../fastq/NA24694*
rm ../../fastq/NA24695*

rm -r ../results/*
rm slurm*
rm -r logs
rm -r benchmarks
```

Run pipelines

```bash
bash dryrun_hpc.sh
bash run_hpc.sh
```

Integrate single sample test data with github repos
