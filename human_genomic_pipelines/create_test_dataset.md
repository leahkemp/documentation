# Create test dataset for pipelines

Created: 2020/10/29 10:07:39
Last modified: 2020/11/09 16:23:01

- **Aim:** Create a test dataset for [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) and [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline) (a trio and a singleton)
- **Prerequisite software:**  [Conda 4.8.5](https://docs.conda.io/projects/conda/en/latest/index.html)
- **OS:** CentOS-7 (ORAC - ESR cluster)

## Table of contents

- [Create test dataset for pipelines](#create-test-dataset-for-pipelines)
  - [Table of contents](#table-of-contents)
  - [Setup](#setup)
  - [Get bams and vcf for publicly available trio](#get-bams-and-vcf-for-publicly-available-trio)
  - [Manually create pedigree file](#manually-create-pedigree-file)
  - [Reduce dataset - subset by exome capture regions](#reduce-dataset---subset-by-exome-capture-regions)
  - [Process vcf](#process-vcf)
    - [SnpSift_filter_proband](#snpsift_filter_proband)
  - [Randomly sub-sample variants](#randomly-sub-sample-variants)
  - [Create a bed file from vcf](#create-a-bed-file-from-vcf)
  - [Pull out fastq reads from bam](#pull-out-fastq-reads-from-bam)
  - [Test running test data through pipelines](#test-running-test-data-through-pipelines)
  - [Other](#other)


## Setup

```bash
mkdir create_test_dataset
cd create_test_dataset
```

## Get bams and vcf for publicly available trio

ChineseTrio

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

See [here](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/analysis/RTG_RTGJointTrio_06062019/GRCh37/README.md) for information on this data

## Manually create pedigree file

Label manually created pedigree file as `NA24694_pedigree.ped` and put in the `.pedigrees/` directory

## Reduce dataset - subset by exome capture regions

Create conda env with GATK4 installed

```bash
conda create -n gatk4 python=3.7.6
conda install -c bioconda gatk4=4.1.9.0
conda activate gatk4
```

Subset vcf by custom exome capture regions

```bash
gatk SelectVariants \
-V family.merged.vcf.gz \
-O family.merged.subset.vcf \
-L custom.bed.gz
```

This takes the vcf file from 6,226,947 variants to 329,619 variants

## Process vcf

- Data/vcf is already filtered

### SnpSift_filter_proband

Filter for variants only found in the proband

Code based on [this version of the pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline/blob/c82f0749a705ebc85d9589cd2b84f014cc8b735f/workflow/rules/SnpSift_filter_proband.smk)

```bash
conda create -n snpsift python=3.7.6
conda activate snpsift
conda install -c bioconda snpsift=4.3.1t

cat family.merged.subset.vcf | SnpSift filter '( isVariant ( GEN[NA24631] ) ) ' > family.merged.subset.probandonly.vcf
```

Now 255,194 variants

## Randomly sub-sample variants

50% of the variants

```bash
conda activate gatk4

gatk SelectVariants \
-V family.merged.subset.probandonly.vcf \
-O family.merged.subset.probandonly.randomsubset.vcf \
--select-random-fraction 0.5
```

Now 1246 variants throughout these chromosomes:

```bash
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22
X
Y
```

## Create a bed file from vcf

Get the positions of the variants - create a bed file

```bash
conda create -n bedops python=3.7.6
conda activate bedops
conda install -c bioconda bedops=2.4.39

vcf2bed < family.merged.subset.probandonly.randomsubset.vcf > family.merged.subset.probandonly.randomsubset.bed
```

Re-sort bed

```bash
sort -V -k1,1 -k2,2 family.merged.subset.probandonly.randomsubset.bed > family.merged.subset.probandonly.randomsubset.sorted.bed
```

Add chr prefix (so it is compatible with bams downstream, note. chrM/MT doesn't need to be fixed, not present in bed)

```bash
awk 'OFS="\t" {if (NR > 5) $1="chr"$1; print}' family.merged.subset.probandonly.randomsubset.sorted.bed > family.merged.subset.probandonly.randomsubset.sorted.chr.bed
```

## Pull out fastq reads from bam

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

Check that the number of unmapped and mapped reads total the number of reads in the original bam file (for the reduced bams)

```bash
for bam in "${bams_to_process[@]}"; do
samtools flagstat $bam.reduced.bam
done
```

```bash
1790779 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
13026 + 0 supplementary
0 + 0 duplicates
1790648 + 0 mapped (99.99% : N/A)
1777753 + 0 paired in sequencing
895801 + 0 read1
881952 + 0 read2
1667792 + 0 properly paired (93.81% : N/A)
1765307 + 0 with itself and mate mapped
12315 + 0 singletons (0.69% : N/A)
70943 + 0 with mate mapped to a different chr
55191 + 0 with mate mapped to a different chr (mapQ>=5)

2034368 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
12248 + 0 supplementary
0 + 0 duplicates
2034223 + 0 mapped (99.99% : N/A)
2022120 + 0 paired in sequencing
1018632 + 0 read1
1003488 + 0 read2
1826298 + 0 properly paired (90.32% : N/A)
2009917 + 0 with itself and mate mapped
12058 + 0 singletons (0.60% : N/A)
81082 + 0 with mate mapped to a different chr
63847 + 0 with mate mapped to a different chr (mapQ>=5)

1902482 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
11682 + 0 supplementary
0 + 0 duplicates
1902354 + 0 mapped (99.99% : N/A)
1890800 + 0 paired in sequencing
952414 + 0 read1
938386 + 0 read2
1709973 + 0 properly paired (90.44% : N/A)
1878819 + 0 with itself and mate mapped
11853 + 0 singletons (0.63% : N/A)
71685 + 0 with mate mapped to a different chr
55869 + 0 with mate mapped to a different chr (mapQ>=5)
```

Compare this to the number of mapped and unmapped reads

```bash
for bam in "${bams_to_process[@]}"; do
# Mapped reads
samtools view -c $bam.reduced_mapped.sort
# Unmapped reads
samtools view -c $bam.reduced_unmapped.sort
done
```

My output:

```bash
1777858
12921
2021739
12629
1890089
12393
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
echo $(zcat fastq/$fastq | wc -l)/4| bc
done
```

My output:

```bash
27730
27730
27758
27758
23746
23746
```

## Test running test data through pipelines

Setup

```bash
# Clone pipelines
git clone https://github.com/ESR-NZ/human_genomics_pipeline.git
cd human_genomics_pipeline
git checkout a4dc43d7557df41f95f8f5963642b53f5cded99e
cd ..
git clone https://github.com/ESR-NZ/vcf_annotation_pipeline.git
cd vcf_annotation_pipeline
git checkout 5390d9320030e6dcba2596395782667d34cc0a9e
```

Manually configure pipelines

Run pipelines

```bash
screen -S test_data_test

conda activate pipeline_run_env

cd ../human_genomics_pipeline/workflow/
bash dryrun_hpc.sh
bash run_hpc.sh
```

```bash
cd ../../vcf_annotation_pipeline/workflow/
bash dryrun_hpc.sh
bash run_hpc.sh
```

## Other

Provide custom exome capture file and bed file

Random seed for random sampling of variants - re-generating this file might produce different results due to randomisation
