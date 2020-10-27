
- [Context](#context)
- [Check files don't load into scout (on Wintermute)](#check-files-dont-load-into-scout-on-wintermute)
- [Reduce dataset](#reduce-dataset)
- [Run through pipeline on ESR cluster](#run-through-pipeline-on-esr-cluster)
- [Run through pipeline on research server (Wintermute)](#run-through-pipeline-on-research-server-wintermute)
- [2. Check for any differences between the dbNSFP databases between Wintermute and the ESR cluster](#2-check-for-any-differences-between-the-dbnsfp-databases-between-wintermute-and-the-esr-cluster)
- [4. Re-move the dbNSFP database from Wintermute to the ESR cluster - see if a fresh move will fix the error](#4-re-move-the-dbnsfp-database-from-wintermute-to-the-esr-cluster---see-if-a-fresh-move-will-fix-the-error)

## Context

I'm having some issues with loading scout cases when analysing clinical genomic data through the pipelines. It appears to be working when I run the pipelines manually (by running the code for each step manually), but running the pipelines manually is confounded with running the pipelines on a single research server (Wintermute) rather than the ESR cluster. I fixed a previous issue with the dbNSFP annotation step running in the pipeline on the ESR cluster - the datatype of each dbNSFP annotation weren't being set properly until I ensured the file (dbNSFPv4.0a.hg19.custombuild.gz.data_types) was sitting alongside my dbNSFP database (which it *was* on Wintermute, and *wasn't* on the ESR cluster). This fixed the datatypes in the dbNSFP annotations in the final annotated vcf. I though this was going to fix the scout loading error, but low and behold it doesn't. Hence the error has is now being referred to as dumb bug.

Error when loading case config into scout:

```bash
2020-10-27 10:27:27 Wintermute scout.adapter.mongo.variant_loader[7030] ERROR unexpected error
Traceback (most recent call last):
  File "/store/lkemp/GA_clinical_genomics/scout/scout/adapter/mongo/variant_loader.py", line 688, in load_variants
    sample_info=sample_info,
  File "/store/lkemp/GA_clinical_genomics/scout/scout/adapter/mongo/variant_loader.py", line 438, in _load_variants
    category=category,
  File "/store/lkemp/GA_clinical_genomics/scout/scout/parse/variant/variant.py", line 307, in parse_variant
    parsed_variant["conservation"] = parse_conservations(variant, parsed_transcripts)
  File "/store/lkemp/GA_clinical_genomics/scout/scout/parse/variant/conservation.py", line 36, in parse_conservations
    result = parse_conservation_info(variant, value, key)
  File "/store/lkemp/GA_clinical_genomics/scout/scout/parse/variant/conservation.py", line 65, in parse_conservation_info
    if score >= CONSERVATION[field_key]["conserved_min"]:
TypeError: '>=' not supported between instances of 'str' and 'float'
2020-10-27 10:27:27 Wintermute scout.adapter.mongo.variant_loader[7030] WARNING Deleting inserted variants
2020-10-27 10:27:27 Wintermute scout.adapter.mongo.variant[7030] INFO Deleting old clinical  variants for case DHB4001
2020-10-27 10:27:27 Wintermute scout.adapter.mongo.variant[7030] INFO 58 variants deleted
2020-10-27 10:27:27 Wintermute scout.commands.load.case[7030] ERROR Something went wrong during loading
2020-10-27 10:27:27 Wintermute scout.commands.load.case[7030] WARNING '>=' not supported between instances of 'str' and 'float'
Aborted!
```

See more info on this dumb bug on my [jira](https://leah-kemp.atlassian.net/jira/software/projects/DHBG/boards/1?selectedIssue=DHBG-214)

## Check files don't load into scout (on Wintermute)

```bash
screen -r -S scout_serve
scout -db dhb-database load case /store/lkemp/GA_clinical_genomics/scout_case_configs/DHB4001_config.yaml
```

Still getting the error

## Reduce dataset

I'll reduce the dataset to speed up the debugging process. To do this, I'll randomly subsample the fastq reads after mapping to the reference genome.

Setup

```bash
cd /NGS/scratch/KSCBIOM/HumanGenomics/
mkdir fix_dumb_bug
cd fix_dumb_bug
mkdir deploy_to_cluster
cd deploy_to_cluster
screen -S fix_dumb_bug
conda activate pipeline_run_env
```

Copy previously run human_genomics_pipeline output to run vcf_annotation_pipeline from

```bash
cd /NGS/scratch/KSCBIOM/HumanGenomics/fix_dumb_bug/deploy_to_cluster/
mkdir -p human_genomics_pipeline/results/called/
cp /NGS/scratch/KSCBIOM/HumanGenomics/DHB_WGS/DHB4177/DHB4177_GRCh37_WGS_trio.combined.gpu.genotyped.vcf.gz ./human_genomics_pipeline/results/called/
gunzip ./human_genomics_pipeline/results/called/DHB4177_GRCh37_WGS_trio.combined.gpu.genotyped.vcf.gz

mkdir pedigrees
cp /NGS/scratch/KSCBIOM/HumanGenomics/DHB_WGS/DHB4177/pedigree_4177.ped ./pedigrees/DHB4177_pedigree.ped
```

Manually edit pedigree file (family name)

Get gatk

```bash
conda create -n gatk4 python=3.7.6
conda activate gatk4
conda install -c bioconda gatk4=4.1.9.0
```

Randomly subsample dataset (0.8% of genome)

```bash
gatk SelectVariants \
-R /NGS/scratch/KSCBIOM/HumanGenomics/publicData/b37/human_g1k_v37_decoy.fasta \
-V ./human_genomics_pipeline/results/called/DHB4177_GRCh37_WGS_trio.combined.gpu.genotyped.vcf \
-O ./human_genomics_pipeline/results/called/DHB4177_raw_snps_indels.g.vcf \
--select-random-fraction 0.008
```

```bash
rm ./human_genomics_pipeline/results/called/DHB4177_GRCh37_WGS_trio.combined.gpu.genotyped.vcf
```

Get rough number of variants in the reduced dataset

```bash
grep -v '#' ./human_genomics_pipeline/results/called/DHB4177_raw_snps_indels.g.vcf | wc -l
```

Output

```bash
55157
```

Get coverage of chromosomes

```bash
grep -v '#' ./human_genomics_pipeline/results/called/DHB4177_raw_snps_indels.g.vcf | awk '{print $1}' | uniq
```

Output

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
GL000226.1
GL000229.1
GL000235.1
GL000201.1
GL000247.1
GL000245.1
GL000203.1
GL000248.1
GL000244.1
GL000202.1
GL000234.1
GL000232.1
GL000240.1
GL000241.1
GL000230.1
GL000237.1
GL000233.1
GL000198.1
GL000208.1
GL000191.1
GL000227.1
GL000228.1
GL000214.1
GL000221.1
GL000209.1
GL000218.1
GL000220.1
GL000211.1
GL000199.1
GL000217.1
GL000216.1
GL000205.1
GL000219.1
GL000224.1
GL000195.1
GL000212.1
GL000222.1
GL000193.1
GL000194.1
GL000225.1
GL000192.1
hs37d5
```

## Run through pipeline on ESR cluster

Get pipelines

```bash
git clone https://github.com/ESR-NZ/vcf_annotation_pipeline.git
cd vcf_annotation_pipeline
git checkout 75774b4ffef3258d0c7a5cff13072971e8553104
```

Manually configure pipeline

Run pipelines on reduced dataset

```bash
cd /NGS/scratch/KSCBIOM/HumanGenomics/fix_dumb_bug/deploy_to_cluster/vcf_annotation_pipeline/workflow/
bash dryrun_hpc.sh
bash run_hpc.sh
```

Transfer files to Wintermute to test loading into scout (catch intermediate files before they are removed by the pipeline)

```bash
# Run this code from wintermute
scp -r orac:/NGS/scratch/KSCBIOM/HumanGenomics/fix_dumb_bug/deploy_to_cluster/vcf_annotation_pipeline/results/filtered/DHB4177_filtered.vcf /store/lkemp/fix_dumb_bug/pipeline_cluster_output/
scp -r orac:/NGS/scratch/KSCBIOM/HumanGenomics/fix_dumb_bug/deploy_to_cluster/vcf_annotation_pipeline/results/annotated/DHB4177_filtered_dbnsfp.vcf.gz /store/lkemp/fix_dumb_bug/pipeline_cluster_output/
scp -r orac:/NGS/scratch/KSCBIOM/HumanGenomics/fix_dumb_bug/deploy_to_cluster/vcf_annotation_pipeline/results/annotated/DHB4177_filtered_dbnsfp_vep.vcf.gz /store/lkemp/fix_dumb_bug/pipeline_cluster_output/
scp -r orac:/NGS/scratch/KSCBIOM/HumanGenomics/fix_dumb_bug/deploy_to_cluster/vcf_annotation_pipeline/results/annotated/DHB4177_filtered_dbnsfp_vep_cadd_dbsnp.vcf /store/lkemp/fix_dumb_bug/pipeline_cluster_output/
```

Move bams

```bash
scp -r orac:/NGS/scratch/KSCBIOM/HumanGenomics/DHB_WGS/DHB4177/DHB4177_GRCh37_WGS.bam /store/lkemp/fix_dumb_bug/bams/
scp -r orac:/NGS/scratch/KSCBIOM/HumanGenomics/DHB_WGS/DHB4177/DHB4177_GRCh37_WGS.bam.bai /store/lkemp/fix_dumb_bug/bams/
scp -r orac:/NGS/scratch/KSCBIOM/HumanGenomics/DHB_WGS/DHB4177/DHB4178_GRCh37_WGS.bam /store/lkemp/fix_dumb_bug/bams/
scp -r orac:/NGS/scratch/KSCBIOM/HumanGenomics/DHB_WGS/DHB4177/DHB4178_GRCh37_WGS.bam.bai /store/lkemp/fix_dumb_bug/bams/
scp -r orac:/NGS/scratch/KSCBIOM/HumanGenomics/DHB_WGS/DHB4177/DHB4179_GRCh37_WGS.bam /store/lkemp/fix_dumb_bug/bams/
scp -r orac:/NGS/scratch/KSCBIOM/HumanGenomics/DHB_WGS/DHB4177/DHB4179_GRCh37_WGS.bam.bai /store/lkemp/fix_dumb_bug/bams/
```

Filter multiallelic sites

```bash
conda activate bcftools

bcftools view -O z --max-alleles 2 --exclude-types indels \
/store/lkemp/fix_dumb_bug/pipeline_cluster_output/DHB4177_filtered.vcf > /store/lkemp/fix_dumb_bug/pipeline_cluster_output/DHB4177_filtered_cleaned.vcf

bcftools view -O z --max-alleles 2 --exclude-types indels \
/store/lkemp/fix_dumb_bug/pipeline_cluster_output/DHB4177_filtered_dbnsfp.vcf.gz > /store/lkemp/fix_dumb_bug/pipeline_cluster_output/DHB4177_filtered_dbnsfp_cleaned.vcf.gz

bcftools view -O z --max-alleles 2 --exclude-types indels \
/store/lkemp/fix_dumb_bug/pipeline_cluster_output/DHB4177_filtered_dbnsfp_vep.vcf.gz > /store/lkemp/fix_dumb_bug/pipeline_cluster_output/DHB4177_filtered_dbnsfp_vep_cleaned.vcf.gz

bcftools view -O z --max-alleles 2 --exclude-types indels \
/store/lkemp/fix_dumb_bug/pipeline_cluster_output/DHB4177_filtered_dbnsfp_vep_cadd_dbsnp.vcf > /store/lkemp/fix_dumb_bug/pipeline_cluster_output/DHB4177_filtered_dbnsfp_vep_cadd_dbsnp_cleaned.vcf
```

Test loading into scout

```bash
scout -db demo-database load case /store/lkemp/fix_dumb_bug/DHB4177_config.yaml

# Run through pipeline on the ESR cluster
DHB4177_filtered_cleaned.vcf # works
DHB4177_filtered_dbnsfp_cleaned.vcf.gz # throws error
DHB4177_filtered_dbnsfp_vep_cleaned.vcf.gz # throws error
DHB4177_filtered_dbnsfp_vep_cadd_dbsnp_cleaned.vcf # throws error
```

## Run through pipeline on research server (Wintermute)

```bash
cd /store/lkemp/
mkdir fix_dumb_bug
cd fix_dumb_bug
screen -S fix_dumb_bug
```

Copy reduced fastq data to working directory (transfer from cluster to Wintermute)

```bash
mkdir -p human_genomics_pipeline/results/called/
mkdir pedigrees

scp -r orac:/NGS/scratch/KSCBIOM/HumanGenomics/fix_dumb_bug/deploy_to_cluster/pedigrees/DHB4177_pedigree.ped /store/lkemp/fix_dumb_bug/pedigrees/
scp -r orac:/NGS/scratch/KSCBIOM/HumanGenomics/fix_dumb_bug/deploy_to_cluster/human_genomics_pipeline/results/called/* /store/lkemp/fix_dumb_bug/human_genomics_pipeline/results/called/
```

Get pipelines

```bash
git clone https://github.com/ESR-NZ/vcf_annotation_pipeline.git
cd vcf_annotation_pipeline
git checkout 75774b4ffef3258d0c7a5cff13072971e8553104
```

Manually configure pipeline

Run pipelines on reduced dataset

```bash
cd /NGS/scratch/KSCBIOM/HumanGenomics/fix_dumb_bug/vcf_annotation_pipeline/workflow/
bash dryrun_hpc.sh
bash run_hpc.sh
```

Transfer files to different dir in Wintermute to test loading into scout (catch intermediate files before they are removed by the pipeline)

```bash
# Run this code from wintermute
cp /store/lkemp/fix_dumb_bug/vcf_annotation_pipeline/results/filtered/DHB4177_filtered.vcf /store/lkemp/fix_dumb_bug/pipeline_wintermute_output/
cp /store/lkemp/fix_dumb_bug/vcf_annotation_pipeline/results/annotated/DHB4177_filtered_dbnsfp.vcf.gz /store/lkemp/fix_dumb_bug/pipeline_wintermute_output/
cp /store/lkemp/fix_dumb_bug/vcf_annotation_pipeline/results/annotated/DHB4177_filtered_dbnsfp_vep.vcf.gz /store/lkemp/fix_dumb_bug/pipeline_wintermute_output/
cp /store/lkemp/fix_dumb_bug/vcf_annotation_pipeline/results/annotated/DHB4177_filtered_dbnsfp_vep_cadd_dbsnp.vcf /store/lkemp/fix_dumb_bug/pipeline_wintermute_output/

```

Filter multiallelic sites

```bash
conda activate bcftools

bcftools view -O z --max-alleles 2 --exclude-types indels \
/store/lkemp/fix_dumb_bug/pipeline_wintermute_output/DHB4177_filtered.vcf > /store/lkemp/fix_dumb_bug/pipeline_wintermute_output/DHB4177_filtered_cleaned.vcf

bcftools view -O z --max-alleles 2 --exclude-types indels \
/store/lkemp/fix_dumb_bug/pipeline_wintermute_output/DHB4177_filtered_dbnsfp.vcf.gz > /store/lkemp/fix_dumb_bug/pipeline_wintermute_output/DHB4177_filtered_dbnsfp_cleaned.vcf.gz

bcftools view -O z --max-alleles 2 --exclude-types indels \
/store/lkemp/fix_dumb_bug/pipeline_wintermute_output/DHB4177_filtered_dbnsfp_vep.vcf.gz > /store/lkemp/fix_dumb_bug/pipeline_wintermute_output/DHB4177_filtered_dbnsfp_vep_cleaned.vcf.gz

bcftools view -O z --max-alleles 2 --exclude-types indels \
/store/lkemp/fix_dumb_bug/pipeline_wintermute_output/DHB4177_filtered_dbnsfp_vep_cadd_dbsnp.vcf > /store/lkemp/fix_dumb_bug/pipeline_wintermute_output/DHB4177_filtered_dbnsfp_vep_cadd_dbsnp_cleaned.vcf

```

Test loading into scout

```bash
scout -db demo-database load case /store/lkemp/fix_dumb_bug/DHB4177_config.yaml

# Run through pipeline on Wintermute
DHB4177_filtered_cleaned.vcf #
DHB4177_filtered_dbnsfp_cleaned.vcf.gz # works
DHB4177_filtered_dbnsfp_vep_cleaned.vcf.gz # works
DHB4177_filtered_dbnsfp_vep_cadd_dbsnp_cleaned.vcf # works
```

**Error seems to introduced at the dbNSFP annotation step - only when the pipeline is run on the ESR cluster**

Possible contributing factors:

- Differences in the dbNSFP database between Wintermute and the ESR cluster
- Differences between deploying to the HPC or not

1. Run the pipeline again on the ESR cluster (but not deployed to the cluster) to see if the same error is introduced at the dbNSFP annotation step
2. Check for any differences between the dbNSFP databases between Wintermute and the ESR cluster
3. Re-move the dbNSFP database from Wintermute to the ESR cluster - see if a fresh move will fix the error
4. Check for differences between `/store/lkemp/fix_dumb_bug/pipeline_cluster_output/DHB4177_filtered_dbnsfp.vcf.gz` and  `/store/lkemp/fix_dumb_bug/pipeline_wintermute_output/DHB4177_filtered_dbnsfp.vcf.gz`

## 2. Check for any differences between the dbNSFP databases between Wintermute and the ESR cluster

- The dbNSFP database on Wintermute has a bunch of associated tabix files such as `dbNSFPv4.0a.hg19.custombuild.gz.tabixIndex_1.txt` where the dbNSFP database on the ESR cluster doesn't
- The size of the two databases are the same (in megabytes)
- Both gzip compressed files

```bash
(base) orac$ file dbNSFPv4.0a.hg19.custombuild.gz
dbNSFPv4.0a.hg19.custombuild.gz: gzip compressed data, extra field
```

## 4. Re-move the dbNSFP database from Wintermute to the ESR cluster - see if a fresh move will fix the error

Move fresh dbNSFP database to the cluster (inlcuding all the associated tabix files)

```bash
# Run code from Wintermute
scp -r /store/lkemp/publicData/dbNSFP/GRCh37/dbNSFPv4.0a.hg19.custombuild.gz* orac:/NGS/scratch/KSCBIOM/HumanGenomics/publicData/dbNSFP_fresh/GRCh37/
```

Setup to run pipeline again with fresh dbNSFP database

```bash
cd /NGS/scratch/KSCBIOM/HumanGenomics/fix_dumb_bug/fresh_dbnsfp/

# Get pipeline inputs
cp -r ../deploy_to_cluster/human_genomics_pipeline/ .
cp -r ../deploy_to_cluster/pedigrees/ .

# Clone same pipeline
git clone https://github.com/ESR-NZ/vcf_annotation_pipeline.git
cd vcf_annotation_pipeline
git checkout 75774b4ffef3258d0c7a5cff13072971e8553104

# Copy configuration files
cp ../../deploy_to_cluster/vcf_annotation_pipeline/config/* ./config/
cp ../../deploy_to_cluster/vcf_annotation_pipeline/workflow/*run*.sh ./workflow/
```

Manually set the dbNSFP database in the config file to use the fresh database

Run the pipeline

```bash
cd /NGS/scratch/KSCBIOM/HumanGenomics/fix_dumb_bug/fresh_dbnsfp/vcf_annotation_pipeline/workflow/
bash dryrun_hpc.sh
bash run_hpc.sh
```

Transfer files to Wintermute to test loading into scout (catch intermediate files before they are removed by the pipeline)

```bash
# Run this code from wintermute
scp -r orac:/NGS/scratch/KSCBIOM/HumanGenomics/fix_dumb_bug/fresh_dbnsfp/vcf_annotation_pipeline/results/filtered/DHB4177_filtered.vcf /store/lkemp/fix_dumb_bug/fresh_dbnsfp/
scp -r orac:/NGS/scratch/KSCBIOM/HumanGenomics/fix_dumb_bug/fresh_dbnsfp/vcf_annotation_pipeline/results/annotated/DHB4177_filtered_dbnsfp.vcf.gz /store/lkemp/fix_dumb_bug/fresh_dbnsfp/
scp -r orac:/NGS/scratch/KSCBIOM/HumanGenomics/fix_dumb_bug/fresh_dbnsfp/vcf_annotation_pipeline/results/annotated/DHB4177_filtered_dbnsfp_vep.vcf.gz /store/lkemp/fix_dumb_bug/fresh_dbnsfp/
scp -r orac:/NGS/scratch/KSCBIOM/HumanGenomics/fix_dumb_bug/fresh_dbnsfp/vcf_annotation_pipeline/results/annotated/DHB4177_filtered_dbnsfp_vep_cadd_dbsnp.vcf /store/lkemp/fix_dumb_bug/fresh_dbnsfp/
```

Filter multiallelic sites

```bash
conda activate bcftools

bcftools view -O z --max-alleles 2 --exclude-types indels \
/store/lkemp/fix_dumb_bug/fresh_dbnsfp/DHB4177_filtered.vcf > /store/lkemp/fix_dumb_bug/fresh_dbnsfp/DHB4177_filtered_cleaned.vcf

bcftools view -O z --max-alleles 2 --exclude-types indels \
/store/lkemp/fix_dumb_bug/fresh_dbnsfp/DHB4177_filtered_dbnsfp.vcf.gz > /store/lkemp/fix_dumb_bug/fresh_dbnsfp/DHB4177_filtered_dbnsfp_cleaned.vcf.gz

bcftools view -O z --max-alleles 2 --exclude-types indels \
/store/lkemp/fix_dumb_bug/fresh_dbnsfp/DHB4177_filtered_dbnsfp_vep.vcf.gz > /store/lkemp/fix_dumb_bug/fresh_dbnsfp/DHB4177_filtered_dbnsfp_vep_cleaned.vcf.gz

bcftools view -O z --max-alleles 2 --exclude-types indels \
/store/lkemp/fix_dumb_bug/fresh_dbnsfp/DHB4177_filtered_dbnsfp_vep_cadd_dbsnp.vcf > /store/lkemp/fix_dumb_bug/fresh_dbnsfp/DHB4177_filtered_dbnsfp_vep_cadd_dbsnp_cleaned.vcf
```

Test loading into scout

```bash
scout -db demo-database load case /store/lkemp/fix_dumb_bug/DHB4177_config.yaml

# Run through pipeline on the ESR cluster
DHB4177_filtered_cleaned.vcf # works
DHB4177_filtered_dbnsfp_cleaned.vcf.gz # 
DHB4177_filtered_dbnsfp_vep_cleaned.vcf.gz # 
DHB4177_filtered_dbnsfp_vep_cadd_dbsnp_cleaned.vcf # 
```
