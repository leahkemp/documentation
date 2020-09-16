# Test somalier as replacement for peddy

Created: 2020/08/31 11:36:14
Last modified: 2020/09/17 12:00:38

- **Aim:** Try/test [somalier](https://github.com/brentp/somalier) and see if the outputs can be integrated with [scout](https://github.com/Clinical-Genomics/scout)
- **Prerequisite software:** [singularity 2.5.2](https://singularity.lbl.gov/)
- **OS:** Ubuntu 16.04 (Wintermute - research server)

In our whole genome sequencing (WGS) workflow (fastq to vcf to preparing data to be ingested into scout) have been using [peddy](https://github.com/brentp/peddy) to evaluate/check the sex/familial relationships for the joint vcf's of our cohort samples. The files that peddy outputs are used and presented in [scout](https://github.com/Clinical-Genomics/scout) to show whether the genomic data supports the familial relationships and sex's defined in the pedigree file. We have had some inconsistency in the results of these checks with peddy depending on how the data is subset (eg. exome capture regions only and whether non-pass variants removed from the vcf). [somalier](https://github.com/brentp/somalier) has been released which ["is a more scalable, faster, replacement for peddy that uses some of the same methods as peddy along with some new ones"](https://github.com/brentp/peddy#fast-pedigreevcf-qc). It also allows you to do these checks based on bam's instead of vcf's. So we want to try it out and see if it's outputs can still be integrated into scout, and therefore included in our pipeline/workflow!

## Table of contents

- [Test somalier as replacement for peddy](#test-somalier-as-replacement-for-peddy)
  - [Table of contents](#table-of-contents)
  - [Try somalier](#try-somalier)
    - [Based on vcf](#based-on-vcf)
    - [Based on bam](#based-on-bam)
  - [Notes](#notes)

## Try somalier

Get data to play with from manual pipeline runs

```bash
cd /store/lkemp/test_somalier/
cp /store/lkemp/manual_pipeline_run/WGS_b37/DHB4205_GRCh37_WGS_trio.combined.gpu.genotyped.vqsr.withPosteriors.denovo_ann.dbnsfp.vep.cadd.models.sorted.dbsnp.onlyproband.dualalleles.rankscore.exomes.clean.vcf.gz* .
cp /store/lkemp/manual_pipeline_run/WGS_b37/DHB4205_pedigree.ped .
```

A docker image for somalier is hosted [here](https://hub.docker.com/r/brentp/somalier)

```bash
singularity pull docker://brentp/somalier:latest
```

Get sites files from [here](https://github.com/brentp/somalier/releases/tag/v0.2.12)

```bash
wget https://github.com/brentp/somalier/files/3412453/sites.hg19.vcf.gz
wget https://github.com/brentp/somalier/files/3412454/sites.hg38.nochr.vcf.gz
wget https://github.com/brentp/somalier/files/3412455/sites.GRCh37.vcf.gz
wget https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz
```

Get ancestry files

```bash
wget https://raw.githubusercontent.com/brentp/somalier/master/scripts/ancestry-labels-1kg.tsv
wget https://zenodo.org/record/3479773/files/1kg.somalier.tar.gz?download=1
tar -xf 1kg.somalier.tar.gz?download=1
```

### Based on vcf

Extract sites from vcf

```bash
BIND_PATHS="/store/"
CONTAINER="/store/lkemp/test_somalier/somalier-latest.simg"

singularity exec \
-B ${BIND_PATHS} \
${CONTAINER} \
somalier extract \
-d vcf_extracted/ \
--sites sites.GRCh37.vcf.gz \
-f /store/lkemp/publicData/b37/human_g1k_v37_decoy.fasta \
DHB4205_GRCh37_WGS_trio.combined.gpu.genotyped.vqsr.withPosteriors.denovo_ann.dbnsfp.vep.cadd.models.sorted.dbsnp.onlyproband.dualalleles.rankscore.exomes.clean.vcf.gz
```

Calculate relatedness on the extracted data

```bash
singularity exec \
-B ${BIND_PATHS} \
${CONTAINER} \
somalier relate \
--ped DHB4205_pedigree.ped \
vcf_extracted/*.somalier \
--output-prefix=vcf_extracted
```

Calculate ancestry

```bash
singularity exec \
-B ${BIND_PATHS} \
${CONTAINER} \
somalier ancestry \
--labels ancestry-labels-1kg.tsv \
1kg-somalier/*.somalier ++ vcf_extracted/*.somalier
```

### Based on bam

Extract sites from the three bams in the cohort

```bash
for f in /store/mbenton/WGS_b37/bams/DHB420*.bam; do
    singularity exec \
    -B ${BIND_PATHS} \
    ${CONTAINER} \
    somalier extract \
    -d bam_extracted/ \
    --sites sites.GRCh37.vcf.gz \
    -f /store/lkemp/publicData/b37/human_g1k_v37_decoy.fasta \
    $f
done
```

Calculate relatedness on the extracted data

```bash
singularity exec \
-B ${BIND_PATHS} \
${CONTAINER} \
somalier relate \
--ped DHB4205_pedigree.ped \
bam_extracted/*.somalier \
--output-prefix=bam_extracted
```

Calculate ancestry

```bash
singularity exec \
-B ${BIND_PATHS} \
${CONTAINER} \
somalier ancestry \
--labels ancestry-labels-1kg.tsv \
1kg-somalier/*.somalier ++ bam_extracted/*.somalier
```

## Notes

- It looks like it will be fairly straightforward to include ethnic representation beyond what is included in the [1000 genomes project](https://www.internationalgenome.org/) for the ancestry estimates. It looks as though all we would need to do is run somalier on samples from ethnicities we want to include in the dataset and add them to the ancestry-labels-1kg.tsv file passed to `--labels`
