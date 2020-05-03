# Testing vg (variation graphs)

Created: 2020-04-30 11:38:42
Last modified: 2020/04/30 14:55:10

- **Aim:** Test the utility of [vg (variation graphs)](https://github.com/vgteam/vg) for use in creating human pan genomes/graph genomes
- **Prerequisite software:**
- **OS:** Ubuntu 16.04 (Wintermute - research server)

## Table of contents

- [Testing vg (variation graphs)](#testing-vg-variation-graphs)
  - [Table of contents](#table-of-contents)
  - [Clone repo](#clone-repo)
  - [Download vg](#download-vg)
  - [Trials](#trials)
    - [Construct graph genome](#construct-graph-genome)

## Clone repo

```bash
git clone git@github.com:vgteam/vg.git
```

## Download vg

Download latest release builds for Linux

```bash
wget https://github.com/vgteam/vg/releases/download/v1.23.0/vg
```

Make executable

```bash
chmod +x vg
```

## Trials

### Construct graph genome

```bash
vg construct \
-r /store/lkemp/publicData/referenceGenome/gatkBundle/GRCh37/ucsc.hg19.fasta \
-v /home/lkemp/pan_genome_project/clinvar_20200415.vcf.gz \
>test.vg
```

Notes:

- vcf file *can't* include any gaps or IUPAC ambiguity codes (the vcf output of both the human_genomics_pipeline and vcf_annotation_pipeline appear to have these)


Get example data

```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive_2.0/2020/clinvar_20200415.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive_2.0/2020/clinvar_20200415.vcf.gz.tbi
```