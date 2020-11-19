# Test mosdepth and seqcover

Created: 2020/11/19 09:42:29
Last modified: 2020/11/19 17:28:02

- **Aim:** Try/test [mosdepth](https://github.com/brentp/mosdepth) and [seqcover](https://github.com/brentp/seqcover) as a method for "viewing and evaluating depth-of-coverage" over many samples and for each gene
- **Prerequisite software:** [conda 4.9.0](https://singularity.lbl.gov/)
- **OS:** Ubuntu 16.04 (Wintermute - research server)

## Table of contents

- [Test mosdepth and seqcover](#test-mosdepth-and-seqcover)
  - [Table of contents](#table-of-contents)
  - [Try mosdepth and seqcover](#try-mosdepth-and-seqcover)
  - [Notes/findings](#notesfindings)

## Try mosdepth and seqcover

```bash
cd /store/lkemp/GA_clinical_genomics/run_1/bams/
```

Create conda env

```bash
conda create -n mosdepth_seqcover python=3.7.6
conda activate mosdepth_seqcover
```

Install tools

```bash
# mosdepth
conda install mosdepth=0.3.1

# seqcover
wget https://github.com/brentp/seqcover/releases/download/v0.0.2/seqcover
chmod +x ./seqcover
```

Get GRCh37 reference genome (b37 doesn't seem to work for seqcover)

```bash
wget https://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.fai
wget https://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
```

Prepare reference genome (seqcover doesn't seem to like it being gzipped)

```bash
gunzip human_g1k_v37.fasta.gz
bgzip human_g1k_v37.fasta
```

Try on existing data

```bash
mkdir -p samples/

# Mosdepth
for b in *.bam; do
  n=$(basename $b .bam)
  mosdepth -x -t 4 samples/$n $b
done

# Generate report
./seqcover report \
--genes PIGA,KCNQ2,ARX,DNM1,SLC25A22,CDKL5,GABRA1,CAD,MDH2,SCN1B,CNPY3,CPLX1,NEB,HNRNPA1,CCDC39,AIFM1,CHCHD10 \
--fasta human_g1k_v37.fasta.gz \
./samples/*.bed.gz \
-r my_genes_report.html
```

Can't get it working, error:

```bash
[seqcover] read 3 sample coverage files
tables.nim(262)          []
Error: unhandled exception: key not found: exons [KeyError]
```

Looks related to [this issue](https://github.com/brentp/seqcover/issues/22)

Will try with the new binary provided

```bash
rm seqcover
wget https://github.com/brentp/seqcover/files/5416692/seqcover.gz
gunzip seqcover.gz
chmod +x ./seqcover
```

Try seqcover again

```bash
./seqcover report \
--genes PIGA,KCNQ2,ARX,DNM1,SLC25A22,CDKL5,GABRA1,CAD,MDH2,SCN1B,CNPY3,CPLX1,NEB,HNRNPA1,CCDC39,AIFM1,CHCHD10 \
--fasta human_g1k_v37.fasta.gz \
./samples/*.bed.gz \
--hg19 \
-r my_genes_report.html
```

Worked!

In [this issue](https://github.com/brentp/seqcover/issues/22), they indicate the new binary I used doesn't work with the `--hg19` flag. I want to use this flag so it doesn't default to hg38 (we analysed with b37). Will try again with this flag

```bash
./seqcover report \
--genes PIGA,KCNQ2,ARX,DNM1,SLC25A22,CDKL5,GABRA1,CAD,MDH2,SCN1B,CNPY3,CPLX1,NEB,HNRNPA1,CCDC39,AIFM1,CHCHD10 \
--fasta human_g1k_v37.fasta.gz \
./samples/*.bed.gz \
--hg19 \
-r my_genes_report.html
```

Also doesn't work for me, but they indicate they are actively working on a fix for this!

I want to try it with the same reference genome we used to analyse the data (b37)

```bash
./seqcover report \
--genes PIGA,KCNQ2,ARX,DNM1,SLC25A22,CDKL5,GABRA1,CAD,MDH2,SCN1B,CNPY3,CPLX1,NEB,HNRNPA1,CCDC39,AIFM1,CHCHD10 \
--fasta /store/lkemp/publicData/b37/human_g1k_v37_decoy.fasta \
./samples/*.bed.gz \
-r my_genes_report.html
```

Worked great! (see the [output report](./test_mosdepth_and_seqcover/my_genes_report.html))

I also want to try seqcover with a transcript file

```bash
# Generate a transcript file
./seqcover save-transcripts \
--genes PIGA,KCNQ2,ARX,DNM1,SLC25A22,CDKL5,GABRA1,CAD,MDH2,SCN1B,CNPY3,CPLX1,NEB,HNRNPA1,CCDC39,AIFM1,CHCHD10 \
--output-path transcripts.json
```

Huh, error:

```bash
unknown program 'save-transcripts'
```

Doesn't seem like that's an available command despite [the README](https://github.com/brentp/seqcover/blob/master/README.md) referring to using this command

```bash
./seqcover --help
```

Output:

```bash
  generate-background:   generate background file(s) from a set of samples
  report             :   create an HTML report from a set of sample coverage files
```

Only `report` and `generate-background` seem to be available commands

Ah well, I'll try out the `generate-background` command

```bash
# Generate background
./seqcover generate-background \
--percentile 5 \
-f /store/lkemp/publicData/b37/human_g1k_v37_decoy.fasta \
-o seqcover/ \
./samples/*.bed.gz

# Re-create report
./seqcover report \
--genes PIGA,KCNQ2,ARX,DNM1,SLC25A22,CDKL5,GABRA1,CAD,MDH2,SCN1B,CNPY3,CPLX1,NEB,HNRNPA1,CCDC39,AIFM1,CHCHD10 \
--background seqcover/seqcover_p5.d4 \
--fasta /store/lkemp/publicData/b37/human_g1k_v37_decoy.fasta \
./samples/*.bed.gz \
-r my_genes_report.html
```

Error:

```bash
[seqcover] need at least 4 samples to generate a background
```

Cool cool, I'll try again on more samples

```bash
cd /store/lkemp/GA_clinical_genomics/run_3/bams/

mkdir -p samples/

# Mosdepth
for b in *.bam; do
  n=$(basename $b .bam)
  mosdepth -x -t 4 samples/$n $b
done

# Re-create report
../../run_1/bams/seqcover report \
--genes PIGA,KCNQ2,ARX,DNM1,SLC25A22,CDKL5,GABRA1,CAD,MDH2,SCN1B,CNPY3,CPLX1,NEB,HNRNPA1,CCDC39,AIFM1,CHCHD10 \
--fasta /store/lkemp/publicData/b37/human_g1k_v37_decoy.fasta \
./samples/*.bed.gz \
-r my_genes_report_2.html
```

Worked great! ([see the output report](./test_mosdepth_and_seqcover/my_genes_report_2.html))

## Notes/findings

- mosdepth take's a little bit of time to compute, but they say you should be able to parallize it
- seqcover is very speedy
- Will have to remember that you should only be running seqcover on exome data produced with the same capture kit (and preferably from the same sequencing centre) - shouldn't be a problem for our clinical exome sequencing/data
- Looks like you can choose the gene's you want to visualise, I'm not sure if it scales to all the genes in a genome since they say the output html report will be responsive up to 200 genes
- There are some early bugs in seqcover (new software!) that might limit how we want to use this software (eg. using the `--hg19 flag` and using the `save-transcripts` command). They look to be in active development though so I suspect it will very soon be ready for use :smile:

