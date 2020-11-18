# Test mosdepth and seqcover

Created: 2020/11/19 09:42:29
Last modified: 2020/11/19 10:25:12

- **Aim:** Try/test [mosdepth](https://github.com/brentp/mosdepth) and [seqcover](https://github.com/brentp/seqcover) as a method for "viewing and evaluating depth-of-coverage" over many samples and for each gene
- **Prerequisite software:** [conda](https://singularity.lbl.gov/)
- **OS:** Ubuntu 16.04 (Wintermute - research server)

## Table of contents

- [Test mosdepth and seqcover](#test-mosdepth-and-seqcover)
  - [Table of contents](#table-of-contents)
  - [Try mosdepth and seqcover](#try-mosdepth-and-seqcover)
  - [Notes](#notes)

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

Try on existing data

```bash
mkdir -p samples/

# Generate report
for b in *.bam; do
  n=$(basename $b .bam)
  mosdepth -x -t 4 samples/$n $b
done

seqcover report --genes PIGA,KCNQ2,ARX,DNM1,SLC25A22,CDKL5,GABRA1,CAD,MDH2,SCN1B,CNPY3,CPLX1,NEB,HNRNPA1,CCDC39,AIFM1,CHCHD10 \
		 --background seqcover/seqcover_p5.d4 \
		 --fasta $fasta samples/*.bed.gz \
		 -r my_genes_report.html

# Generate a background level
seqcover generate-background --percentile 5 -f $fasta -o seqcover/ d4s/HG00*.d4

# Generate a transcript file
seqcover save-transcripts --genes PIGA,KCNQ2,ARX,DNM1,SLC25A22,CDKL5,GABRA1,CAD,MDH2,SCN1B,CNPY3,CPLX1,NEB,HNRNPA1,CCDC39,AIFM1,CHCHD10 \
		 --output-path transcripts.json \
		 --hg19

seqcover report --genes PIGA,KCNQ2,ARX,DNM1,SLC25A22,CDKL5,GABRA1,CAD,MDH2,SCN1B,CNPY3,CPLX1,NEB,HNRNPA1,CCDC39,AIFM1,CHCHD10 \
		 --background seqcover/seqcover_p5.d4 \
		 --fasta $fasta samples/*.bed.gz \
		 -r my_genes_report.html \
		 --transcripts-file transcripts.json
```


## Notes

- It take's a little bit of time to compute, but they say you should be able to parallize it
- Will have to remember that you should only be running seqcover on exome data produced with the same capture kit (and prefereably from the same sequencing center) - shouldn't be a problem for our clinical exome sequencing/data

