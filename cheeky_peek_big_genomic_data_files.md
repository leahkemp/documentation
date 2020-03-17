# A cheeky peek at big genomic data files for human genomic pipelines

Created: 2020/03/18 09:47:43
Last modified: 2020/03/18 12:07:32

The reference human genomes and other databases files differ when downloaded from different sources (eg. chromosome labelling). Unfortunately mixing files sourced from different hosts within a pipeline can cause computer says no's. The typical advice is to [use data sources from the same place](https://gatkforums.broadinstitute.org/gatk/discussion/11359/input-files-reference-and-features-have-incompatible-contigs).

However, you may need to combine data files sourced from different hosts within a given genomic pipeline. For example I need a more recent build of dbSNP for GRCh37 which is available from [NCBI](https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/), but I also need to annotate my VCF files later in my pipeline with databases such as 1000G indels and 1000G snps available in the [GATK resource bundle](https://gatkforums.broadinstitute.org/gatk/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it).

Don't want to be stuck using only the GATK resource bundle? You can use this document as a reference to manually make your data files compatible when they sourced from different hosts (for example by editing the chromosome labelling).

On another note, I think it is helpful and general good practice to have a look at the data you're working with, especially when dealing with big data.

So here are the heads and tails of the reference human genomes and other databases (excluding the info section at the beginning of the files) used by ['human_genomics_pipeline'](https://github.com/ESR-NZ/human_genomics_pipeline) and [vcf_annotation_pipeline](https://github.com/leahkemp/vcf_annotation_pipeline.git).

*note. this document does not describe all reference genome/databases available*

## Table of contents

- [A cheeky peek at big genomic data files for human genomic pipelines](#a-cheeky-peek-at-big-genomic-data-files-for-human-genomic-pipelines)
  - [Table of contents](#table-of-contents)
  - [GRCh37/hg19](#grch37hg19)
    - [Reference human genome](#reference-human-genome)
      - [Sourced from NCBI](#sourced-from-ncbi)
      - [Sourced from the GATK resource bundle](#sourced-from-the-gatk-resource-bundle)
      - [Sourced from UCSC](#sourced-from-ucsc)
    - [dbSNP](#dbsnp)
      - [Sourced from NCBI](#sourced-from-ncbi-1)
      - [Sourced from the GATK resource bundle](#sourced-from-the-gatk-resource-bundle-1)
    - [Hapmap](#hapmap)
      - [Sourced from the GATK resource bundle](#sourced-from-the-gatk-resource-bundle-2)
    - [OMNI](#omni)
      - [Sourced from the GATK resource bundle](#sourced-from-the-gatk-resource-bundle-3)
    - [1000G indels](#1000g-indels)
      - [Sourced from the GATK resource bundle](#sourced-from-the-gatk-resource-bundle-4)
    - [1000G snps](#1000g-snps)
      - [Sourced from the GATK resource bundle](#sourced-from-the-gatk-resource-bundle-5)
    - [Mills](#mills)
      - [Sourced from the GATK resource bundle](#sourced-from-the-gatk-resource-bundle-6)
    - [Axiom exome plus](#axiom-exome-plus)
      - [Sourced from the GATK resource bundle](#sourced-from-the-gatk-resource-bundle-7)
    - [CADD](#cadd)
      - [Sourced from CADD](#sourced-from-cadd)
  - [GRCh38/hg38](#grch38hg38)
    - [Reference human genome](#reference-human-genome-1)
      - [Sourced from NCBI](#sourced-from-ncbi-2)
      - [Sourced from the GATK resource bundle](#sourced-from-the-gatk-resource-bundle-8)
      - [Sourced from UCSC](#sourced-from-ucsc-1)
    - [dbSNP](#dbsnp-1)
      - [Sourced from NCBI](#sourced-from-ncbi-3)
      - [Sourced from the GATK resource bundle](#sourced-from-the-gatk-resource-bundle-9)
    - [Hapmap](#hapmap-1)
      - [Sourced from the GATK resource bundle](#sourced-from-the-gatk-resource-bundle-10)
    - [OMNI](#omni-1)
      - [Sourced from the GATK resource bundle](#sourced-from-the-gatk-resource-bundle-11)
    - [1000G indels](#1000g-indels-1)
      - [Sourced from the GATK resource bundle](#sourced-from-the-gatk-resource-bundle-12)
    - [1000G snps](#1000g-snps-1)
      - [Sourced from the GATK resource bundle](#sourced-from-the-gatk-resource-bundle-13)
    - [Mills](#mills-1)
      - [Sourced from the GATK resource bundle](#sourced-from-the-gatk-resource-bundle-14)
    - [Axiom exome plus](#axiom-exome-plus-1)
      - [Sourced from the GATK resource bundle](#sourced-from-the-gatk-resource-bundle-15)
    - [CADD](#cadd-1)
      - [Sourced from CADD](#sourced-from-cadd-1)

## GRCh37/hg19

### Reference human genome

#### Sourced from [NCBI](https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/)

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.fna.gz
```

```bash
NC_000001.10 Homo sapiens chromosome 1, GRCh37.p13 Primary Assembly
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
.
.
.
>NC_012920.1 Homo sapiens mitochondrion, complete genome
GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTATGCACGC
GATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGATTCCTGCCTCATCCTATTATTTA
TCGCACCTACGTTCAATATTACAGGCGAACATACTTACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATA
ACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTC
TGGCCACAGCACTTAAACACATCTCTGCCAAACCCCAAAAACAAAGAACCCTAACACCAGCCTAACCAGATTTCAAATTT
TATCTTTTGGCGGTATGCACTTTTAACAGTCACCCCCCAACTAACACATTATTTTCCCCTCCCACTCCCATACTACTAAT
CTCATCAATACAACCCCCGCCCATCCTACCCAGCACACACACACCGCTGCTAACCCCATACCCCGAACCAACCAAACCCC
AAAGACACCCCCCACAGTTTATGTAGCTTACCTCCTCAAAGCAATACACTGAAAATGTTTAGACGGGCTCACATCACCCC
ATAAACAAATAGGTTTGGTCCTAGCCTTTCTATTAGCTCTTAGTAAGATTACACATGCAAGCATCCCCGTTCCAGTGAGT
TCACCCTCTAAATCACCACGATCAAAAGGAACAAGCATCAAGCACGCAGCAATGCAGCTCAAAACGCTTAGCCTAGCCAC
ACCCCCACGGGAAACAGCAGTGATTAACCTTTAGCAATAAACGAAAGTTTAACTAAGCTATACTAACCCCAGGGTTGGTC
.
.
.
```

#### Sourced from the [GATK resource bundle](https://gatkforums.broadinstitute.org/gatk/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it)

```bash
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.fasta.gz
```

```bash
>chrM
GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCAT
TTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTG
GAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGATTCCTGCCTCATT
CTATTATTTATCGCACCTACGTTCAATATTACAGGCGAACATACCTACTA
AAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAAT
GTCTGCACAGCCGCTTTCCACACAGACATCATAACAAAAAATTTCCACCA
AACCCCCCCCTCCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGC
CAAACCCCAAAAACAAAGAACCCTAACACCAGCCTAACCAGATTTCAAAT
TTTATCTTTAGGCGGTATGCACTTTTAACAGTCACCCCCCAACTAACACA
TTATTTTCCCCTCCCACTCCCATACTACTAATCTCATCAATACAACCCCC
GCCCATCCTACCCAGCACACACACACCGCTGCTAACCCCATACCCCGAAC
CAACCAAACCCCAAAGACACCCCCCACAGTTTATGTAGCTTACCTCCTCA
AAGCAATACACTGAAAATGTTTAGACGGGCTCACATCACCCCATAAACAA
ATAGGTTTGGTCCTAGCCTTTCTATTAGCTCTTAGTAAGATTACACATGC
AAGCATCCCCGTTCCAGTGAGTTCACCCTCTAAATCACCACGATCAAAAG
GGACAAGCATCAAGCACGCAGCAATGCAGCTCAAAACGCTTAGCCTAGCC
ACACCCCCACGGGAAACAGCAGTGATTAACCTTTAGCAATAAACGAAAGT
TTAACTAAGCTATACTAACCCCAGGGTTGGTCAATTTCGTGCCAGCCACC
GCGGTCACACGATTAACCCAAGTCAATAGAAGCCGGCGTAAAGAGTGTTT
TAGATCACCCCCTCCCCAATAAAGCTAAAACTCACCTGAGTTGTAAAAAA
.
.
.
chrUn_gl000249
GATCACCAAGGCTGGGAAATATTAGGGAAGGGGGGTTGCCAGACTCCGGA
CCTGCAGAGTTAAGCCTTCCTCCCCCACCCCCACCACCCCATGTCATGTG
GCTGGTGGCAATTCCCTCCAGGACAAAAGGCCCGATTTAATCCAGCCCAC
CATCACCACTGTCGCCACTGGGACACAATGCGGCAGGTTTGTGGCCAACG
TCTGCCGATCTGGTTTCGTGTAACATCCCTGCCAGCCTGCCCGGGCCAGC
AGACAAAGGCCTCTTTGTTGCAAATATGTTTTTTAAATCCCTGAAGATAT
TAGCAGTGCGGGTGAACTCACACTGTGAAACAGTTCAGAAATTGTTTAAG
GACATGTTTCAAACTGGGGGTGATCATTTAAATGGAATCTGCCCTCCTGC
TTTCTTATCGAGAGCAAGATTCCTCAGAGCCAGCTTGGGCCCTGGACCTG
GGCAGGGAAGTTTCCGAGGCCAAATAACCCTAACACTCATCGTAACAACA
ATGCCAGCAGCCATTTATGGAGTCTGCAAGGCCAGACGCCGGCCCACTTG
TGCATTTGTGCATGTGACCGCCAAACCTCGGCActtagagagattcagat
gagcagtcccattctacagatggaaaatctggggtctagggaagtgTTGG
CGCACGTCACAGCTGGCCTTGAACAAGGTTTGTTGGATCCTCCACACAGC
CCCTGTGCCCTCTTCATGCTGTTGATGGCCACAGGCACGGAACACATCCT
CACCTCTGAACTCTGTCCCATCAGGCACCCAGGAGGGCATCCGACGGGCA
CCTGGTCCCCTCACATCCTGTTCCTTCAAGTTCACACAAGAAGGCCTCCT
CCAACACATCTCTCCACCTACTTctctgttacctgggctgagtgcagtgg
tcacaatcatggctcactgcagcctcaacgtcccaagctccaccaatcct
.
.
.
```

#### Sourced from [UCSC](https://hgdownload.soe.ucsc.edu/downloads.html)

```bash
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
```

```bash
>chr1
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
.
.
.
>chr18_gl000207_random
CATACATTTAATATACCCTCACCATACAGAATGTTCTTTCCCTATTACAT
AAGGAGTATATGTATTAAGCACTAAATCTTTGGAATAATAAAAGACTATA
TTCATATTTGGTAACTTATTTAATCAGAACAGGTTTACAACATAAATAAA
TAGATTGAACTTACTTTGTATAAAAATTGTATTATCTAAAACTTCACAGA
GAAAATAGGAAAACTCCACGTTAAGTGTTCAGGTCTGCTAAGTCATATTT
TTAAATCTAAGTGTAGGCATGAAAATGACTTCATTGTCATGGTAGTAATA
GCTGAACTGGTCAATAAGTAAGTCTTGTTTTAGAAGTGGTTGTTTCTGGT
AACCACTATAACCCCACGTTTTTGGAGTTATATGTTGGCACTGATACTGG
CCATAGAATTCCCTATGGTATTGATGGTTGATACAGAATCTGGTATATCT
GGTGTAGGAATGAGGATAGTACATGTTAGGATTTAATTTCTAGGTATACT
AACCACAGTTGCCCCAGGAATAGGTCATTAGAAAGTTGTCATTTTTGTCA
AGTGAGAGTAGTAATCTGACTTAGTCAATTTACTAGAAAAAAATGATAGT
GTTGCATTTATGAATAAGCACTGAATGACTTTGAATAACATGCTAATAGA
AATTGTGATGCTTGTAGATTCCAGTGATAGGGTACTAAAAGGTTCAGTAA
AGTTTTGTCTGTTGTAGGAGTACATAGAGTAGAAATGTGGTCTCTATTGT
GAAATTTCAGTAGATGTCATGGGAGTAGCACATTTTGCTGATGTGGCATG
GTTGCTGTGGCTATGAAGGGACTCATTTGAAGGGGCCATGGAAGAGCTGG
CAGTATTATGATAGCAGATTATGTAGAAACTGTAATTGCATTTTCTATTA
ACAGGTAATGATGCTCATGGAGATTGTAAACTGAGCTGACATTGCTAGGG
.
.
.
```

### dbSNP

#### Sourced from [NCBI](https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/)

```bash
# Build 153
wget ftp://ftp.ncbi.nlm.nih.gov:21/snp/latest_release/VCF/GCF_000001405.25.gz
wget ftp://ftp.ncbi.nlm.nih.gov:21/snp/latest_release/VCF/GCF_000001405.25.gz.tbi
```

```bash
#CHROM  POS             ID              REF     ALT     QUAL    FILTER  INFO
chrM    64              rs3883917       C       T       .       .       ASP;CAF=[0.9373,0.06268];COMMON=1;OTHERKG;R5;RS=3883917;RSPOS=64;SAO=0;SSR=0;VC=SNV;VP=0x050000020005000002000100;WGT=1;dbSNPBuildID=108
chrM    72              rs370271105     T       C       .       .       ASP;CAF=[0.9878,0.01123];COMMON=1;HD;OTHERKG;R5;RS=370271105;RSPOS=72;SAO=0;SSR=0;VC=SNV;VP=0x050000020005000402000100;WGT=1;dbSNPBuildID=138
chrM    93              rs369034419     A       G       .       .       ASP;CAF=[0.9663,0.03368];COMMON=1;HD;OTHERKG;R5;RS=369034419;RSPOS=93;SAO=0;SSR=0;VC=SNV;VP=0x050000020005000402000100;WGT=1;dbSNPBuildID=138
chrM    97              rs147830800     G       A       .       .       ASP;CAF=[0.9972,0.002806];COMMON=1;OTHERKG;R5;RS=147830800;RSPOS=97;SAO=0;SSR=0;VC=SNV;VP=0x050000020005000002000100;WGT=1;dbSNPBuildID=134
chrM    103             rs369070397     G       A       .       .       ASP;CAF=[0.9906,0.009355];COMMON=1;OTHERKG;R5;RS=369070397;RSPOS=103;SAO=0;SSR=0;VC=SNV;VP=0x050000020005000002000100;WGT=1;dbSNPBuildID=138
chrM    125             rs144402189     T       C       .       .       ASP;CAF=[0.9991,0.0009355];COMMON=0;HD;OTHERKG;R5;RS=144402189;RSPOS=125;SAO=0;SSR=0;VC=SNV;VP=0x050000020005000402000100;WGT=1;dbSNPBuildID=134
chrM    127             rs139684161     T       C       .       .       ASP;CAF=[0.9991,0.0009355];COMMON=0;OTHERKG;R5;RS=139684161;RSPOS=127;SAO=0;SSR=0;VC=SNV;VP=0x050000020005000002000100;WGT=1;dbSNPBuildID=134
chrM    143             rs375589100     G       A       .       .       ASP;CAF=[0.9673,0.03274];COMMON=1;OTHERKG;R5;RS=375589100;RSPOS=143;SAO=0;SSR=0;VC=SNV;VP=0x050000020005000002000100;WGT=1;dbSNPBuildID=138
.
.
.
chrY    59362755        rs372555462     AG      A       .       .       ASP;CFL;OTHERKG;RS=372555462;RSPOS=59362756;SAO=0;SSR=0;VC=DIV;VP=0x05000000000d000002000200;WGT=1;dbSNPBuildID=138
chrY    59362881        rs368008987     T       G       .       .       ASP;CFL;OTHERKG;RS=368008987;RSPOS=59362881;SAO=0;SSR=0;VC=SNV;VP=0x05000000000d000002000100;WGT=1;dbSNPBuildID=138
chrY    59362903        rs374874946     GGT     G       .       .       ASP;CFL;OTHERKG;RS=374874946;RSPOS=59362904;SAO=0;SSR=0;VC=DIV;VP=0x05000000000d000002000200;WGT=1;dbSNPBuildID=138
chrY    59362926        rs369293409     AG      A       .       .       ASP;CFL;OTHERKG;RS=369293409;RSPOS=59362927;SAO=0;SSR=0;VC=DIV;VP=0x05000000000d000002000200;WGT=1;dbSNPBuildID=138
chrY    59362952        rs372923705     AG      A       .       .       ASP;CFL;OTHERKG;RS=372923705;RSPOS=59362953;SAO=0;SSR=0;VC=DIV;VP=0x05000000000d000002000200;WGT=1;dbSNPBuildID=138
chrY    59363034        rs374134835     TAG     T,TG    .       .       CFL;NOC;OTHERKG;RS=374134835;RSPOS=59363035;SAO=0;SSR=0;VC=DIV;VP=0x050000000009000002000210;WGT=1;dbSNPBuildID=138
chrY    59363052        rs367673729     TAG     T,TG    .       .       CFL;NOC;OTHERKG;RS=367673729;RSPOS=59363053;SAO=0;SSR=0;VC=DIV;VP=0x050000000009000002000210;WGT=1;dbSNPBuildID=138
chrY    59363407        rs376709892     AG      A       .       .       ASP;CFL;OTHERKG;RS=376709892;RSPOS=59363408;SAO=0;SSR=0;VC=DIV;VP=0x05000000000d000002000200;WGT=1;dbSNPBuildID=138
chrY    59363442        rs371324725     AG      A       .       .       ASP;CFL;OTHERKG;RS=371324725;RSPOS=59363443;SAO=0;SSR=0;VC=DIV;VP=0x05000000000d000002000200;WGT=1;dbSNPBuildID=138
```

#### Sourced from the [GATK resource bundle](https://gatkforums.broadinstitute.org/gatk/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it)

```bash
# Build 138
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/dbsnp_138.hg19.vcf.gz
```

```bash
#CHROM          POS     ID              REF     ALT                                           QUAL    FILTER  INFO
NC_000001.10    10019   rs775809821     TA      T                                             .       .       RS=775809821;dbSNPBuildID=144;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=INDEL
NC_000001.10    10039   rs978760828     A       C                                             .       .       RS=978760828;dbSNPBuildID=150;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV
NC_000001.10    10043   rs1008829651    T       A                                             .       .       RS=1008829651;dbSNPBuildID=150;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV
NC_000001.10    10051   rs1052373574    A       G                                             .       .       RS=1052373574;dbSNPBuildID=150;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV
NC_000001.10    10051   rs1326880612    A       AC                                            .       .       RS=1326880612;dbSNPBuildID=151;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=INDEL
NC_000001.10    10055   rs768019142     T       TA                                            .       .       RS=768019142;dbSNPBuildID=144;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=INDEL
NC_000001.10    10055   rs892501864     T       A                                             .       .       RS=892501864;dbSNPBuildID=150;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV
NC_000001.10    10063   rs1010989343    A       C                                             .       .       RS=1010989343;dbSNPBuildID=150;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV
NC_000001.10    10067   rs1489251879    T       TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC    .       .       RS=1489251879;dbSNPBuildID=151;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=INDEL
NC_000001.10    10077   rs1022805358    C       G                                             .       .       RS=1022805358;dbSNPBuildID=150;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV
NC_000001.10    10108   rs62651026      C       T                                             .       .       RS=62651026;dbSNPBuildID=129;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV
NC_000001.10    10108   rs1322538365    C       CT                                            .       .       RS=1322538365;dbSNPBuildID=151;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=INS
.
.
.
NW_004775435.1  307206  rs1322535525    GCT     G                                             .       .       RS=1322535525;dbSNPBuildID=151;SSR=0;GENEINFO=LOC101928689:101928689;VC=INDEL;GNO;FREQ=TOPMED:.,7.964e-06|Vietnamese:.,0.00463
NW_004775435.1  307213  rs1555965343    G       A                                             .       .       RS=1555965343;dbSNPBuildID=152;SSR=0;GENEINFO=LOC101928689:101928689;VC=SNV;GNO;FREQ=GnomAD:1,3.187e-05
NW_004775435.1  307217  rs1292345196    C       T                                             .       .       RS=1292345196;dbSNPBuildID=151;SSR=0;GENEINFO=LOC101928689:101928689;VC=SNV;GNO;FREQ=TOPMED:1,7.964e-06
NW_004775435.1  307223  rs577356758     AG      A                                             .       .       RS=577356758;dbSNPBuildID=142;SSR=0;GENEINFO=LOC101928689:101928689;VC=INDEL;GNO;FREQ=1000Genomes:.,0.001198|GnomAD:.,0.000255|TOPMED:.,0.0009158
NW_004775435.1  307229  rs565797977     G       A                                             .       .       RS=565797977;dbSNPBuildID=142;SSR=0;GENEINFO=LOC101928689:101928689;VC=SNV;GNO;FREQ=ALSPAC:0.9995,0.0005189|GnomAD:0.9999,9.564e-05|TOPMED:0.9998,0.000215|TWINSUK:1,0
NW_004775435.1  307231  rs572672074     C       G,T                                           .       .       RS=572672074;dbSNPBuildID=142;SSR=0;GENEINFO=LOC101928689:101928689;VC=SNV;GNO;FREQ=1000Genomes:0.9998,.,0.0001997|GnomAD:0.9999,.,6.375e-05|TOPMED:0.9998,0.000223,2.389e-05
NW_004775435.1  307232  rs961157725     G       A                                             .       .       RS=961157725;dbSNPBuildID=150;SSR=0;GENEINFO=LOC101928689:101928689;VC=SNV;GNO;FREQ=GnomAD:1,3.187e-05|TOPMED:0.9999,5.575e-05
NW_004775435.1  307237  rs1292431636    G       T                                             .       .       RS=1292431636;dbSNPBuildID=151;SSR=0;GENEINFO=LOC101928689:101928689;VC=SNV;GNO;FREQ=TOPMED:1,1.593e-05
NW_004775435.1  307246  rs1458020546    G       A                                             .       .       RS=1458020546;dbSNPBuildID=151;SSR=0;GENEINFO=LOC101928689:101928689;VC=SNV;GNO;FREQ=TOPMED:1,2.389e-05
NW_004775435.1  307248  rs1351236822    G       A                                             .       .       RS=1351236822;dbSNPBuildID=151;SSR=0;GENEINFO=LOC101928689:101928689;VC=SNV;GNO;FREQ=TOPMED:1,1.593e-05
```

### Hapmap

#### Sourced from the [GATK resource bundle](https://gatkforums.broadinstitute.org/gatk/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it)

```bash
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/hapmap_3.3.hg19.sites.vcf.gz
```

```bash
#CHROM  POS             ID              REF     ALT     QUAL    FILTER  INFO
chrM    1189            rs28358571      T       C       .       PASS    AC=6;AF=0.011;AN=536
chrM    1442            rs28358573      G       A       .       PASS    AC=16;AF=0.030;AN=534
chrM    1888            rs28358577      G       A       .       PASS    AC=22;AF=0.041;AN=540
chrM    3450            rs28358583      C       T       .       PASS    AC=30;AF=0.056;AN=540
chrM    3480            rs28358584      A       G       .       PASS    AC=6;AF=0.011;AN=532
chrM    3666            rs28357968      G       A       .       PASS    AC=30;AF=0.056;AN=540
chrM    4883            rs28357979      C       T       .       PASS    AC=52;AF=0.097;AN=538
chrM    4917            rs28357980      A       G       .       PASS    AC=17;AF=0.032;AN=534
chrM    5285            rs28357986      A       G       .       PASS    AC=4;AF=7.407e-03;AN=540
chrM    5465            rs3902405       T       C       .       PASS    AC=6;AF=0.011;AN=540
chrM    5773            rs35855595      G       A       .       PASS    AC=26;AF=0.050;AN=524
chrM    6071            rs28358868      T       .       .       PASS    AN=530
.
.
.
chrY    23993156        rs1276032       C       A       .       PASS    AC=222;AF=0.771;AN=288
chrY    24069943        rs4025796       T       .       .       PASS    AN=288
chrY    24069946        rs3951517       C       .       .       PASS    AN=288
chrY    24070105        rs324937        G       .       .       PASS    AN=284
chrY    24359931        rs9786095       C       T       .       PASS    AC=140;AF=0.486;AN=288
chrY    24401940        rs9786258       T       C       .       PASS    AC=100;AF=0.347;AN=288
chrY    24475669        rs17316910      G       T       .       PASS    AC=4;AF=0.014;AN=288
chrY    28538621        rs2890638       C       .       .       PASS    AN=284
chrY    28582937        rs16981830      A       G       .       PASS    AC=100;AF=0.347;AN=288
chrY    28606269        rs9306850       C       T       .       PASS    AC=226;AF=0.785;AN=288
chrY    28612323        rs4141927       G       A       .       PASS    AC=224;AF=0.783;AN=286
chrY    28733101        rs1358368       G       C       .       PASS    AC=224;AF=0.778;AN=288
chrY    28758193        rs9786773       A       G       .       PASS    AC=224;AF=0.778;AN=288
```

### OMNI

#### Sourced from the [GATK resource bundle](https://gatkforums.broadinstitute.org/gatk/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it)

```bash
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/1000G_omni2.5.hg19.sites.vcf.gz
```

```bash
#CHROM  POS     ID      REF             ALT     QUAL    FILTER  INFO
chr1    534247  SNP1-524110             C       T       .       PASS    CR=99.95311;GentrainScore=0.7423;HW=1.0
chr1    565286  SNP1-555149             C       T       .       PASS    CR=99.16279;GentrainScore=0.7029;HW=1.0
chr1    569624  SNP1-559487             T       C       .       PASS    CR=97.697975;GentrainScore=0.8070;HW=1.0
chr1    689186  rs4000335               G       A       .       NOT_POLY_IN_1000G       CR=99.906586;GentrainScore=0.7934;HW=1.0
chr1    723918  SNP1-713781             G       A       .       PASS    CR=99.9051;GentrainScore=0.4541;HW=0.2935265
chr1    729632  SNP1-719495             C       T       .       PASS    CR=99.485985;GentrainScore=0.6870;HW=0.002817232
chr1    752566  rs3094315               G       A       .       PASS    CR=99.925735;GentrainScore=0.8141;HW=8.3064125E-11
chr1    752721  rs3131972               A       G       .       PASS    CR=99.859055;GentrainScore=0.8578;HW=1.7165839E-7
chr1    754063  SNP1-743926             G       T       .       PASS    CR=99.95258;GentrainScore=0.5893;HW=0.29127777
chr1    756652  SNP1-746515             T       G       .       NOT_POLY_IN_1000G       CR=100.0;GentrainScore=0.6899;HW=1.0
.
.
.
chrY    28815253        rs9786787       G       A       .       NOT_POLY_IN_1000G       CR=48.482018;GentrainScore=0.4067;HW=1.0
chrY    28817442        rs9786160       A       C       .       NOT_POLY_IN_1000G       CR=31.013542;GentrainScore=0.0000;HW=1.0
chrY    28817458        rs9786224       C       A       .       NOT_POLY_IN_1000G       CR=47.45446;GentrainScore=0.4222;HW=1.0
chrY    28817636        rs7067479       G       A       .       PASS    CR=31.794395;GentrainScore=0.4201;HW=1.0
chrY    28817799        rs9782828       A       G       .       NOT_POLY_IN_1000G       CR=31.106964;GentrainScore=0.0000;HW=1.0
chrY    28817936        rs7067236       G       T       .       NOT_POLY_IN_1000G       CR=66.417564;GentrainScore=0.4993;HW=1.0
chrY    58855394        rs2880672       A       G       .       NOT_POLY_IN_1000G       CR=29.472206;GentrainScore=0.0000;HW=1.0
chrY    58865960        rs1832424       G       A       .       NOT_POLY_IN_1000G       CR=49.41616;GentrainScore=0.3950;HW=1.0
chrY    58883690        rs9786720       T       C       .       NOT_POLY_IN_1000G       CR=31.293785;GentrainScore=0.0000;HW=1.0
chrY    58969089        rs11799070      T       C       .       NOT_POLY_IN_1000G       CR=97.24428;GentrainScore=0.5542;HW=1.0
chrY    59032197        rs28715603      T       C       .       PASS    CR=44.988186;GentrainScore=0.0000;HW=1.0
```

### 1000G indels

#### Sourced from the [GATK resource bundle](https://gatkforums.broadinstitute.org/gatk/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it)

```bash
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.gz
```

```bash
#CHROM  POS             ID      REF     ALT     QUAL    FILTER  INFO
chr1    13957           .       TC      T       28      PASS    AC=35;AF=0.02;AFR_AF=0.02;AMR_AF=0.02;AN=2184;ASN_AF=0.01;AVGPOST=0.8711;ERATE=0.0065;EUR_AF=0.02;LDAF=0.0788;RSQ=0.2501;THETA=0.0100;VT=INDEL
chr1    46402           .       C       CTGT    31      PASS    AC=8;AF=0.0037;AFR_AF=0.01;AN=2184;ASN_AF=0.0017;AVGPOST=0.8325;ERATE=0.0072;LDAF=0.0903;RSQ=0.0960;THETA=0.0121;VT=INDEL
chr1    47190           .       G       GA      192     PASS    AC=29;AF=0.01;AFR_AF=0.06;AMR_AF=0.0028;AN=2184;AVGPOST=0.9041;ERATE=0.0041;LDAF=0.0628;RSQ=0.2883;THETA=0.0153;VT=INDEL
chr1    52185           .       TTAA    T       244     PASS    AC=10;AF=0.0046;AFR_AF=0.0020;AMR_AF=0.02;AN=2184;ASN_AF=0.0035;AVGPOST=0.9840;ERATE=0.0037;LDAF=0.0124;RSQ=0.4271;THETA=0.0232;VT=INDEL
chr1    53234           .       CAT     C       227     PASS    AC=10;AF=0.0046;AFR_AF=0.02;AMR_AF=0.0028;AN=2184;AVGPOST=0.9936;ERATE=0.0007;LDAF=0.0074;RSQ=0.6237;THETA=0.0119;VT=INDEL
chr1    55249           .       C       CTATGG  443     PASS    AC=151;AF=0.07;AFR_AF=0.03;AMR_AF=0.08;AN=2184;ASN_AF=0.16;AVGPOST=0.9073;ERATE=0.0063;EUR_AF=0.02;LDAF=0.0968;RSQ=0.5891;THETA=0.0038;VT=INDEL
chr1    63735           .       CCTA    C       455     PASS    AC=829;AF=0.38;AFR_AF=0.13;AMR_AF=0.33;AN=2184;ASN_AF=0.69;AVGPOST=0.7654;ERATE=0.0047;EUR_AF=0.34;LDAF=0.4128;RSQ=0.6424;THETA=0.0062;VT=INDEL
chr1    72119           .       G       GTA     158     PASS    AC=8;AF=0.0037;AMR_AF=0.0028;AN=2184;ASN_AF=0.01;AVGPOST=0.9589;ERATE=0.0026;EUR_AF=0.0013;LDAF=0.0243;RSQ=0.2268;THETA=0.0016;VT=INDEL
chr1    72297           .       G       GTAT    160     PASS    AC=19;AF=0.01;AMR_AF=0.02;AN=2184;ASN_AF=0.01;AVGPOST=0.9383;ERATE=0.0055;EUR_AF=0.01;LDAF=0.0399;RSQ=0.3194;THETA=0.0064;VT=INDEL
chr1    84005           .       AG      A       78      PASS    AC=52;AF=0.02;AFR_AF=0.02;AMR_AF=0.03;AN=2184;ASN_AF=0.01;AVGPOST=0.9360;ERATE=0.0049;EUR_AF=0.04;LDAF=0.0514;RSQ=0.4690;THETA=0.0005;VT=INDEL

.
.
.
chrX    155207762       .       C       CT      833     PASS    AC=2174;AF=1.00;AFR_AF=0.99;AMR_AF=1.00;AN=2184;ASN_AF=1.00;AVGPOST=0.9951;ERATE=0.0006;EUR_AF=1.00;LDAF=0.9953;RSQ=0.6062;THETA=0.0004;VT=INDEL
chrX    155209404       .       T       TA      291     PASS    AC=1893;AF=0.87;AFR_AF=0.96;AMR_AF=0.86;AN=2184;ASN_AF=0.90;AVGPOST=0.8395;ERATE=0.0290;EUR_AF=0.79;LDAF=0.8202;RSQ=0.5841;THETA=0.0004;VT=INDEL
chrX    155215152       .       TC      T       670     PASS    AC=2075;AF=0.95;AFR_AF=0.95;AMR_AF=0.95;AN=2184;ASN_AF=0.94;AVGPOST=0.9298;ERATE=0.0138;EUR_AF=0.96;LDAF=0.9292;RSQ=0.5851;THETA=0.0007;VT=INDEL
chrX    155219875       .       GT      G       172     PASS    AC=61;AF=0.03;AFR_AF=0.11;AMR_AF=0.01;AN=2184;AVGPOST=0.9982;ERATE=0.0003;LDAF=0.0286;RSQ=0.9758;THETA=0.0005;VT=INDEL
chrX    155222251       .       G       GC      703     PASS    AC=2156;AF=0.99;AFR_AF=0.97;AMR_AF=0.99;AN=2184;ASN_AF=0.99;AVGPOST=0.9829;ERATE=0.0013;EUR_AF=1.00;LDAF=0.9856;RSQ=0.5953;THETA=0.0004;VT=INDEL
chrX    155224341       .       G       GC      641     PASS    AC=2170;AF=0.99;AFR_AF=1.00;AMR_AF=1.00;AN=2184;ASN_AF=0.99;AVGPOST=0.9772;ERATE=0.0032;EUR_AF=1.00;LDAF=0.9843;RSQ=0.4063;THETA=0.0005;VT=INDEL
chrX    155231142       .       GT      G       541     PASS    AC=119;AF=0.05;AFR_AF=0.22;AMR_AF=0.02;AN=2184;AVGPOST=0.9960;ERATE=0.0006;EUR_AF=0.0013;LDAF=0.0543;RSQ=0.9723;THETA=0.0006;VT=INDEL
chrX    155231738       .       T       TG      884     PASS    AC=2172;AF=0.99;AFR_AF=0.99;AMR_AF=1.00;AN=2184;ASN_AF=1.00;AVGPOST=0.9871;ERATE=0.0032;EUR_AF=0.99;LDAF=0.9907;RSQ=0.4336;THETA=0.0007;VT=INDEL
chrX    155235935       .       CTG     C       522     PASS    AC=37;AF=0.02;AFR_AF=0.07;AMR_AF=0.01;AN=2184;AVGPOST=0.9931;ERATE=0.0004;LDAF=0.0188;RSQ=0.8674;THETA=0.0036;VT=INDEL
chrX    155247518       .       TC      T       245     PASS    AC=2163;AF=0.99;AFR_AF=0.99;AMR_AF=0.99;AN=2184;ASN_AF=0.99;AVGPOST=0.7146;ERATE=0.0154;EUR_AF=0.99;LDAF=0.8371;RSQ=0.0835;THETA=0.0162;VT=INDEL
chrX    155248401       .       TG      T       152     PASS    AC=2176;AF=1.00;AFR_AF=0.99;AMR_AF=0.99;AN=2184;ASN_AF=1.00;AVGPOST=0.8914;ERATE=0.0049;EUR_AF=1.00;LDAF=0.9409;RSQ=0.1274;THETA=0.0065;VT=INDEL
chrX    155259685       .       TA      T       318     PASS    AC=150;AF=0.07;AFR_AF=0.06;AMR_AF=0.09;AN=2184;ASN_AF=0.08;AVGPOST=0.6898;ERATE=0.0309;EUR_AF=0.06;LDAF=0.2296;RSQ=0.2641;THETA=0.0069;VT=INDEL
```

### 1000G snps

#### Sourced from the [GATK resource bundle](https://gatkforums.broadinstitute.org/gatk/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it)

```bash
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
```

```bash
#CHROM  POS     ID              REF     ALT     QUAL            FILTER  INFO
chr1    51479   rs116400033     T       A       11726.81        PASS    AC=229;AF=0.3253;AN=704;BaseQRankSum=-6.949;DB;DP=1570;Dels=0.00;FS=3.130;HRun=0;HaplotypeScore=0.1377;InbreedingCoeff=0.2907;MQ=34.37;MQ0=174;MQRankSum=1.476;QD=16.08;ReadPosRankSum=-0.202;SB=-4317.78;VQSLOD=5.1635;pop=EUR.admix
chr1    55367   .               G       A       207.20          PASS    AC=2;AF=0.00117;AN=1714;BaseQRankSum=2.243;DP=4926;Dels=0.00;FS=3.005;HRun=0;HaplotypeScore=0.1382;InbreedingCoeff=-0.0188;MQ=45.57;MQ0=365;MQRankSum=0.185;QD=21.22;ReadPosRankSum=0.136;SB=-111.01;VQSLOD=6.3979;pop=ALL
chr1    55388   .               C       T       95.61           PASS    AC=1;AF=0.00056;AN=1792;BaseQRankSum=-0.038;DP=5282;Dels=0.00;FS=0.000;HRun=2;HaplotypeScore=0.1980;InbreedingCoeff=-0.0278;MQ=48.19;MQ0=20;MQRankSum=-0.397;QD=18.13;ReadPosRankSum=-0.945;SB=-59.46;VQSLOD=5.7297;pop=ALL
chr1    55852   .               G       C       503.17          PASS    AC=6;AF=0.0080;AN=748;BaseQRankSum=1.052;DP=1772;Dels=0.00;FS=1.623;HRun=0;HaplotypeScore=0.1542;InbreedingCoeff=-0.0547;MQ=39.86;MQ0=114;MQRankSum=1.243;QD=9.15;ReadPosRankSum=1.657;SB=-289.09;VQSLOD=5.8416;pop=EUR.admix
chr1    61462   rs56992750      T       A       2671.41         PASS    AC=70;AF=0.04679;AN=1496;BaseQRankSum=-0.165;DB;DP=3256;Dels=0.00;FS=1.789;HRun=1;HaplotypeScore=0.1091;InbreedingCoeff=-0.0562;MQ=36.49;MQ0=395;MQRankSum=2.142;QD=6.50;ReadPosRankSum=0.789;SB=-1462.83;VQSLOD=6.2119;pop=ALL
chr1    62157   rs10399597      G       A       285.13          PASS    AC=5;AF=0.00274;AN=1826;BaseQRankSum=3.510;DB;DP=5819;Dels=0.00;FS=0.802;HRun=0;HaplotypeScore=0.2357;InbreedingCoeff=-0.0660;MQ=39.50;MQ0=993;MQRankSum=2.274;QD=7.11;ReadPosRankSum=0.978;SB=-213.74;VQSLOD=5.1283;pop=ALL
chr1    82609   .               C       G       1821.30         PASS    AC=53;AF=0.0716;AN=740;BaseQRankSum=-2.473;DP=1820;Dels=0.00;FS=4.066;HRun=1;HaplotypeScore=0.1207;InbreedingCoeff=0.1324;MQ=53.52;MQ0=320;MQRankSum=0.750;QD=9.59;ReadPosRankSum=-0.216;SB=-1153.72;VQSLOD=6.8223;pop=EUR.admix
chr1    82734   rs4030331       T       C       12360.27        PASS    AC=316;AF=0.20654;AN=1530;BaseQRankSum=-7.774;DB;DP=3665;Dels=0.00;FS=3.305;HRun=0;HaplotypeScore=0.1481;InbreedingCoeff=0.0578;MQ=52.42;MQ0=272;MQRankSum=0.086;QD=8.82;ReadPosRankSum=2.068;SB=-3526.93;VQSLOD=6.9968;pop=ALL
chr1    83084   .               T       A       29167.41        PASS    AC=726;AF=0.8462;AN=858;BaseQRankSum=3.752;DP=1533;Dels=0.00;FS=6.435;HRun=3;HaplotypeScore=0.0799;InbreedingCoeff=0.2112;MQ=46.84;MQ0=108;MQRankSum=1.586;QD=24.11;ReadPosRankSum=1.138;SB=-17136.27;VQSLOD=5.5460;pop=ALL
.
.
.
chr1    61462   rs56992750      T       A       2671.41         PASS    AC=70;AF=0.04679;AN=1496;BaseQRankSum=-0.165;DB;DP=3256;Dels=0.00;FS=1.789;HRun=1;HaplotypeScore=0.1091;InbreedingCoeff=-0.0562;MQ=36.49;MQ0=395;MQRankSum=2.142;QD=6.50;ReadPosRankSum=0.789;SB=-1462.83;VQSLOD=6.2119;pop=ALL
chr1    62157   rs10399597      G       A       285.13          PASS    AC=5;AF=0.00274;AN=1826;BaseQRankSum=3.510;DB;DP=5819;Dels=0.00;FS=0.802;HRun=0;HaplotypeScore=0.2357;InbreedingCoeff=-0.0660;MQ=39.50;MQ0=993;MQRankSum=2.274;QD=7.11;ReadPosRankSum=0.978;SB=-213.74;VQSLOD=5.1283;pop=ALL
chr1    82609   .               C       G       1821.30         PASS    AC=53;AF=0.0716;AN=740;BaseQRankSum=-2.473;DP=1820;Dels=0.00;FS=4.066;HRun=1;HaplotypeScore=0.1207;InbreedingCoeff=0.1324;MQ=53.52;MQ0=320;MQRankSum=0.750;QD=9.59;ReadPosRankSum=-0.216;SB=-1153.72;VQSLOD=6.8223;pop=EUR.admix
chr1    82734   rs4030331       T       C       12360.27        PASS    AC=316;AF=0.20654;AN=1530;BaseQRankSum=-7.774;DB;DP=3665;Dels=0.00;FS=3.305;HRun=0;HaplotypeScore=0.1481;InbreedingCoeff=0.0578;MQ=52.42;MQ0=272;MQRankSum=0.086;QD=8.82;ReadPosRankSum=2.068;SB=-3526.93;VQSLOD=6.9968;pop=ALL
chr1    83084   .               T       A       29167.41        PASS    AC=726;AF=0.8462;AN=858;BaseQRankSum=3.752;DP=1533;Dels=0.00;FS=6.435;HRun=3;HaplotypeScore=0.0799;InbreedingCoeff=0.2112;MQ=46.84;MQ0=108;MQRankSum=1.586;QD=24.11;ReadPosRankSum=1.138;SB=-17136.27;VQSLOD=5.5460;pop=ALL
chr1    83771   .               T       G       663.46          PASS    AC=22;AF=0.01391;AN=1582;BaseQRankSum=-4.181;DP=4131;Dels=0.00;FS=3.802;HRun=0;HaplotypeScore=0.1154;InbreedingCoeff=-0.0062;MQ=38.11;MQ0=947;MQRankSum=-0.593;QD=6.27;ReadPosRankSum=0.038;SB=-368.75;VQSLOD=5.6475;pop=ALL
chr1    86028   rs114608975     T       C       2389.90         PASS    AC=55;AF=0.0688;AN=800;BaseQRankSum=-15.243;DB;DP=1996;Dels=0.00;FS=0.000;HRun=0;HaplotypeScore=0.4883;InbreedingCoeff=0.1564;MQ=61.98;MQ0=6;MQRankSum=1.939;QD=10.30;ReadPosRankSum=-0.141;SB=-1659.01;VQSLOD=5.8942;pop=EUR.admix
chr1    86065   rs116504101     G       C       2600.45         PASS    AC=54;AF=0.0705;AN=766;BaseQRankSum=1.016;DB;DP=2073;Dels=0.00;FS=3.978;HRun=0;HaplotypeScore=0.1110;InbreedingCoeff=0.1498;MQ=54.67;MQ0=50;MQRankSum=2.517;QD=11.41;ReadPosRankSum=-1.355;SB=-1245.00;VQSLOD=5.4933;pop=EUR.admix
chr1    86282   .               T       G       82.19           PASS    AC=1;AF=0.00066;AN=1516;BaseQRankSum=1.121;DP=4113;Dels=0.00;FS=4.450;HRun=0;HaplotypeScore=0.1230;InbreedingCoeff=-0.0713;MQ=50.47;MQ0=377;MQRankSum=-0.472;QD=15.00;ReadPosRankSum=-0.047;SB=-31.80;VQSLOD=5.1323;pop=ALL
chr1    86331   rs115209712     A       G       15414.55        PASS    AC=218;AF=0.12616;AN=1728;BaseQRankSum=-4.790;DB;DP=4696;Dels=0.00;FS=0.922;HRun=2;HaplotypeScore=0.1808;InbreedingCoeff=0.0185;MQ=58.01;MQ0=1;MQRankSum=-7.246;QD=11.72;ReadPosRankSum=1.185;SB=-5186.86;VQSLOD=5.4426;pop=ALL
chr1    87021   .               T       C       351.74          PASS    AC=8;AF=0.00503;AN=1590;BaseQRankSum=-2.753;DP=4900;Dels=0.00;FS=4.736;HRun=0;HaplotypeScore=0.1953;InbreedingCoeff=-0.0172;MQ=40.99;MQ0=570;MQRankSum=-2.250;QD=9.02;ReadPosRankSum=0.212;SB=-183.08;VQSLOD=5.1704;pop=ALL
```

### Mills

#### Sourced from the [GATK resource bundle](https://gatkforums.broadinstitute.org/gatk/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it)

```bash
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
```

```bash
#CHROM  POS     ID              REF          ALT               QUAL      FILTER  INFO
chr1    55249   .               C            CTATGG            6160.83   PASS    set=Intersect1000GMinusBI
chr1    87114   .               CT           C                 666.86    PASS    set=Intersect1000GAll
chr1    233586  .               CCT          C                 1083.97   PASS    set=Intersect1000GMinusBI
chr1    691567  .               CACAG        C                 2260.77   PASS    set=Intersect1000GMinusBI
chr1    715309  .               AG           A                 108.86    PASS    set=Intersect1000GMinusOX
chr1    722074  .               A            ATT               361.06    PASS    set=Intersect1000GMinusDI
chr1    732134  .               T            TA                1722.44   PASS    set=Intersect1000GAll
chr1    752307  .               AT           A                 2106.07   PASS    set=Intersect1000GAll
chr1    758176  .               CCT          C                 436.88    PASS    set=Intersect1000GAll
chr1    759444  rs78982110      ATATT        A                 6659.29   PASS    set=Intersect1000GAll
chr1    765084  .               TTTTG        T                 1998.11   PASS    set=Intersect1000GMinusSI
chr1    765155  rs71576582      TTGAC        T                 1466.09   PASS    set=Intersect1000GAll
chr1    769138  rs59306077      CAT          C                 29441.73  PASS    set=Intersect1000GMinusSI
.
.
.
chrY    59004091        742403  C            CGT               .         PASS    set=MillsDoubleCenter
chrY    59004094        742404  ATG          A                 .         PASS    set=MillsDoubleCenter
chrY    59006822        856660  A            AC                .         PASS    set=MillsDoubleCenter
chrY    59006878        856662  CAA          C                 .         PASS    set=MillsDoubleCenter
chrY    59014273        90884   C            CATT              .         PASS    set=MillsDoubleCenter
chrY    59014594        269468  A            ATTAAGT,ATTGT     .         PASS    set=MillsDoubleCenter-MillsTracesUnknown
chrY    59014682        269469  AAT          A                 .         PASS    set=MillsDoubleCenter
chrY    59014794        269470  T            TG                .         PASS    set=MillsDoubleCenter
chrY    59015470        717022  T            TG                .         PASS    set=MillsDoubleCenter
chrY    59016882        722048  G            GCTTCTC           .         PASS    set=MillsDoubleCenter
chrY    59030173        843507  GT           G                 .         PASS    set=MillsDoubleCenter
chrY    59030478        125411  AAAACAAAC    A,AAAAC           .         PASS    set=MillsCenterUnknown-MillsTracesUnknown
```

### Axiom exome plus

#### Sourced from the [GATK resource bundle](https://gatkforums.broadinstitute.org/gatk/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it)

```bash
## Not avaialble?
```

```bash

```

### CADD

#### Sourced from [CADD](https://cadd.gs.washington.edu/download)

```bash
aria2c -c -s 10 https://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/whole_genome_SNVs.tsv.gz
```

```bash
#Chr    Pos             Ref     Alt     RawScore        PHRED
1       10001           T       A       0.118631        4.575
1       10001           T       C       0.135541        4.848
1       10001           T       G       0.111762        4.462
1       10002           A       C       0.111699        4.461
1       10002           A       G       0.133604        4.817
1       10002           A       T       0.117529        4.557
1       10003           A       C       0.111849        4.463
1       10003           A       G       0.133753        4.819
1       10003           A       T       0.117678        4.559
1       10004           C       A       0.008692        2.756
1       10004           C       G       -0.000601       2.614
1       10004           C       T       0.055241        3.516
1       10005           C       A       0.008520        2.754
1       10005           C       G       -0.000772       2.611
.
.
.
Y       59034045        T       C       0.338613        7.696
Y       59034045        T       G       0.314834        7.405
Y       59034046        T       A       0.321877        7.493
Y       59034046        T       C       0.338786        7.698
Y       59034046        T       G       0.315007        7.407
Y       59034047        T       A       0.321730        7.491
Y       59034047        T       C       0.338639        7.696
Y       59034047        T       G       0.314860        7.405
Y       59034048        T       A       0.321903        7.493
Y       59034048        T       C       0.338813        7.698
Y       59034048        T       G       0.315034        7.407
Y       59034049        C       A       0.212010        6.017
Y       59034049        C       G       0.202718        5.882
Y       59034049        C       T       0.258559        6.667
```

## GRCh38/hg38

### Reference human genome

#### Sourced from [NCBI](https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/)

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
```

```bash
>NC_000001.11 Homo sapiens chromosome 1, GRCh38.p13 Primary Assembly
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
.
.
.
>NC_012920.1 Homo sapiens mitochondrion, complete genome
GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTATGCACGC
GATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGATTCCTGCCTCATCCTATTATTTA
TCGCACCTACGTTCAATATTACAGGCGAACATACTTACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATA
ACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCAAACCCCCCCTCCCCCGCTTC
TGGCCACAGCACTTAAACACATCTCTGCCAAACCCCAAAAACAAAGAACCCTAACACCAGCCTAACCAGATTTCAAATTT
TATCTTTTGGCGGTATGCACTTTTAACAGTCACCCCCCAACTAACACATTATTTTCCCCTCCCACTCCCATACTACTAAT
CTCATCAATACAACCCCCGCCCATCCTACCCAGCACACACACACCGCTGCTAACCCCATACCCCGAACCAACCAAACCCC
AAAGACACCCCCCACAGTTTATGTAGCTTACCTCCTCAAAGCAATACACTGAAAATGTTTAGACGGGCTCACATCACCCC
ATAAACAAATAGGTTTGGTCCTAGCCTTTCTATTAGCTCTTAGTAAGATTACACATGCAAGCATCCCCGTTCCAGTGAGT
TCACCCTCTAAATCACCACGATCAAAAGGAACAAGCATCAAGCACGCAGCAATGCAGCTCAAAACGCTTAGCCTAGCCAC
ACCCCCACGGGAAACAGCAGTGATTAACCTTTAGCAATAAACGAAAGTTTAACTAAGCTATACTAACCCCAGGGTTGGTC
AATTTCGTGCCAGCCACCGCGGTCACACGATTAACCCAAGTCAATAGAAGCCGGCGTAAAGAGTGTTTTAGATCACCCCC
.
.
.
```

#### Sourced from the [GATK resource bundle](https://gatkforums.broadinstitute.org/gatk/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it)

```bash
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.gz
```

```bash
>chr1  AC:CM000663.2  gi:568336023  LN:248956422  rl:Chromosome  M5:6aef897c3d6ff0c78aff06ac189178dd  AS:GRCh38
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
.
.
.
>HLA-DRB1*16:02:01      HLA00878 11005 bp
GATGCTGATTGGTTCTCCAACACGAGATTACCCAACCCAGGAGCAAGGAAATCAGTAACTTCCTCCCTATAACTTGGAATGTGGGTGGAGGGGTTCATAG
TTCTCCCTGAGTGAGACTTGCCTGCTTCTCTGGCCCCTGGTCCTGTCCTGTTCTCCAGCATGGTGTGTCTGAAGCTCCCTGGAGGCTCCTGCATGACAGC
GCTGACAGTGACACTGATGGTGCTGAGCTCCCCACTGGCTTTGGCTGGGGACACCCGACGTAAGTGCACATTGCGGGTGCTGAGCTACTATGGGGTGGGG
AAAATAGGGAGTTTTGTTAACATTGTGCCCAGGCCATGTCCCTTAAGAAATTGTGACGTTTTCTTCAGAGATTGCCCATCTTTATCATTGGATCCCAAAT
TATTTCCTCCATAAAAGGAGCTTGGGTACTTGCCCTCTTCATGAGACTTGTGTAAGGGGCCTTTGCACAAGTCATTTCTTTTCAAATCTCCACCAATAAA
ACCTTTGCCTCACATGTCCTCAGGGTCTTTAGAGGATTTAGAAATAAGGATTCTAAAATAAATTCCCCATACAGCACTTCCCTTTATTATGTTGACTTAT
GTCAGACAAAAGGAGGTTCTTACTGAAAATTTTGTGGGAGTCAAGGGAATTCAAAGGGTCTCTCCTAGACGATCCAGTGTTAGGTTCCCCACAGGACCTT
TGGTGTTGGCCATAGTCCTCATATGTGAGGATGGACCCAGTGGCCTCCCCATTATCTCCTTTCTTTTCTTGCTGAACTCCAATGTTTATAAGGCCTGTAT
CCCTGTAGCGTATGTAGGTTCTCTGACAGAAGTTATACTTAGTGCTCTTCCTTTCTTGTGGGGAAAAGTCCCTGGAACTGAAGCTGAGATTGTTAGTACT
TGGAGTCACCTTACAGATACAGAGCATTTATGAGGTATTCTTTGGTGCCTAAAGAACTTAAGGCATCCTCTGAAAAACTAGCCCAGGTTCGTGTTCATTA
TGAATCTTTTTTAACCTTTCTGTACTTGTTTCTCTTGCATCTCCTATGTGCTCTAACTAGACATGACAGAAGAGATTTAACTAATGTATAAATTATATGA
AATTCTATTTTTTAAGTCAAAAATAATCAACTATCAGAAATTTAATAATGTTCAAACTATATACTCTGTGTGGGGTTACCGAGATGATGTGAACATTGTT
CACGTCTCATAGGGCTGAAAGTCAATGGGCAAGTCTTGGGAACTCATTGTCTTACTGGGGTCTTGTCCTAAATTTCCTAGGTTCACCCATCATGCCCTCA
GCTTTCCTTAACTAGCCATGTCTGCTTACCTCTTCCTCCAGTTTCTATTTTTCCCCAGCTATGTTGTCATCATTTCCAGAAATCTCTAAAGCTTGCACAG
.
.
.

```

#### Sourced from [UCSC](https://hgdownload.soe.ucsc.edu/downloads.html)

```bash
wget ftp://hgdownload.soe.ucsc.edu:21/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
```

```bash
>chr1
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
.
.
.
>chrY_KI270740v1_random
TAATAAATTTTGAAGAAAATGAAGACTGTGTTCTCAGTTCCAGGTGCTTC
ATCAGGCTCATTGTGGATCCAGACTACCAGACACAAGACATTACACATTG
TAATGCATTAAATGCATAGTTTTAACAGTAATAATTTAAAAGAGATTTAG
AATTTTATAATGTTTGGAAAAATACATAGAGGCTTACTTTTTATTTTATT
TTTTTGAGATAGGAAGCCtttttttttgtttttgtttttgtttctgtttt
tgttttttgagacagagtctcaccatgtcacccagactggagtgcagtgg
tgcaatatcggcccattgcaagctccacatcccaggttcacaccattctc
ctgcctcagcctcccaagtagctgggactacaggtgcccgccaccacatc
cagctaatttttttttgtacttttagtagagacggggtatcaccatgtga
gccaagatggtctccatctcctgacctcgtgatctgcccaccttggcctc
ccaaagtgctgggattacaggggtgagccaccacgcccagGCATAGAGGC
ACTTTTAACCATAAATGAACACTGTTATGATTTGTATTACCACAGTATCA
TTATTCTGTCCTGTTTGCCTTACAttttatttatttattatactgtaagt
tctgggatacatgtgcagaatgtgcaggtttgttacagagatatatgctt
gtttgctgcacctgtcagtttttcatctacattaggtatttctcctaatg
ctattccctgttaggtccccaccctccaacagtctccagtgtttgatgtt
cccctccctatgtccatgtattctcattttacaactcccacctatgagtg
agaaattgcagtgtttgTGtgtttggaacttattccttccagtgggtttg
.
.
.
```

### dbSNP

#### Sourced from [NCBI](https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/)

```bash
# Build 153
wget ftp://ftp.ncbi.nlm.nih.gov:21/snp/latest_release/VCF/GCF_000001405.38.gz
```

```bash
#CHROM          POS     ID              REF     ALT                                             QUAL    FILTER  INFO
NC_000001.11    10019   rs775809821     TA      T                                               .       .       RS=775809821;dbSNPBuildID=144;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=INDEL
NC_000001.11    10039   rs978760828     A       C                                               .       .       RS=978760828;dbSNPBuildID=150;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV
NC_000001.11    10043   rs1008829651    T       A                                               .       .       RS=1008829651;dbSNPBuildID=150;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV
NC_000001.11    10051   rs1052373574    A       G                                               .       .       RS=1052373574;dbSNPBuildID=150;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV
NC_000001.11    10051   rs1326880612    A       AC                                              .       .       RS=1326880612;dbSNPBuildID=151;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=INDEL
NC_000001.11    10055   rs768019142     T       TA                                              .       .       RS=768019142;dbSNPBuildID=144;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=INDEL
NC_000001.11    10055   rs892501864     T       A                                               .       .       RS=892501864;dbSNPBuildID=150;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV
NC_000001.11    10063   rs1010989343    A       C                                               .       .       RS=1010989343;dbSNPBuildID=150;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV
NC_000001.11    10067   rs1489251879    T       TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC      .       .       RS=1489251879;dbSNPBuildID=151;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=INDEL
NC_000001.11    10077   rs1022805358    C       G                                               .       .       RS=1022805358;dbSNPBuildID=150;SSR=0;PSEUDOGENEINFO=DDX11L1:100287102;VC=SNV
.
.
.
NW_019805503.1  163143  rs1236207813    AC      A                                               .       .       RS=1236207813;dbSNPBuildID=151;SSR=0;VC=DEL;GNO;FREQ=GnomAD:1,3.186e-05|TOPMED:1,7.964e-06
NW_019805503.1  163145  rs778248020     T       A,C                                             .       .       RS=778248020;dbSNPBuildID=144;SSR=0;VC=SNV;GNO;FREQ=ALSPAC:1,.,0|TOPMED:1,7.964e-06,.|TWINSUK:0.9997,.,0.0002697
NW_019805503.1  163153  rs1305013879    G       T                                               .       .       RS=1305013879;dbSNPBuildID=151;SSR=0;VC=SNV;GNO;FREQ=TOPMED:1,7.964e-06
NW_019805503.1  163159  rs573488209     C       T                                               .       .       RS=573488209;dbSNPBuildID=142;SSR=0;VC=SNV;GNO;FREQ=TOPMED:0.9999,5.575e-05
NW_019805503.1  163160  rs566808919     G       A                                               .       .       RS=566808919;dbSNPBuildID=142;SSR=0;VC=SNV;GNO;FREQ=1000Genomes:0.9996,0.0003994|ALSPAC:0.9997,0.0002595|Estonian:0.9996,0.0004464|GnomAD:0.9991,0.0008924|NorthernSweden:0.9833,0.01667|TOPMED:0.9992,0.0008043|TWINSUK:0.9997,0.0002697
NW_019805503.1  163164  rs539135532     C       A,T                                             .       .       RS=539135532;dbSNPBuildID=142;SSR=0;VC=SNV;GNO;FREQ=1000Genomes:0.9998,.,0.0001997
NW_019805503.1  163167  rs1567898395    T       C                                               .       .       RS=1567898395;dbSNPBuildID=153;SSR=0;VC=SNV;GNO;FREQ=Vietnamese:0.9953,0.004673
NW_019805503.1  163169  rs948909810     T       A                                               .       .       RS=948909810;dbSNPBuildID=150;SSR=0;VC=SNV
NW_019805503.1  163176  rs559138705     G       A                                               .       .       RS=559138705;dbSNPBuildID=142;SSR=0;VC=SNV;GNO;FREQ=1000Genomes:0.9998,0.0001997|ALSPAC:0.9995,0.0005189|GnomAD:1,3.186e-05|TOPMED:1,3.186e-05|TWINSUK:0.9995,0.0005394
NW_019805503.1  163177  rs542557269     A       G                                               .       .       RS=542557269;dbSNPBuildID=142;SSR=0;VC=SNV;GNO;FREQ=ALSPAC:0.9995,0.0005189|TOPMED:1,4.778e-05|TWINSUK:1,0
NW_019805503.1  163183  rs1159174310    G       C                                               .       .       RS=1159174310;dbSNPBuildID=151;SSR=0;VC=SNV;GNO;FREQ=TOPMED:1,7.964e-06
NW_019805503.1  163185  rs1418195656    T       A                                               .       .       RS=1418195656;dbSNPBuildID=151;SSR=0;VC=SNV;GNO;FREQ=TOPMED:1,7.964e-06
```

#### Sourced from the [GATK resource bundle](https://gatkforums.broadinstitute.org/gatk/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it)

```bash
# Build 146
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/dbsnp_146.hg38.vcf.gz
```

```bash
#CHROM  POS     ID              REF     ALT     QUAL    FILTER  INFO
chr1    10019   rs775809821     TA      T       .       .       RS=775809821;RSPOS=10020;dbSNPBuildID=144;SSR=0;SAO=0;VP=0x050000020005000002000200;GENEINFO=DDX11L1:100287102;WGT=1;VC=DIV;R5;ASP
chr1    10055   rs768019142     T       TA      .       .       RS=768019142;RSPOS=10055;dbSNPBuildID=144;SSR=0;SAO=0;VP=0x050000020005000002000200;GENEINFO=DDX11L1:100287102;WGT=1;VC=DIV;R5;ASP
chr1    10108   rs62651026      C       T       .       .       RS=62651026;RSPOS=10108;dbSNPBuildID=129;SSR=0;SAO=0;VP=0x050000020005000002000100;GENEINFO=DDX11L1:100287102;WGT=1;VC=SNV;R5;ASP
chr1    10109   rs376007522     A       T       .       .       RS=376007522;RSPOS=10109;dbSNPBuildID=138;SSR=0;SAO=0;VP=0x050000020005000002000100;GENEINFO=DDX11L1:100287102;WGT=1;VC=SNV;R5;ASP
chr1    10128   rs796688738     A       AC      .       .       RS=796688738;RSPOS=10128;dbSNPBuildID=146;SSR=0;SAO=0;VP=0x050000020005000002000200;GENEINFO=DDX11L1:100287102;WGT=1;VC=DIV;R5;ASP
chr1    10139   rs368469931     A       T       .       .       RS=368469931;RSPOS=10139;dbSNPBuildID=138;SSR=0;SAO=0;VP=0x050000020005000002000100;GENEINFO=DDX11L1:100287102;WGT=1;VC=SNV;R5;ASP
chr1    10144   rs144773400     TA      T       .       .       RS=144773400;RSPOS=10145;dbSNPBuildID=134;SSR=0;SAO=0;VP=0x050000020005000002000200;GENEINFO=DDX11L1:100287102;WGT=1;VC=DIV;R5;ASP
chr1    10146   rs779258992     AC      A       .       .       RS=779258992;RSPOS=10147;dbSNPBuildID=144;SSR=0;SAO=0;VP=0x050000020005000002000200;GENEINFO=DDX11L1:100287102;WGT=1;VC=DIV;R5;ASP
chr1    10150   rs371194064     C       T       .       .       RS=371194064;RSPOS=10150;dbSNPBuildID=138;SSR=0;SAO=0;VP=0x050000020005000002000100;GENEINFO=DDX11L1:100287102;WGT=1;VC=SNV;R5;ASP
.
.
.
chrM    16362   rs62581341      T       C       .       .       RS=62581341;RSPOS=16362;dbSNPBuildID=129;SSR=0;SAO=0;VP=0x050000000005000102000100;WGT=1;VC=SNV;ASP;GNO
chrM    16390   rs41378955      G       A       .       .       RS=41378955;RSPOS=16390;dbSNPBuildID=127;SSR=0;SAO=0;VP=0x050000000005000402000100;WGT=1;VC=SNV;ASP;HD
chrM    16391   rs34301918      G       C,T     .       .       RS=34301918;RSPOS=16391;RV;dbSNPBuildID=126;SSR=0;SAO=0;VP=0x050000000005000402000110;WGT=1;VC=SNV;ASP;HD;NOC
chrM    16399   rs139001869     A       G       .       .       RS=139001869;RSPOS=16399;dbSNPBuildID=134;SSR=0;SAO=0;VP=0x050000000005000002000100;WGT=1;VC=SNV;ASP
chrM    16445   rs371960162     T       C       .       .       RS=371960162;RSPOS=16445;dbSNPBuildID=138;SSR=0;SAO=0;VP=0x050000000005000002000100;WGT=1;VC=SNV;ASP
chrM    16499   rs376846509     A       G       .       .       RS=376846509;RSPOS=16499;dbSNPBuildID=138;SSR=0;SAO=0;VP=0x050000020005000002000100;GENEINFO=ND6:4541;WGT=1;VC=SNV;R5;ASP
chrM    16512   rs373943637     T       C       .       .       RS=373943637;RSPOS=16512;dbSNPBuildID=138;SSR=0;SAO=0;VP=0x050000020005000002000100;GENEINFO=ND6:4541;WGT=1;VC=SNV;R5;ASP
chrM    16519   rs3937033       T       C       .       .       RS=3937033;RSPOS=16519;dbSNPBuildID=108;SSR=0;SAO=0;VP=0x050100020005000502000100;GENEINFO=ND6:4541;WGT=1;VC=SNV;SLO;R5;ASP;HD;GNO
chrM    16528   rs386829315     C       T       .       .       RS=386829315;RSPOS=16528;RV;dbSNPBuildID=138;SSR=0;SAO=0;VP=0x050000020005000002000100;GENEINFO=ND6:4541;WGT=1;VC=SNV;R5;ASP
chrM    16529   rs370705831     T       C       .       .       RS=370705831;RSPOS=16529;dbSNPBuildID=138;SSR=0;SAO=0;VP=0x050000020005000002000100;GENEINFO=ND6:4541;WGT=1;VC=SNV;R5;ASP
```

### Hapmap

#### Sourced from the [GATK resource bundle](https://gatkforums.broadinstitute.org/gatk/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it)

```bash
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/hapmap_3.3.hg38.vcf.gz
```

```bash
#CHROM                  POS     ID              REF     ALT     QUAL    FILTER  INFO
chr1                    55299   rs10399749      C       .       .       PASS    AN=510
chr1                    55394   rs2949420       T       .       .       PASS    AN=178
chr1                    55550   rs2949421       A       T       .       PASS    AC=173;AF=0.972;AN=178
chr1                    56981   rs2691310       C       .       .       PASS    AN=180
chr1                    82571   rs4030303       G       A       .       PASS    AC=1;AF=1.859e-03;AN=538
chr1                    82652   rs4030300       C       .       .       PASS    AN=524
chr1                    87826   rs3855952       A       .       .       PASS    AN=318
chr1                    88169   rs940550        C       T       .       PASS    AC=355;AF=0.992;AN=358
chr1                    91605   rs13328714      C       .       .       PASS    AN=534
chr1                    107352  rs4124251       T       A,G     .       PASS    AC=416,194;AF=0.682,0.318;AN=610
chr1                    262463  rs11490937      G       .       .       PASS    AN=536
chr1                    264562  rs8179466       C       T       .       PASS    AC=57;AF=0.343;AN=166
chr1                    417474  rs9439424       T       .       .       PASS    AN=176
chr1                    599203  rs6683466       C       .       .       PASS    AN=538
chr1                    605060  rs1538941       T       .       .       PASS    AN=358
.
.
.
chr2_KI270894v1_alt     199811  rs7594880       A       .       .       PASS    AN=540
chr2_KI270894v1_alt     200051  rs6421783       G       A       .       PASS    AC=905;AF=0.999;AN=906
chr2_KI270894v1_alt     202146  rs6737667       T       C       .       PASS    AC=180;AF=1.00;AN=180
chr2_KI270894v1_alt     202303  rs677393        A       C       .       PASS    AC=538;AF=1.00;AN=538
chr2_KI270894v1_alt     203906  rs2601500       G       C       .       PASS    AC=538;AF=1.00;AN=538
chr2_KI270894v1_alt     210247  rs6421792       A       G       .       PASS    AC=180;AF=1.00;AN=180
chr19_KI270938v1_alt    195720  rs16985497      C       A       .       PASS    AC=1;AF=1.969e-03;AN=508
chr19_KI270938v1_alt    272313  rs11574615      C       T       .       PASS    AC=753;AF=0.639;AN=1178
chr19_KI270938v1_alt    272589  rs7245918       T       G       .       PASS    AC=540;AF=1.00;AN=540
chr19_KI270938v1_alt    273219  rs1975046       C       T       .       PASS    AC=136;AF=0.384;AN=354
chr19_KI270938v1_alt    274493  rs8100319       C       T       .       PASS    AC=13;AF=0.016;AN=808
chr19_KI270938v1_alt    274618  rs6509862       A       C       .       PASS    AC=171;AF=0.950;AN=180
chr19_KI270938v1_alt    275374  rs11574605      A       G       .       PASS    AC=83;AF=0.034;AN=2432
chr19_KI270938v1_alt    276439  rs3764631       A       .       .       PASS    AN=176
```

### OMNI

#### Sourced from the [GATK resource bundle](https://gatkforums.broadinstitute.org/gatk/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it)

```bash
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/1000G_omni2.5.hg38.vcf.gz
```

```bash
#CHROM                  POS     ID              REF     ALT     QUAL    FILTER                  INFO
chr1                    598867  SNP1-524110     C       T       .       PASS                    CR=99.95311;GentrainScore=0.7423;HW=1.0
chr1                    629906  SNP1-555149     C       T       .       PASS                    CR=99.16279;GentrainScore=0.7029;HW=1.0
chr1                    634244  SNP1-559487     T       C       .       PASS                    CR=97.697975;GentrainScore=0.8070;HW=1.0
chr1                    753806  rs4000335       G       A       .       NOT_POLY_IN_1000G       CR=99.906586;GentrainScore=0.7934;HW=1.0
chr1                    788538  SNP1-713781     G       A       .       PASS                    CR=99.9051;GentrainScore=0.4541;HW=0.2935265
chr1                    794252  SNP1-719495     C       T       .       PASS                    CR=99.485985;GentrainScore=0.6870;HW=0.002817232
chr1                    817186  rs3094315       G       A       .       PASS                    CR=99.925735;GentrainScore=0.8141;HW=8.3064125E-11
chr1                    817341  rs3131972       A       G       .       PASS                    CR=99.859055;GentrainScore=0.8578;HW=1.7165839E-7
chr1                    818683  SNP1-743926     G       T       .       PASS                    CR=99.95258;GentrainScore=0.5893;HW=0.29127777
chr1                    821272  SNP1-746515     T       G       .       NOT_POLY_IN_1000G       CR=100.0;GentrainScore=0.6899;HW=1.0
chr1                    822311  SNP1-747554     T       C       .       PASS                    CR=99.85697;GentrainScore=0.5544;HW=0.04430997
chr1                    823439  SNP1-748682     C       T       .       NOT_POLY_IN_1000G       CR=99.95329;GentrainScore=0.5732;HW=1.0
chr1                    823656  SNP1-748899     G       A       .       PASS                    CR=92.56258;GentrainScore=0.6754;HW=0.037904903
.
.
.
chr2_KI270894v1_alt     191123  SNP2-90979776   C       T       .       NOT_POLY_IN_1000G       CR=99.95329;GentrainScore=0.6799;HW=1.0
chr2_KI270894v1_alt     191201  SNP2-90979698   A       G       .       NOT_POLY_IN_1000G       CR=99.67305;GentrainScore=0.7260;HW=1.0
chr2_KI270894v1_alt     195846  SNP2-90975053   T       C       .       PASS                    CR=99.953285;GentrainScore=0.7153;HW=1.0
chr2_KI270894v1_alt     199502  SNP2-90971397   G       A       .       NOT_POLY_IN_1000G       CR=99.85988;GentrainScore=0.7007;HW=1.0
chr2_KI270894v1_alt     202673  SNP2-90968226   G       A       .       PASS                    CR=99.44108;GentrainScore=0.7609;HW=1.0
chr2_KI270894v1_alt     204093  SNP2-90966806   G       A       .       NOT_POLY_IN_1000G       CR=99.71976;GentrainScore=0.6774;HW=1.0
chr2_KI270894v1_alt     206575  SNP2-90964324   G       A       .       PASS                    CR=99.25547;GentrainScore=0.6906;HW=1.0
chr2_KI270894v1_alt     208601  SNP2-90962298   T       A       .       PASS                    CR=99.53402;GentrainScore=0.7552;HW=1.0
chr17_KI270909v1_alt    224256  rs4796184       T       G       .       NOT_POLY_IN_1000G       CR=99.95329;GentrainScore=0.7755;HW=1.0
chr19_KI270938v1_alt    273807  SNP19-59494505  C       T       .       PASS                    CR=99.904854;GentrainScore=0.7221;HW=0.11059104
chr19_KI270938v1_alt    273962  SNP19-59494660  C       T       .       PASS                    CR=99.370766;GentrainScore=0.7625;HW=2.5549047E-5
chr19_KI270938v1_alt    276643  SNP19-59497341  T       C       .       PASS                    CR=89.78102;GentrainScore=0.5594;HW=4.737169E-5
chr19_KI270938v1_alt    277913  SNP19-59498611  C       T       .       PASS                    CR=90.28072;GentrainScore=0.7765;HW=4.5566884E-4
chr19_KI270938v1_alt    327468  SNP19-59548166  G       A       .       PASS                    CR=99.83753;GentrainScore=0.8322;HW=0.008843731
```

### 1000G indels

#### Sourced from the [GATK resource bundle](https://gatkforums.broadinstitute.org/gatk/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it)

```bash
## Not available?
```

```bash

```

### 1000G snps

#### Sourced from the [GATK resource bundle](https://gatkforums.broadinstitute.org/gatk/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it)

```bash
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz
```

```bash
#CHROM                  POS     ID              REF     ALT     QUAL            FILTER  INFO
chr1                    51479   rs116400033     T       A       11726.81        PASS    AC=229;AF=0.3253;AN=704;BaseQRankSum=-6.949;DB;DP=1570;Dels=0.00;FS=3.130;HRun=0;HaplotypeScore=0.1377;InbreedingCoeff=0.2907;MQ=34.37;MQ0=174;MQRankSum=1.476;QD=16.08;ReadPosRankSum=-0.202;SB=-4317.78;VQSLOD=5.1635;pop=EUR.admix
chr1                    55367   .               G       A       207.20          PASS    AC=2;AF=0.00117;AN=1714;BaseQRankSum=2.243;DP=4926;Dels=0.00;FS=3.005;HRun=0;HaplotypeScore=0.1382;InbreedingCoeff=-0.0188;MQ=45.57;MQ0=365;MQRankSum=0.185;QD=21.22;ReadPosRankSum=0.136;SB=-111.01;VQSLOD=6.3979;pop=ALL
chr1                    55388   .               C       T       95.61           PASS    AC=1;AF=0.00056;AN=1792;BaseQRankSum=-0.038;DP=5282;Dels=0.00;FS=0.000;HRun=2;HaplotypeScore=0.1980;InbreedingCoeff=-0.0278;MQ=48.19;MQ0=20;MQRankSum=-0.397;QD=18.13;ReadPosRankSum=-0.945;SB=-59.46;VQSLOD=5.7297;pop=ALL
chr1                    55852   .               G       C       503.17          PASS    AC=6;AF=0.0080;AN=748;BaseQRankSum=1.052;DP=1772;Dels=0.00;FS=1.623;HRun=0;HaplotypeScore=0.1542;InbreedingCoeff=-0.0547;MQ=39.86;MQ0=114;MQRankSum=1.243;QD=9.15;ReadPosRankSum=1.657;SB=-289.09;VQSLOD=5.8416;pop=EUR.admix
chr1                    61462   rs56992750      T       A       2671.41         PASS    AC=70;AF=0.04679;AN=1496;BaseQRankSum=-0.165;DB;DP=3256;Dels=0.00;FS=1.789;HRun=1;HaplotypeScore=0.1091;InbreedingCoeff=-0.0562;MQ=36.49;MQ0=395;MQRankSum=2.142;QD=6.50;ReadPosRankSum=0.789;SB=-1462.83;VQSLOD=6.2119;pop=ALL
chr1                    62157   rs10399597      G       A       285.13          PASS    AC=5;AF=0.00274;AN=1826;BaseQRankSum=3.510;DB;DP=5819;Dels=0.00;FS=0.802;HRun=0;HaplotypeScore=0.2357;InbreedingCoeff=-0.0660;MQ=39.50;MQ0=993;MQRankSum=2.274;QD=7.11;ReadPosRankSum=0.978;SB=-213.74;VQSLOD=5.1283;pop=ALL
chr1                    82609   .               C       G       1821.30         PASS    AC=53;AF=0.0716;AN=740;BaseQRankSum=-2.473;DP=1820;Dels=0.00;FS=4.066;HRun=1;HaplotypeScore=0.1207;InbreedingCoeff=0.1324;MQ=53.52;MQ0=320;MQRankSum=0.750;QD=9.59;ReadPosRankSum=-0.216;SB=-1153.72;VQSLOD=6.8223;pop=EUR.admix
chr1                    82734   rs4030331       T       C       12360.27        PASS    AC=316;AF=0.20654;AN=1530;BaseQRankSum=-7.774;DB;DP=3665;Dels=0.00;FS=3.305;HRun=0;HaplotypeScore=0.1481;InbreedingCoeff=0.0578;MQ=52.42;MQ0=272;MQRankSum=0.086;QD=8.82;ReadPosRankSum=2.068;SB=-3526.93;VQSLOD=6.9968;pop=ALL
chr1                    83084   .               T       A       29167.41        PASS    AC=726;AF=0.8462;AN=858;BaseQRankSum=3.752;DP=1533;Dels=0.00;FS=6.435;HRun=3;HaplotypeScore=0.0799;InbreedingCoeff=0.2112;MQ=46.84;MQ0=108;MQRankSum=1.586;QD=24.11;ReadPosRankSum=1.138;SB=-17136.27;VQSLOD=5.5460;pop=ALL
chr1                    83771   .               T       G       663.46          PASS    AC=22;AF=0.01391;AN=1582;BaseQRankSum=-4.181;DP=4131;Dels=0.00;FS=3.802;HRun=0;HaplotypeScore=0.1154;InbreedingCoeff=-0.0062;MQ=38.11;MQ0=947;MQRankSum=-0.593;QD=6.27;ReadPosRankSum=0.038;SB=-368.75;VQSLOD=5.6475;pop=ALL
.
.
.

chr19_KI270938v1_alt    278726  .               T       C       10308.33        PASS    AC=208;AF=0.15686;AN=1326;BaseQRankSum=-14.899;DP=2535;Dels=0.00;FS=2.612;HRun=2;HaplotypeScore=0.2625;InbreedingCoeff=0.2486;MQ=51.71;MQ0=238;MQRankSum=5.526;QD=17.38;ReadPosRankSum=2.087;SB=-3101.33;VQSLOD=6.4418;pop=ALL
chr19_KI270938v1_alt    282341  .               A       T       55.72           PASS    AC=1;AF=0.0021;AN=484;BaseQRankSum=1.500;DP=852;Dels=0.00;FS=2.440;HRun=0;HaplotypeScore=0.1413;InbreedingCoeff=-0.0404;MQ=84.29;MQ0=0;MQRankSum=-0.552;QD=16.38;ReadPosRankSum=1.718;SB=-30.44;VQSLOD=6.0772;pop=ALL
chr19_KI270938v1_alt    282516  .               T       A       79.95           PASS    AC=2;AF=0.0169;AN=118;BaseQRankSum=-3.135;DP=260;Dels=0.00;FS=5.123;HRun=0;HaplotypeScore=0.2289;InbreedingCoeff=0.1421;MQ=53.81;MQ0=11;MQRankSum=-1.180;QD=11.42;ReadPosRankSum=-0.830;SB=-28.57;VQSLOD=5.2958;pop=AMR.admix
chr19_KI270938v1_alt    327349  rs73938674      C       G,T     208.04          PASS    AC=4;AF=0.00207;AN=1934;BaseQRankSum=0.600;DB;DP=4841;Dels=0.00;FS=1.984;HRun=1;HaplotypeScore=0.3034;InbreedingCoeff=-0.0783;MQ=71.96;MQ0=116;MQRankSum=0.964;QD=4.62;ReadPosRankSum=1.798;SB=-137.57;VQSLOD=5.9874;pop=ALL
chr19_KI270938v1_alt    327468  rs73061008      G       A       44068.81        PASS    AC=574;AF=0.28027;AN=2048;BaseQRankSum=4.666;DB;DP=5439;Dels=0.00;FS=26.323;HRun=0;HaplotypeScore=0.2635;InbreedingCoeff=0.0888;MQ=130.46;MQ0=0;MQRankSum=11.224;QD=17.17;ReadPosRankSum=0.342;SB=-20820.53;VQSLOD=9.4649;pop=ALL
chr19_KI270938v1_alt    327664  rs10405172      C       T       49419.30        PASS    AC=1229;AF=0.81283;AN=1512;BaseQRankSum=15.379;DB;DP=2905;Dels=0.00;FS=52.475;HRun=1;HaplotypeScore=0.1735;InbreedingCoeff=0.1852;MQ=54.07;MQ0=115;MQRankSum=-2.803;QD=20.69;ReadPosRankSum=0.683;SB=-25146.93;VQSLOD=5.4197;pop=ALL
chr19_KI270938v1_alt    327688  rs57525193      C       T       497.98          PASS    AC=20;AF=0.01302;AN=1536;BaseQRankSum=5.755;DB;DP=2781;Dels=0.00;FS=2.019;HRun=0;HaplotypeScore=0.2344;InbreedingCoeff=-0.0767;MQ=55.74;MQ0=32;MQRankSum=3.445;QD=6.02;ReadPosRankSum=0.805;SB=-318.88;VQSLOD=7.2107;pop=ALL
chr19_KI270938v1_alt    327719  .               T       G       319.14          PASS    AC=23;AF=0.0431;AN=534;BaseQRankSum=-4.655;DP=1101;Dels=0.00;FS=5.942;HRun=2;HaplotypeScore=0.2028;InbreedingCoeff=-0.0033;MQ=50.03;MQ0=45;MQRankSum=-0.040;QD=4.66;ReadPosRankSum=1.291;SB=-189.57;VQSLOD=5.1866;pop=AFR.admix
chr19_KI270938v1_alt    327775  .               A       G       171.95          PASS    AC=11;AF=0.0281;AN=392;BaseQRankSum=-3.234;DP=1020;Dels=0.00;FS=0.000;HRun=0;HaplotypeScore=0.1918;InbreedingCoeff=-0.0171;MQ=39.01;MQ0=164;MQRankSum=-1.160;QD=4.97;ReadPosRankSum=-1.388;SB=-104.49;VQSLOD=5.2143;pop=AFR.admix
```

### Mills

#### Sourced from the [GATK resource bundle](https://gatkforums.broadinstitute.org/gatk/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it)

```bash
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
```

```bash
#CHROM                  POS     ID              REF     ALT             QUAL            FILTER  INFO
chr1                    55249   .               C       CTATGG          6160.83         PASS    set=Intersect1000GMinusBI
chr1                    87114   .               CT      C               666.86          PASS    set=Intersect1000GAll
chr1                    263835  .               CCT     C               1083.97         PASS    set=Intersect1000GMinusBI
chr1                    756187  .               CACAG   C               2260.77         PASS    set=Intersect1000GMinusBI
chr1                    779929  .               AG      A               108.86          PASS    set=Intersect1000GMinusOX
chr1                    786694  .               A       ATT             361.06          PASS    set=Intersect1000GMinusDI
chr1                    796754  .               T       TA              1722.44         PASS    set=Intersect1000GAll
chr1                    816927  .               AT      A               2106.07         PASS    set=Intersect1000GAll
chr1                    822796  .               CCT     C               436.88          PASS    set=Intersect1000GAll
chr1                    824064  rs78982110      ATATT   A               6659.29         PASS    set=Intersect1000GAll
chr1                    829704  .               TTTTG   T               1998.11         PASS    set=Intersect1000GMinusSI
chr1                    829775  rs71576582      TTGAC   T               1466.09         PASS    set=Intersect1000GAll
chr1                    833758  rs59306077      CAT     C               29441.73        PASS    set=Intersect1000GMinusSI
chr1                    836930  .               G       GTT             826.84          PASS    set=Intersect1000GAll
chr1                    839006  .               CAG     C               95.90           PASS    set=Intersect1000GAll
.
.
.
chr17_KI270909v1_alt    232448  1769495         C       CTAAT           .               PASS    set=MillsDoubleCenter
chr17_KI270909v1_alt    240559  1396155         T       TGG,TGGGG,TGTG  .               PASS    set=MillsTracesUnknown
chr17_KI270909v1_alt    244544  1571409         G       GC              14846.09        PASS    set=MillsAlleleMatch1000G-MillsDoubleCenter
chr17_KI270909v1_alt    256951  1990797         C       CGAAA           8875.86         PASS    set=MillsCenterUnknown-MillsAlleleMatch1000G-MillsTracesUnknown
chr17_KI270909v1_alt    259929  199091          A       ATTAT           .               PASS    set=MillsDoubleCenter
chr17_KI270909v1_alt    271425  1000179         A       AGT             .               PASS    set=MillsDoubleCenter
chr17_KI270909v1_alt    289794  1396523         CA      C,CAA,CAAA      .               PASS    set=MillsTracesUnknown
chr19_KI270938v1_alt    274597  .               AG      A               487.76          PASS    set=Intersect1000GMinusOX
chr19_KI270938v1_alt    275383  .               CCTT    C               4874.67         PASS    set=Intersect1000GMinusSI
chr19_KI270938v1_alt    276880  .               TTCTCCC T               1316.97         PASS    set=Intersect1000GMinusSI
chr19_KI270938v1_alt    535096  .               TTGTA   T               1979.57         PASS    set=Intersect1000GMinusSI
chr19_KI270938v1_alt    535122  .               CGT     C               5759            PASS    set=Intersect1000GMinusSI
chr19_KI270938v1_alt    535130  .               TGA     T               4101.79         PASS    set=Intersect1000GMinusSI
chr19_KI270938v1_alt    535158  .               TGA     T               3176.81         PASS    set=Intersect1000GMinusSI
chr19_KI270938v1_alt    535173  .               GTA     G               6850.34         PASS    set=Intersect1000GMinusSI
```

### Axiom exome plus

#### Sourced from the [GATK resource bundle](https://gatkforums.broadinstitute.org/gatk/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it)

```bash
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org:21/bundle/hg38/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
```

```bash
#CHROM                  POS     ID              REF     ALT     QUAL    FILTER  INFO
chr1                    69496   rs150690004     G       A       .       PASS    AC=34;AF=0.014;AN=2466;DB;NS=1233
chr1                    69761   .               A       T       .       PASS    AC=51;AF=0.021;AN=2462;NS=1231
chr1                    786344  .               A       AT      .       PASS    AC=27;AF=0.011;AN=2498;NS=1249
chr1                    827105  rs148989274     C       A       .       PASS    AC=467;AF=0.188;AN=2488;DB;NS=1244
chr1                    930165  .               G       A       .       PASS    AC=1;AF=4.010e-04;AN=2494;NS=1247
chr1                    930204  rs148711625     G       A       .       PASS    AC=21;AF=8.407e-03;AN=2498;DB;NS=1249
chr1                    930245  rs146327803     G       A       .       PASS    AC=1;AF=4.003e-04;AN=2498;DB;NS=1249
chr1                    930248  rs41285790      G       A       .       PASS    AC=9;AF=3.603e-03;AN=2498;DB;NS=1249
chr1                    930314  rs9988179       C       T       .       PASS    AC=147;AF=0.059;AN=2490;DB;NS=1245
chr1                    930320  rs116730894     C       T       .       PASS    AC=6;AF=2.408e-03;AN=2492;DB;NS=1246
chr1                    935779  rs143282473     G       A       .       PASS    AC=2;AF=8.026e-04;AN=2492;DB;NS=1246
chr1                    935849  rs149944086     G       C       .       PASS    AC=3;AF=1.201e-03;AN=2498;DB;NS=1249
.
.
.
chr7_KI270803v1_alt     453542  rs73742263      T       G       .       PASS    AC=2470;AF=0.990;AN=2494;DB;NS=1247
chr7_KI270803v1_alt     491373  rs361452        G       A       .       PASS    AC=2314;AF=0.928;AN=2494;DB;NS=1247
chr7_KI270803v1_alt     495615  rs361387        C       T       .       PASS    AC=11;AF=4.432e-03;AN=2482;DB;NS=1241
chr7_KI270803v1_alt     516048  rs71545367      C       A       .       PASS    AC=53;AF=0.021;AN=2486;DB;NS=1243
chr7_KI270803v1_alt     775980  rs57725689      T       C       .       PASS    AC=2391;AF=0.958;AN=2496;DB;NS=1248
chr8_KI270821v1_alt     553680  rs4361783       T       G       .       PASS    AC=1140;AF=0.458;AN=2490;DB;NS=1245
chr15_KI270850v1_alt    14657   .               C       T       .       PASS    AC=428;AF=0.173;AN=2468;NS=1234
chr22_KI270879v1_alt    270759  rs144686326     G       A       .       PASS    AC=1261;AF=0.506;AN=2492;DB;NS=1246
chr22_KI270879v1_alt    275948  .               C       T       .       PASS    AC=1;AF=4.010e-04;AN=2494;NS=1247
chr19_KI270938v1_alt    195749  rs1052983       C       T       .       PASS    AC=79;AF=0.032;AN=2470;DB;NS=1235
chr19_KI270938v1_alt    273648  rs142893391     G       C       .       PASS    AC=2;AF=8.097e-04;AN=2470;DB;NS=1235
chr19_KI270938v1_alt    274346  rs141580607     C       G       .       PASS    AC=1;AF=4.023e-04;AN=2486;DB;NS=1243
```

### CADD

#### Sourced from [CADD](https://cadd.gs.washington.edu/download)

```bash
aria2c -c -s 10 https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/whole_genome_SNVs.tsv.gz
```

```bash
#Chrom  Pos     Ref     Alt     RawScore        PHRED
1       10001   T       A       0.719644        9.191
1       10001   T       C       0.764391        9.588
1       10001   T       G       0.739988        9.367
1       10002   A       C       0.733459        9.310
1       10002   A       G       0.756198        9.512
1       10002   A       T       0.718562        9.182
1       10003   A       C       0.733963        9.314
1       10003   A       G       0.756702        9.517
1       10003   A       T       0.719066        9.186
1       10004   C       A       0.583042        8.052
1       10004   C       G       0.604244        8.231
1       10004   C       T       0.636603        8.500
1       10005   C       A       0.583634        8.057
1       10005   C       G       0.604836        8.236
1       10005   C       T       0.637195        8.505
1       10006   C       A       0.584746        8.066
.
.
.
Y       57217411        G       A       0.060657        2.348
Y       57217411        G       C       0.035815        2.091
Y       57217411        G       T       0.022050        1.956
Y       57217412        T       A       0.188344        3.911
Y       57217412        T       C       0.233091        4.487
Y       57217412        T       G       0.208688        4.175
Y       57217413        G       A       0.061304        2.355
Y       57217413        G       C       0.036462        2.097
Y       57217413        G       T       0.022697        1.962
Y       57217414        G       A       0.060826        2.350
Y       57217414        G       C       0.035984        2.093
Y       57217414        G       T       0.022219        1.958
Y       57217415        T       A       0.188583        3.914
Y       57217415        T       C       0.233330        4.490
Y       57217415        T       G       0.208926        4.178
```
