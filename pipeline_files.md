# Description of the files output by two genomic pipelines

Created: 2020-04-22 17:35:19
Last modified: 2020/04/23 15:16:17

This document show how the data is manipulated/changed as the input exome data (fastq files) data goes through both the [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) and the [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline), resulting in the final annotated vcf file. Use this to help understand what happens to the data after key stages in the pipelines. See [this file](https://samtools.github.io/hts-specs/VCFv4.2.pdf) for an overview of vcf files and descriptions on what the data fields represent.

---

This example was run on the Genome In A Bottle (GIAB) sample 'NA12878' (see [here](https://github.com/ESR-NZ/ESR-Parabricks) for more information on how this data was prepared).

Both pipelines were run against the GRCh38 reference genome (Homo_sapiens_assembly38.fasta) and known indels (Homo_sapiens_assembly38.known_indels.vcf.gz") from the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle))

---

## Table of contents

- [Description of the files output by two genomic pipelines](#description-of-the-files-output-by-two-genomic-pipelines)
  - [Table of contents](#table-of-contents)
  - [human_genomics_pipeline](#humangenomicspipeline)
    - [Input fastq files](#input-fastq-files)
    - [After trimming](#after-trimming)
    - [After mapping, sorting, marking duplicates, assigning reads and recalibration](#after-mapping-sorting-marking-duplicates-assigning-reads-and-recalibration)
    - [After haplotype calling](#after-haplotype-calling)
  - [vcf_annotation_pipeline](#vcfannotationpipeline)
    - [After genotyping](#after-genotyping)
    - [After modelling and applying recalibration](#after-modelling-and-applying-recalibration)
    - [After annotation with SnpSift](#after-annotation-with-snpsift)
    - [After annotation with vep](#after-annotation-with-vep)
    - [After annotation with genmod](#after-annotation-with-genmod)

## human_genomics_pipeline

### Input fastq files

*Input file:* NA12878_NIST_R1.fastq.gz

```bash
zcat NA12878_NIST_R1.fastq.gz | head -n 10
```

```bash
@HWI-D00119:50:H7AP8ADXX:1:1101:1213:2058 1:N:0:TAAGGCGA
ACTCCAGCCTGGGCAACAGAGCAAGGCTCGGTCTCCCAAAAAAAAAAAAAAAAAAAAAAAATTGGAACTCATTTAAAAACACTTATGAAGAGTTCATTTCT
+
@@@D?BD?A>CBDCED;EFGF;@B3?::8))0)8?B>B@FGCFEEBC######################################################
@HWI-D00119:50:H7AP8ADXX:1:1101:1242:2178 1:N:0:TAAGGCGA
ACATAGTGGTCTGTCTTCTGTTTATTACAGTACCTGTAATAATTCTTGATGCTGTCCAGACTTGGAATATTGAATGCATCTAATTTAAAAAAAAAAAAAAA
+
@@@DDBAAB=AACEEEFEEAHFEEECHAH>CH4A<+:CF<FFE<<DF@BFB9?4?<*?9D;DCFCEEFIII<=FFF4=FDEECDFIFFEECBCBBBBBB>5
@HWI-D00119:50:H7AP8ADXX:1:1101:1326:2059 1:N:0:TAAGGCGA
ATCTCTAGATCTATCTATCGTCTACCATTATCCATCGACCTATCTGTCTCTATTATATACACCTATCATCTACCTAATCTATCTGTCAATGATTTTCTGTC
[...]
```

*Input file:* NA12878_NIST_R2.fastq.gz

```bash
zcat NA12878_NIST_R2.fastq.gz | head -n 10
```

```
@HWI-D00119:50:H7AP8ADXX:1:1101:1213:2058 2:N:0:TAAGGCGA
NNNNNAGGGNGGACNCNNTGANAGGAAGAAATGACNTCTTCCTAAGTGTTTTTAAATTACTTCCATGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
+
#####-28@#4==@#2##32=#2:9==?<????=?#07=??=??>?>?>??<;?????>?@@=>?####################################
@HWI-D00119:50:H7AP8ADXX:1:1101:1242:2178 2:N:0:TAAGGCGA
CCCTTAGTCATTAAAGAAAATGTTTGTATTCTTTTTTTTAATGTGGATGAATGTTTGCTTTTTTTTTTTTTTTAAATTAGATGCATTCAATATTCCAAGTC
+
???DDFDFDHHHDIJJJJGIGIIIJIIGIIEIIIJJJJJHHIIIHIDGIICGHJGJHBGHGIIHFDDDDDDDB?>3:@:+:44((4@:AA:3>C:@CD::3
@HWI-D00119:50:H7AP8ADXX:1:1101:1326:2059 2:N:0:TAAGGCGA
CNNNNATACATAGANANNTAANAGAGATGAGAGTTGGATAGAGAAGTAGGTAGAAAGATAGATAGATANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
[...]
```

### After trimming

*Rule:* trim_galore_pe (trim_galore)

*Output file:* NA12878_NIST_R1_val_1.fq.gz

```bash
zcat NA12878_NIST_R1_val_1.fq.gz | head -n 10
```

```
@HWI-D00119:50:H7AP8ADXX:1:1101:1213:2058 1:N:0:TAAGGCGA
ACTCCAGCCTGGGCAACAGAGCAAGGCTCGGTCTCCCAAAAAAAAA
+
@@@D?BD?A>CBDCED;EFGF;@B3?::8))0)8?B>B@FGCFEEB
@HWI-D00119:50:H7AP8ADXX:1:1101:1242:2178 1:N:0:TAAGGCGA
ACATAGTGGTCTGTCTTCTGTTTATTACAGTACCTGTAATAATTCTTGATGCTGTCCAGACTTGGAATATTGAATGCATCTAATTTAAAAAAAAAAAAAA
+
@@@DDBAAB=AACEEEFEEAHFEEECHAH>CH4A<+:CF<FFE<<DF@BFB9?4?<*?9D;DCFCEEFIII<=FFF4=FDEECDFIFFEECBCBBBBBB>
@HWI-D00119:50:H7AP8ADXX:1:1101:1326:2059 1:N:0:TAAGGCGA
ATCTCTAGATCTATCTATCGTCTACCATTATCCATCGACCTATCTGTCTCTATTATATACACCTATCATCTACCTAATCTATCTGTCAATGATTTTCTGTC
[...]
```

*Output file:* NA12878_NIST_R1_val_2.fq.gz

```bash
zcat NA12878_NIST_R2_val_2.fq.gz | head -n 10
```

```
@HWI-D00119:50:H7AP8ADXX:1:1101:1213:2058 2:N:0:TAAGGCGA
NNNNNAGGGNGGACNCNNTGANAGGAAGAAATGACNTCTTCCTAAGTGTTTTTAAATTACTTCC
+
#####-28@#4==@#2##32=#2:9==?<????=?#07=??=??>?>?>??<;?????>?@@=>
@HWI-D00119:50:H7AP8ADXX:1:1101:1242:2178 2:N:0:TAAGGCGA
CCCTTAGTCATTAAAGAAAATGTTTGTATTCTTTTTTTTAATGTGGATGAATGTTTGCTTTTTTTTTTTTTTTAAATTAGATGCATTCAATATTCCAAGT
+
???DDFDFDHHHDIJJJJGIGIIIJIIGIIEIIIJJJJJHHIIIHIDGIICGHJGJHBGHGIIHFDDDDDDDB?>3:@:+:44((4@:AA:3>C:@CD::
@HWI-D00119:50:H7AP8ADXX:1:1101:1326:2059 2:N:0:TAAGGCGA
CNNNNATACATAGANANNTAANAGAGATGAGAGTTGGATAGAGAAGTAGGTAGAAAGATAGATAGAT
[...]
```

### After mapping, sorting, marking duplicates, assigning reads and recalibration

*Rule:* bwa_map (bwa mem)
*Rule:* sambamba_sort (sambamba sort)
*Rule:* gatk4_readgroup_add (gatk AddOrReplaceReadGroups)
*Rule:* gatk4_recal (gatk ApplyBQSR)

*Output file:* NA12878_NIST_bwa_recal.bam

```bash
samtools view NA12878_NIST_bwa_recal.bam | head -n 10
```

```
HWI-D00119:50:H7AP8ADXX:2:2106:14134:70021      99      chr1    10010   0       3M1I84M         =       10058   128   CCCTTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC                 <@@?>=?>@@@=>=???<?>@@@=?>@@@=?>@@@=?>@@@=?>@?@=?>@??><?A=<<-<?64<7?@:>==?>:><:>>>?=7><@               MC:Z:80M          MD:Z:87    RG:Z:4  NM:i:1  AS:i:84  XS:i:84
HWI-D00119:50:H7AP8ADXX:1:2208:5166:2843        99      chr1    10014   0       85M             =       10064   100   AACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA                    ;?68?>;?=?>@=><>>@<>>@?@<>>@@?=>???@=>?A@=>?>???<@?@?@>?>@@@<?=?=>;?????<?>=?=<<=>??;                  MC:Z:50M          MD:Z:85    RG:Z:4  NM:i:0  AS:i:85  XS:i:85
HWI-D00119:50:H7AP8ADXX:2:2103:8894:14988       99      chr1    10016   0       15M1I6M1I72M    =       10256   304   CCCTAACCCTAACCCTTAACCCTTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC          <@>?=?>@@?=>>>??><?>@???=?>>@@=?>@??>??@@@=?>@@@>?>??@=???@@<@?AA><?9>>?=??@?@<><:<=;@??<>6<?>>        MC:Z:10S14M1I50M  MD:Z:93    RG:Z:4  NM:i:2  AS:i:79  XS:i:87
HWI-D00119:50:H7AP8ADXX:1:1112:11788:8995       163     chr1    10046   0       63M             =       10343   397   CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC                                          =@@?=?>@@@=>>???=?=@@@=?>@@@=?=??@>?=A??<?=@??>??>A@;@>@A><?=>>                                        MC:Z:10M1I90M     MD:Z:63    RG:Z:4  NM:i:0  AS:i:63  XS:i:63
HWI-D00119:50:H7AP8ADXX:1:1103:12425:88782      99      chr1    10052   0       40M61S          =       10111   97    CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTGTCTCTTATACACATCTCCGAGCCCACGAGACTAAGGCGAATCTCGTATGCCGTCTTCTGC    <@@?=?>@@@=>>???<>>@@@=?>@@@=?>@@@=?>@@@A=@@@@?===>>>>=@@@@9>AAAA??;>A=??<?AA@;>?<@@A;<<<B@@;<A?>A?B@  MC:Z:63S38M       MD:Z:40    RG:Z:4  NM:i:0  AS:i:40  XS:i:41
HWI-D00119:50:H7AP8ADXX:2:2106:14134:70021      147     chr1    10058   5       80M             =       10010   -128  CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAACCCTAACCCTAACCCTAACCCTAACCC                         AA@<<<AB?<9=A@@><<BA@<?<BBB;>=@BA>>=B@@<=;BA>:</AAA><AA@<?<AAA<><@@A=><AAA==;AA<                       MC:Z:3M1I84M      MD:Z:80    RG:Z:4  NM:i:0  AS:i:80  XS:i:73
HWI-D00119:50:H7AP8ADXX:1:2208:5166:2843        147     chr1    10064   0       50M             =       10014   -100  CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAACCC                                                       BA=:>>ABA>?;BAB>>9A@A<74A@?5?=A@@=><AAA;=<?@?=;AA<                                                     MC:Z:85M          MD:Z:50    RG:Z:4  NM:i:0  AS:i:50  XS:i:50
HWI-D00119:50:H7AP8ADXX:2:1201:2571:99329       83      chr1    10088   0       21M1I60M        =       10135   -34   CCCTAACCCTAACCCTAACCCGAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCC                       @>::<;@=:-;7@=@>=<@>8096@?B,<=B@@<86>BA<=<@@@+=<@@B<?<AAA<>=A@@@=><@@A==;AAA;><AA<                     MC:Z:43M1I14M     MD:Z:81    RG:Z:4  NM:i:1  AS:i:74  XS:i:77
HWI-D00119:50:H7AP8ADXX:2:2214:6410:98782       147     chr1    10100   0       11M1I38M        =       10117   -32   CCCTAACCCAAGCCCTACCCCTAACCCTAACCCTAACCCTAACCCTAACC                                                       ABA<:;682(=(AA@*<&AAA5>;AA@=><A@A<><@A@;>;@@@;<2?<                                                     MC:Z:44M          MD:Z:16A32 RG:Z:4  NM:i:2  AS:i:37  XS:i:38
HWI-D00119:50:H7AP8ADXX:1:1103:12425:88782      147     chr1    10111   0       63S38M          =       10052   -97   CGGCGACCACCGAGATCTACACCTCTCTATTCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACC    ;@@;B<AB<A;A@A=>A<=B<AA>A>A<<?>;??:AAB@@9>>A@@=>>>>===?@@@@=A@@AAA=?<AAA=><AAA=?<AA@=><@AA;>=AAA==;A<  MC:Z:40M61S       MD:Z:38    RG:Z:4  NM:i:0  AS:i:38  XS:i:40
[...]
```

### After haplotype calling

*Rule:* gatk4_HaplotypeCaller (gatk HaplotypeCaller)

*Output file:* NA12878_NIST.raw.snps.indels.AS.g.vcf

```bash
cat NA12878_NIST.raw.snps.indels.AS.g.vcf | grep -v '^##' | head -n 10
```

```bash
#CHROM  POS     ID      REF     ALT             QUAL    FILTER  INFO            FORMAT                  20
chr1    1       .       N       <NON_REF>       .       .       END=10353       GT:DP:GQ:MIN_DP:PL      0/0:0:0:0:0,0,0
chr1    10354   .       C       <NON_REF>       .       .       END=10357       GT:DP:GQ:MIN_DP:PL      0/0:1:3:1:0,3,30
chr1    10358   .       A       <NON_REF>       .       .       END=10378       GT:DP:GQ:MIN_DP:PL      0/0:2:6:2:0,6,38
chr1    10379   .       C       <NON_REF>       .       .       END=10379       GT:DP:GQ:MIN_DP:PL      0/0:0:0:0:0,0,0
chr1    10380   .       C       <NON_REF>       .       .       END=10383       GT:DP:GQ:MIN_DP:PL      0/0:2:6:2:0,6,47
chr1    10384   .       C       <NON_REF>       .       .       END=10431       GT:DP:GQ:MIN_DP:PL      0/0:1:3:1:0,3,14
chr1    10432   .       A       <NON_REF>       .       .       END=12915       GT:DP:GQ:MIN_DP:PL      0/0:0:0:0:0,0,0
chr1    12916   .       T       <NON_REF>       .       .       END=13005       GT:DP:GQ:MIN_DP:PL      0/0:1:3:1:0,3,21
chr1    13006   .       G       <NON_REF>       .       .       END=13057       GT:DP:GQ:MIN_DP:PL      0/0:0:0:0:0,0,0
[...]
```

## vcf_annotation_pipeline

### After genotyping

*Rule:* gatk4_GenotypeGVCFs (gatk GenotypeGVCFs)

*Output file:* NA12878_NIST.genotype.vcf

```bash
cat NA12878_NIST.genotype.vcf | grep -v '^##' | head -n 10
```

```bash
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO                                                                                                                                                                                     FORMAT                           20
chr1    13684   .       C       T       65.84   .       AC=2;AF=1.00;AN=2;AS_QD=22.00;DP=3;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=21.67;QD=21.95;SOR=2.833                                                                             GT:AD:DP:GQ:PL                   1/1:0,3:3:9:79,9,0
chr1    13813   .       T       G       148.14  .       AC=2;AF=1.00;AN=2;AS_QD=25.36;DP=4;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=21.75;QD=28.73;SOR=3.258                                                                              GT:AD:DP:GQ:PGT:PID:PL:PS        1|1:0,4:4:12:1|1:13813_T_G:162,12,0:13813
chr1    13838   .       C       T       121.84  .       AC=2;AF=1.00;AN=2;AS_QD=30.97;DP=3;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=21.67;QD=27.24;SOR=2.833                                                                             GT:AD:DP:GQ:PGT:PID:PL:PS        1|1:0,3:3:9:1|1:13813_T_G:135,9,0:13813
chr1    14397   .       CTGT    C       702.60  .       AC=1;AF=0.500;AN=2;AS_QD=11.16;BaseQRankSum=0.877;DP=63;ExcessHet=3.0103;FS=4.983;MLEAC=1;MLEAF=0.500;MQ=22.44;MQRankSum=-4.600e-02;QD=11.15;ReadPosRankSum=2.69;SOR=0.050               GT:AD:DP:GQ:PL                   0/1:43,20:63:99:710,0,1745
chr1    14522   .       G       A       140.64  .       AC=1;AF=0.500;AN=2;AS_QD=9.40;BaseQRankSum=1.62;DP=15;ExcessHet=3.0103;FS=3.310;MLEAC=1;MLEAF=0.500;MQ=22.83;MQRankSum=-2.200e-02;QD=9.38;ReadPosRankSum=0.340;SOR=2.546                 GT:AD:DP:GQ:PL                   0/1:8,7:15:99:148,0,175
chr1    14542   .       A       G       228.64  .       AC=1;AF=0.500;AN=2;AS_QD=16.36;BaseQRankSum=2.91;DP=14;ExcessHet=3.0103;FS=5.441;MLEAC=1;MLEAF=0.500;MQ=22.89;MQRankSum=1.11;QD=16.33;ReadPosRankSum=0.573;SOR=3.442                     GT:AD:DP:GQ:PL                   0/1:4,10:14:66:236,0,66
chr1    14574   .       A       G       199.64  .       AC=1;AF=0.500;AN=2;AS_QD=18.18;BaseQRankSum=-2.362e+00;DP=15;ExcessHet=3.0103;FS=7.404;MLEAC=1;MLEAF=0.500;MQ=22.30;MQRankSum=-8.750e-01;QD=18.15;ReadPosRankSum=-7.780e-01;SOR=4.615    GT:AD:DP:GQ:PL                   0/1:2,9:11:21:207,0,21
chr1    14590   .       G       A       95.64   .       AC=1;AF=0.500;AN=2;AS_QD=9.60;BaseQRankSum=2.20;DP=10;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=21.82;MQRankSum=-2.510e-01;QD=9.56;ReadPosRankSum=0.921;SOR=1.911                 GT:AD:DP:GQ:PL                   0/1:5,5:10:99:103,0,103
chr1    14599   .       T       A       236.97  .       AC=2;AF=1.00;AN=2;AS_QD=28.20;DP=6;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=21.17;QD=25.00;SOR=3.912                                                                              GT:AD:DP:GQ:PGT:PID:PL:PS        1|1:0,6:6:18:1|1:14599_T_A:251,18,0:14599
[...]
```

### After modelling and applying recalibration

*Rule:* gatk4_VariantRecalibrator_indel (gatk *VariantRecalibrator)
*Rule:* gatk4_VariantRecalibrator_snp (gatk *VariantRecalibrator)
*Rule:* gatk4_VQSR_indel (gatk ApplyVQSR)
*Rule:* gatk4_VQSR_snp (gatk ApplyVQSR)

*Output file:* NA12878_NIST.recal.indels

```bash
cat NA12878_NIST.recal.indels | grep -v '^##' | head -n 10
```

```bash
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
chr1    14397   .       N       <VQSR>  .       .       END=14400;NEGATIVE_TRAIN_SITE;VQSLOD=-2.0599;culprit=ReadPosRankSum
chr1    15903   .       N       <VQSR>  .       .       END=15903;NEGATIVE_TRAIN_SITE;VQSLOD=-2.1970;culprit=ReadPosRankSum
chr1    129010  .       N       <VQSR>  .       .       END=129013;VQSLOD=1.7222;culprit=SOR
chr1    136644  .       N       <VQSR>  .       .       END=136645;VQSLOD=-102.7024;culprit=FS
chr1    187264  .       N       <VQSR>  .       .       END=187264;VQSLOD=2.3150;culprit=MQRankSum
chr1    188420  .       N       <VQSR>  .       .       END=188420;NEGATIVE_TRAIN_SITE;VQSLOD=-0.9928;culprit=MQRankSum
chr1    267777  .       N       <VQSR>  .       .       END=267777;NEGATIVE_TRAIN_SITE;VQSLOD=-4.8233;culprit=ReadPosRankSum
chr1    514446  .       N       <VQSR>  .       .       END=514447;NEGATIVE_TRAIN_SITE;VQSLOD=-4.5079;culprit=MQRankSum
chr1    826577  .       N       <VQSR>  .       .       END=826577;VQSLOD=-2.8863;culprit=DP
[...]
```

*Output file:* NA12878_NIST.recal.snps

```bash
cat NA12878_NIST.recal.snps | grep -v '^##' | head -n 10
```

```bash
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
chr1    13684   .       N       <VQSR>  .       .       END=13684;VQSLOD=-12.0689;culprit=MQ
chr1    13813   .       N       <VQSR>  .       .       END=13813;VQSLOD=-9.4555;culprit=MQ
chr1    13838   .       N       <VQSR>  .       .       END=13838;VQSLOD=-12.0964;culprit=MQ
chr1    14522   .       N       <VQSR>  .       .       END=14522;VQSLOD=-11.6279;culprit=MQ
chr1    14542   .       N       <VQSR>  .       .       END=14542;VQSLOD=-12.6650;culprit=MQ
chr1    14574   .       N       <VQSR>  .       .       END=14574;VQSLOD=-14.2163;culprit=MQ
chr1    14590   .       N       <VQSR>  .       .       END=14590;VQSLOD=-11.0787;culprit=MQ
chr1    14599   .       N       <VQSR>  .       .       END=14599;VQSLOD=-10.5777;culprit=MQ
chr1    14604   .       N       <VQSR>  .       .       END=14604;VQSLOD=-11.1508;culprit=MQ
[...]
```

*Output file:* NA12878_NIST.vqsr.recal.vcf

```bash
cat NA12878_NIST.vqsr.recal.vcf | grep -v '^##' | head -n 10
```

```bash
#CHROM  POS     ID      REF     ALT     QUAL    FILTER                          INFO                                                                                                                                                                                                                                     FORMAT                          20
chr1    13684   .       C       T       65.84   VQSRTrancheSNP99.95to100.00     AC=2;AF=1.00;AN=2;AS_QD=22.00;DP=3;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=21.67;QD=21.95;SOR=2.833;VQSLOD=-1.207e+01;culprit=MQ                                                                                                GT:AD:DP:GQ:PL                  1/1:0,3:3:9:79,9,
chr1    13813   .       T       G       148.14  VQSRTrancheSNP99.90to99.95      AC=2;AF=1.00;AN=2;AS_QD=25.36;DP=4;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=21.75;QD=28.73;SOR=3.258;VQSLOD=-9.456e+00;    culprit=MQ                                                                                             GT:AD:DP:GQ:PGT:PID:PL:PS       1|1:0,4:4:12:1|1:13813_T_G:162,12,0:13813
chr1    13838   .       C       T       121.84  VQSRTrancheSNP99.95to100.00     AC=2;AF=1.00;AN=2;AS_QD=30.97;DP=3;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=21.67;QD=27.24;SOR=2.833;VQSLOD=-1.210e+01;culprit=MQ                                                                                                GT:AD:DP:GQ:PGT:PID:PL:PS       1|1:0,3:3:9:1|1:13813_T_G:135,9,0:13813
chr1    14397   .       CTGT    C       702.60  PASS                            AC=1;AF=0.500;AN=2;AS_QD=11.16;BaseQRankSum=0.877;DP=63;ExcessHet=3.0103;FS=4.983;MLEAC=1;MLEAF=0.500;MQ=22.44;MQRankSum=-4.600e-02;NEGATIVE_TRAIN_SITE;QD=11.15;ReadPosRankSum=2.69;SOR=0.050;VQSLOD=-2.060e+00;culprit=ReadPosRankSum  GT:AD:DP:GQ:PL                  0/1:43,20:63:99:710,0,1745
chr1    14522   .       G       A       140.64  VQSRTrancheSNP99.95to100.00     AC=1;AF=0.500;AN=2;AS_QD=9.40;BaseQRankSum=1.62;DP=15;ExcessHet=3.0103;FS=3.310;MLEAC=1;MLEAF=0.500;MQ=22.83;MQRankSum=-2.200e-02;QD=9.38;ReadPosRankSum=0.340;SOR=2.546;VQSLOD=-1.163e+01;culprit=MQ                                    GT:AD:DP:GQ:PL                  0/1:8,7:15:99:148,0,175
chr1    14542   .       A       G       228.64  VQSRTrancheSNP99.95to100.00     AC=1;AF=0.500;AN=2;AS_QD=16.36;BaseQRankSum=2.91;DP=14;ExcessHet=3.0103;FS=5.441;MLEAC=1;MLEAF=0.500;MQ=22.89;MQRankSum=1.11;QD=16.33;ReadPosRankSum=0.573;SOR=3.442;VQSLOD=-1.267e+01;culprit=MQ                                        GT:AD:DP:GQ:PL                  0/1:4,10:14:66:236,0,66
chr1    14574   .       A       G       199.64  VQSRTrancheSNP99.95to100.00     AC=1;AF=0.500;AN=2;AS_QD=18.18;BaseQRankSum=-2.362e+00;DP=15;ExcessHet=3.0103;FS=7.404;MLEAC=1;MLEAF=0.500;MQ=22.30;MQRankSum=-8.750e-01;QD=18.15;ReadPosRankSum=-7.780e-01;SOR=4.615;VQSLOD=-1.422e+01;culprit=MQ                       GT:AD:DP:GQ:PL                  0/1:2,9:11:21:207,0,21
chr1    14590   .       G       A       95.64   VQSRTrancheSNP99.95to100.00     AC=1;AF=0.500;AN=2;AS_QD=9.60;BaseQRankSum=2.20;DP=10;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=21.82;MQRankSum=-2.510e-01;QD=9.56;ReadPosRankSum=0.921;SOR=1.911;VQSLOD=-1.108e+01;culprit=MQ                                    GT:AD:DP:GQ:PL                  0/1:5,5:10:99:103,0,103
chr1    14599   .       T       A       236.97  VQSRTrancheSNP99.95to100.00     AC=2;AF=1.00;AN=2;AS_QD=28.20;DP=6;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=21.17;QD=25.00;SOR=3.912;VQSLOD=-1.058e+01;culprit=MQ                                                                                                 GT:AD:DP:GQ:PGT:PID:PL:PS       1|1:0,6:6:18:1|1:14599_T_A:251,18,0:14599
[...]
```

### After annotation with SnpSift

*Rule:* SnpSift (SnpSift dbnsfp)

*Output file:* NA12878_NIST.vqsr.recal.dbnsfp.vcf

```bash
cat NA12878_NIST.vqsr.recal.dbnsfp.vcf | grep -v '^##' | head -n 10
```

```bash

#CHROM  POS     ID      REF     ALT     QUAL    FILTER                          INFO                                                                                                                                                                                                                                     FORMAT                          20
chr1    13684   .       C       T       65.84   VQSRTrancheSNP99.95to100.00     AC=2;AF=1.00;AN=2;AS_QD=22.00;DP=3;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=21.67;QD=21.95;SOR=2.833;VQSLOD=-1.207e+01;culprit=MQ                                                                                                GT:AD:DP:GQ:PL                  1/1:0,3:3:9:79,9,0
chr1    13813   .       T       G       148.14  VQSRTrancheSNP99.90to99.95      AC=2;AF=1.00;AN=2;AS_QD=25.36;DP=4;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=21.75;QD=28.73;SOR=3.258;VQSLOD=-9.456e+00;culprit=MQ                                                                                                 GT:AD:DP:GQ:PGT:PID:PL:PS       1|1:0,4:4:12:1|1:13813_T_G:162,12,0:13813
chr1    13838   .       C       T       121.84  VQSRTrancheSNP99.95to100.00     AC=2;AF=1.00;AN=2;AS_QD=30.97;DP=3;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=21.67;QD=27.24;SOR=2.833;VQSLOD=-1.210e+01;culprit=MQ                                                                                                GT:AD:DP:GQ:PGT:PID:PL:PS       1|1:0,3:3:9:1|1:13813_T_G:135,9,0:13813
chr1    14397   .       CTGT    C       702.6   PASS                            AC=1;AF=0.500;AN=2;AS_QD=11.16;BaseQRankSum=0.877;DP=63;ExcessHet=3.0103;FS=4.983;MLEAC=1;MLEAF=0.500;MQ=22.44;MQRankSum=-4.600e-02;NEGATIVE_TRAIN_SITE;QD=11.15;ReadPosRankSum=2.69;SOR=0.050;VQSLOD=-2.060e+00;culprit=ReadPosRankSum  GT:AD:DP:GQ:PL                  0/1:43,20:63:99:710,0,1745
chr1    14522   .       G       A       140.64  VQSRTrancheSNP99.95to100.00     AC=1;AF=0.500;AN=2;AS_QD=9.40;BaseQRankSum=1.62;DP=15;ExcessHet=3.0103;FS=3.310;MLEAC=1;MLEAF=0.500;MQ=22.83;MQRankSum=-2.200e-02;QD=9.38;ReadPosRankSum=0.340;SOR=2.546;VQSLOD=-1.163e+01;culprit=MQ                                    GT:AD:DP:GQ:PL                  0/1:8,7:15:99:148,0,175
chr1    14542   .       A       G       228.64  VQSRTrancheSNP99.95to100.00     AC=1;AF=0.500;AN=2;AS_QD=16.36;BaseQRankSum=2.91;DP=14;ExcessHet=3.0103;FS=5.441;MLEAC=1;MLEAF=0.500;MQ=22.89;MQRankSum=1.11;QD=16.33;ReadPosRankSum=0.573;SOR=3.442;VQSLOD=-1.267e+01;culprit=MQ                                        GT:AD:DP:GQ:PL                  0/1:4,10:14:66:236,0,66
chr1    14574   .       A       G       199.64  VQSRTrancheSNP99.95to100.00     AC=1;AF=0.500;AN=2;AS_QD=18.18;BaseQRankSum=-2.362e+00;DP=15;ExcessHet=3.0103;FS=7.404;MLEAC=1;MLEAF=0.500;MQ=22.30;MQRankSum=-8.750e-01;QD=18.15;ReadPosRankSum=-7.780e-01;SOR=4.615;VQSLOD=-1.422e+01;culprit=MQ                       GT:AD:DP:GQ:PL                  0/1:2,9:11:21:207,0,21
chr1    14590   .       G       A       95.64   VQSRTrancheSNP99.95to100.00     AC=1;AF=0.500;AN=2;AS_QD=9.60;BaseQRankSum=2.20;DP=10;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=21.82;MQRankSum=-2.510e-01;QD=9.56;ReadPosRankSum=0.921;SOR=1.911;VQSLOD=-1.108e+01;culprit=M                                     GT:AD:DP:GQ:PL                  0/1:5,5:10:99:103,0,103
chr1    14599   .       T       A       236.97  VQSRTrancheSNP99.95to100.00     AC=2;AF=1.00;AN=2;AS_QD=28.20;DP=6;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=21.17;QD=25.00;SOR=3.912;VQSLOD=-1.058e+01;culprit=MQ                                                                                                 GT:AD:DP:GQ:PGT:PID:PL:PS       1|1:0,6:6:18:1|1:14599_T_A:251,18,0:14599
[...]
```

### After annotation with vep

*Rule:* VEP (vep)

*Output file:* NA12878_NIST.vqsr.recal.dbnsfp.vep.vcf

```bash
cat NA12878_NIST.vqsr.recal.dbnsfp.vep.vcf | grep -v '^##' | head -n 10
```

```bash
#CHROM  POS     ID      REF     ALT     QUAL    FILTER                          INFO                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  FORMAT                      20
chr1    13684   .       C       T       65.84   VQSRTrancheSNP99.95to100.00     AC=2;AF=1.00;AN=2;AS_QD=22.00;DP=3;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=21.67;QD=21.95;SOR=2.833;VQSLOD=-1.207e+01;culprit=MQ;CSQ=T|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene||||||||||rs71260404|14|1||SNV|HGNC|HGNC:37102||||||||||||||||||||||||||||||||||||||||||,T|non_coding_transcript_exon_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript|3/3||ENST00000456328.2:n.932C>T||932|||||rs71260404||1||SNV|HGNC|HGNC:37102|YES||1|||||||||||||||||||||||||||||||||||||||,T|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene||||||||||rs71260404|720|-1||SNV|HGNC|HGNC:38034|YES|||||||||||||||||||||||||||||||||||||||||,T|downstream_gene_variant|MODIFIER|MIR6859-1|ENSG00000278267|Transcript|ENST00000619216|miRNA||||||||||rs71260404|3685|-1||SNV|HGNC|HGNC:50039|YES|||||||||||||||||||||||||||||||||||||||||                                                                                                                                                                                   GT:AD:DP:GQ:PL              1/1:0,3:3:9:79,9,0
chr1    13813   .       T       G       148.14  VQSRTrancheSNP99.90to99.95      AC=2;AF=1.00;AN=2;AS_QD=25.36;DP=4;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=21.75;QD=28.73;SOR=3.258;VQSLOD=-9.456e+00;culprit=MQ;CSQ=G|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene||||||||||rs1213979446|143|1||SNV|HGNC|HGNC:37102||||||||||||||||||||||||||||||||||||||||||,G|non_coding_transcript_exon_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript|3/3||ENST00000456328.2:n.1061T>G||1061|||||rs1213979446||1||SNV|HGNC|HGNC:37102|YES||1|||||||||||||||||||||||||||||||||||||||,G|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene||||||||||rs1213979446|591|-1||SNV|HGNC|HGNC:38034|YES|||||||||||||||||||||||||||||||||||||||||,G|downstream_gene_variant|MODIFIER|MIR6859-1|ENSG00000278267|Transcript|ENST00000619216|miRNA||||||||||rs1213979446|3556|-1||SNV|HGNC|HGNC:50039|YES|||||||||||||||||||||||||||||||||||||||||                                                                                                                                                                         GT:AD:DP:GQ:PGT:PID:PL:PS   1|1:0,4:4:12:1|1:13813_T_G:162,12,0:13813
chr1    13838   .       C       T       121.84  VQSRTrancheSNP99.95to100.00     AC=2;AF=1.00;AN=2;AS_QD=30.97;DP=3;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=21.67;QD=27.24;SOR=2.833;VQSLOD=-1.210e+01;culprit=MQ;CSQ=T|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene||||||||||rs28428499|168|1||SNV|HGNC|HGNC:37102||||||||||||||||||||||||||||||||||||||||||,T|non_coding_transcript_exon_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript|3/3||ENST00000456328.2:n.1086C>T||1086|||||rs28428499||1||SNV|HGNC|HGNC:37102|YES||1|||||||||||||||||||||||||||||||||||||||,T|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene||||||||||rs28428499|566|-1||SNV|HGNC|HGNC:38034|YES|||||||||||||||||||||||||||||||||||||||||,T|downstream_gene_variant|MODIFIER|MIR6859-1|ENSG00000278267|Transcript|ENST00000619216|miRNA||||||||||rs28428499|3531|-1||SNV|HGNC|HGNC:50039|YES|||||||||||||||||||||||||||||||||||||||||                                                                                                                                                                                GT:AD:DP:GQ:PGT:PID:PL:PS   1|1:0,3:3:9:1|1:13813_T_G:135,9,0:13813
chr1    14397   .       CTGT    C       702.6   PASS                            AC=1;AF=0.500;AN=2;AS_QD=11.16;BaseQRankSum=0.877;DP=63;ExcessHet=3.0103;FS=4.983;MLEAC=1;MLEAF=0.500;MQ=22.44;MQRankSum=-4.600e-02;NEGATIVE_TRAIN_SITE;QD=11.15;ReadPosRankSum=2.69;SOR=0.050;VQSLOD=-2.060e+00;culprit=ReadPosRankSum;CSQ=-|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene||||||||||rs756427959|728|1||deletion|HGNC|HGNC:37102||||||||||||||||||||||||||||||||||||||||||,-|non_coding_transcript_exon_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript|3/3||ENST00000456328.2:n.1648_1650del||1646-1648|||||rs756427959||1||deletion|HGNC|HGNC:37102|YES||1||||||||||||2|||||||||||||||||||||||||||,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene||||||||||rs756427959|4|-1||deletion|HGNC|HGNC:38034|YES|||||||||||||||||||||||||||||||||||||||||,-|downstream_gene_variant|MODIFIER|MIR6859-1|ENSG00000278267|Transcript|ENST00000619216|miRNA||||||||||rs756427959|2969|-1||deletion|HGNC|HGNC:50039|YES|||||||||||||||||||||||||||||||||||||||||                                                 GT:AD:DP:GQ:PL              0/1:43,20:63:99:710,0,1745
chr1    14522   .       G       A       140.64  VQSRTrancheSNP99.95to100.00     AC=1;AF=0.500;AN=2;AS_QD=9.40;BaseQRankSum=1.62;DP=15;ExcessHet=3.0103;FS=3.310;MLEAC=1;MLEAF=0.500;MQ=22.83;MQRankSum=-2.200e-02;QD=9.38;ReadPosRankSum=0.340;SOR=2.546;VQSLOD=-1.163e+01;culprit=MQ;CSQ=A|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene||||||||||rs1441808061|852|1||SNV|HGNC|HGNC:37102||||||||||||||||||||||||||||||||||||||||||,A|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript||||||||||rs1441808061|113|1||SNV|HGNC|HGNC:37102|YES||1|||||||||||||||||||||||||||||||||||||||,A|intron_variant&non_coding_transcript_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene||10/10|ENST00000488147.1:n.1254-21C>T|||||||rs1441808061||-1||SNV|HGNC|HGNC:38034|YES|||||||||||||||||||||||||||||||||||||||||,A|downstream_gene_variant|MODIFIER|MIR6859-1|ENSG00000278267|Transcript|ENST00000619216|miRNA||||||||||rs1441808061|2847|-1||SNV|HGNC|HGNC:50039|YES|||||||||||||||||||||||||||||||||||||||||                                                                                                 GT:AD:DP:GQ:PL              0/1:8,7:15:99:148,0,175
chr1    14542   .       A       G       228.64  VQSRTrancheSNP99.95to100.00     AC=1;AF=0.500;AN=2;AS_QD=16.36;BaseQRankSum=2.91;DP=14;ExcessHet=3.0103;FS=5.441;MLEAC=1;MLEAF=0.500;MQ=22.89;MQRankSum=1.11;QD=16.33;ReadPosRankSum=0.573;SOR=3.442;VQSLOD=-1.267e+01;culprit=MQ;CSQ=G|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene||||||||||rs1045833|872|1||SNV|HGNC|HGNC:37102||||||||||||||||||||||||||||||||||||||||||,G|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript||||||||||rs1045833|133|1||SNV|HGNC|HGNC:37102|YES||1|||||||||||||||||||||||||||||||||||||||,G|intron_variant&non_coding_transcript_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene||10/10|ENST00000488147.1:n.1254-41T>C|||||||rs1045833||-1||SNV|HGNC|HGNC:38034|YES|||||||||||||||||||||||||||||||||||||||||,G|downstream_gene_variant|MODIFIER|MIR6859-1|ENSG00000278267|Transcript|ENST00000619216|miRNA||||||||||rs1045833|2827|-1||SNV|HGNC|HGNC:50039|YES|||||||||||||||||||||||||||||||||||||||||                                                                                                                 GT:AD:DP:GQ:PL              0/1:4,10:14:66:236,0,66
chr1    14574   .       A       G       199.64  VQSRTrancheSNP99.95to100.00     AC=1;AF=0.500;AN=2;AS_QD=18.18;BaseQRankSum=-2.362e+00;DP=15;ExcessHet=3.0103;FS=7.404;MLEAC=1;MLEAF=0.500;MQ=22.30;MQRankSum=-8.750e-01;QD=18.15;ReadPosRankSum=-7.780e-01;SOR=4.615;VQSLOD=-1.422e+01;culprit=MQ;CSQ=G|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene||||||||||rs28503599|904|1||SNV|HGNC|HGNC:37102||||||||||||||||||||||||||||||||||||||||||,G|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript||||||||||rs28503599|165|1||SNV|HGNC|HGNC:37102|YES||1|||||||||||||||||||||||||||||||||||||||,G|intron_variant&non_coding_transcript_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene||10/10|ENST00000488147.1:n.1254-73T>C|||||||rs28503599||-1||SNV|HGNC|HGNC:38034|YES|||||||||||||||||||||||||||||||||||||||||,G|downstream_gene_variant|MODIFIER|MIR6859-1|ENSG00000278267|Transcript|ENST00000619216|miRNA||||||||||rs28503599|2795|-1||SNV|HGNC|HGNC:50039|YES|||||||||||||||||||||||||||||||||||||||||                                                                                            GT:AD:DP:GQ:PL              0/1:2,9:11:21:207,0,21
chr1    14590   .       G       A       95.64   VQSRTrancheSNP99.95to100.00     AC=1;AF=0.500;AN=2;AS_QD=9.60;BaseQRankSum=2.20;DP=10;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=21.82;MQRankSum=-2.510e-01;QD=9.56;ReadPosRankSum=0.921;SOR=1.911;VQSLOD=-1.108e+01;culprit=MQ;CSQ=A|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene||||||||||rs707679|920|1||SNV|HGNC|HGNC:37102||||||||||||||||||||||||||||||||||||||||||,A|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript||||||||||rs707679|181|1||SNV|HGNC|HGNC:37102|YES||1|||||||||||||||||||||||||||||||||||||||,A|intron_variant&non_coding_transcript_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene||10/10|ENST00000488147.1:n.1254-89C>T|||||||rs707679||-1||SNV|HGNC|HGNC:38034|YES|||||||||||||||||||||||||||||||||||||||||,A|downstream_gene_variant|MODIFIER|MIR6859-1|ENSG00000278267|Transcript|ENST00000619216|miRNA||||||||||rs707679|2779|-1||SNV|HGNC|HGNC:50039|YES|||||||||||||||||||||||||||||||||||||||||                                                                                                                 GT:AD:DP:GQ:PL              0/1:5,5:10:99:103,0,103
chr1    14599   .       T       A       236.97  VQSRTrancheSNP99.95to100.00     AC=2;AF=1.00;AN=2;AS_QD=28.20;DP=6;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=21.17;QD=25.00;SOR=3.912;VQSLOD=-1.058e+01;culprit=MQ;CSQ=A|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene||||||||||rs707680|929|1||SNV|HGNC|HGNC:37102||||||||||||||||0.1476|0.121|0.1758|0.0893|0.161|0.2096||||||||||||0.2096|SAS||||||||,A|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript||||||||||rs707680|190|1||SNV|HGNC|HGNC:37102|YES||1|||||||||||||0.1476|0.121|0.1758|0.0893|0.161|0.2096||||||||||||0.2096|SAS||||||||,A|intron_variant&non_coding_transcript_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene||10/10|ENST00000488147.1:n.1254-98A>T|||||||rs707680||-1||SNV|HGNC|HGNC:38034|YES|||||||||||||||0.1476|0.121|0.1758|0.0893|0.161|0.2096||||||||||||0.2096|SAS||||||||,A|downstream_gene_variant|MODIFIER|MIR6859-1|ENSG00000278267|Transcript|ENST00000619216|miRNA||||||||||rs707680|2770|-1||SNV|HGNC|HGNC:50039|YES|||||||||||||||0.1476|0.121|0.1758|0.0893|0.161|0.2096||||||||||||0.2096|SAS||||||||  GT:AD:DP:GQ:PGT:PID:PL:PS   1|1:0,6:6:18:1|1:14599_T_A:251,18,0:14599
[...]
```

### After annotation with genmod

*Rule:* GENMOD (genmod annotate)

*Output file:* NA12878_NIST.vqsr.recal.dbnsfp.vep.genmod.vcf

```bash
cat NA12878_NIST.vqsr.recal.dbnsfp.vep.genmod.vcf | grep -v '^##' | head -n 10
```

```bash
#CHROM  POS     ID      REF     ALT     QUAL    FILTER                          INFO                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  FORMAT                       20
chr1    13684   .       C       T       65.84   VQSRTrancheSNP99.95to100.00     AC=2;AF=1.00;AN=2;AS_QD=22.00;DP=3;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=21.67;QD=21.95;SOR=2.833;VQSLOD=-1.207e+01;culprit=MQ;CSQ=T|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene||||||||||rs71260404|14|1||SNV|HGNC|HGNC:37102||||||||||||||||||||||||||||||||||||||||||,T|non_coding_transcript_exon_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript|3/3||ENST00000456328.2:n.932C>T||932|||||rs71260404||1||SNV|HGNC|HGNC:37102|YES||1|||||||||||||||||||||||||||||||||||||||,T|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene||||||||||rs71260404|720|-1||SNV|HGNC|HGNC:38034|YES|||||||||||||||||||||||||||||||||||||||||,T|downstream_gene_variant|MODIFIER|MIR6859-1|ENSG00000278267|Transcript|ENST00000619216|miRNA||||||||||rs71260404|3685|-1||SNV|HGNC|HGNC:50039|YES|||||||||||||||||||||||||||||||||||||||||;Annotation=DDX11L1,WASH7P;CADD=2.058                                                                                                                                                                              GT:AD:DP:GQ:PL               1/1:0,3:3:9:79,9,0
chr1    13813   .       T       G       148.14  VQSRTrancheSNP99.90to99.95      AC=2;AF=1.00;AN=2;AS_QD=25.36;DP=4;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=21.75;QD=28.73;SOR=3.258;VQSLOD=-9.456e+00;culprit=MQ;CSQ=G|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene||||||||||rs1213979446|143|1||SNV|HGNC|HGNC:37102||||||||||||||||||||||||||||||||||||||||||,G|non_coding_transcript_exon_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript|3/3||ENST00000456328.2:n.1061T>G||1061|||||rs1213979446||1||SNV|HGNC|HGNC:37102|YES||1|||||||||||||||||||||||||||||||||||||||,G|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene||||||||||rs1213979446|591|-1||SNV|HGNC|HGNC:38034|YES|||||||||||||||||||||||||||||||||||||||||,G|downstream_gene_variant|MODIFIER|MIR6859-1|ENSG00000278267|Transcript|ENST00000619216|miRNA||||||||||rs1213979446|3556|-1||SNV|HGNC|HGNC:50039|YES|||||||||||||||||||||||||||||||||||||||||;Annotation=DDX11L1,WASH7P;CADD=5.197                                                                                                                                                                    GT:AD:DP:GQ:PGT:PID:PL:PS    1|1:0,4:4:12:1|1:13813_T_G:162,12,0:13813
chr1    13838   .       C       T       121.84  VQSRTrancheSNP99.95to100.00     AC=2;AF=1.00;AN=2;AS_QD=30.97;DP=3;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=21.67;QD=27.24;SOR=2.833;VQSLOD=-1.210e+01;culprit=MQ;CSQ=T|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene||||||||||rs28428499|168|1||SNV|HGNC|HGNC:37102||||||||||||||||||||||||||||||||||||||||||,T|non_coding_transcript_exon_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript|3/3||ENST00000456328.2:n.1086C>T||1086|||||rs28428499||1||SNV|HGNC|HGNC:37102|YES||1|||||||||||||||||||||||||||||||||||||||,T|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene||||||||||rs28428499|566|-1||SNV|HGNC|HGNC:38034|YES|||||||||||||||||||||||||||||||||||||||||,T|downstream_gene_variant|MODIFIER|MIR6859-1|ENSG00000278267|Transcript|ENST00000619216|miRNA||||||||||rs28428499|3531|-1||SNV|HGNC|HGNC:50039|YES|||||||||||||||||||||||||||||||||||||||||;Annotation=DDX11L1,WASH7P;CADD=3.372                                                                                                                                                                           GT:AD:DP:GQ:PGT:PID:PL:PS    1|1:0,3:3:9:1|1:13813_T_G:135,9,0:13813
chr1    14397   .       CTGT    C       702.6   PASS                            AC=1;AF=0.500;AN=2;AS_QD=11.16;BaseQRankSum=0.877;DP=63;ExcessHet=3.0103;FS=4.983;MLEAC=1;MLEAF=0.500;MQ=22.44;MQRankSum=-4.600e-02;NEGATIVE_TRAIN_SITE;QD=11.15;ReadPosRankSum=2.69;SOR=0.050;VQSLOD=-2.060e+00;culprit=ReadPosRankSum;CSQ=-|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene||||||||||rs756427959|728|1||deletion|HGNC|HGNC:37102||||||||||||||||||||||||||||||||||||||||||,-|non_coding_transcript_exon_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript|3/3||ENST00000456328.2:n.1648_1650del||1646-1648|||||rs756427959||1||deletion|HGNC|HGNC:37102|YES||1||||||||||||2|||||||||||||||||||||||||||,-|downstream_gene_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene||||||||||rs756427959|4|-1||deletion|HGNC|HGNC:38034|YES|||||||||||||||||||||||||||||||||||||||||,-|downstream_gene_variant|MODIFIER|MIR6859-1|ENSG00000278267|Transcript|ENST00000619216|miRNA||||||||||rs756427959|2969|-1||deletion|HGNC|HGNC:50039|YES|||||||||||||||||||||||||||||||||||||||||;Annotation=DDX11L1,WASH7P                                                       GT:AD:DP:GQ:PL               0/1:43,20:63:99:710,0,1745
chr1    14522   .       G       A       140.64  VQSRTrancheSNP99.95to100.00     AC=1;AF=0.500;AN=2;AS_QD=9.40;BaseQRankSum=1.62;DP=15;ExcessHet=3.0103;FS=3.310;MLEAC=1;MLEAF=0.500;MQ=22.83;MQRankSum=-2.200e-02;QD=9.38;ReadPosRankSum=0.340;SOR=2.546;VQSLOD=-1.163e+01;culprit=MQ;CSQ=A|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene||||||||||rs1441808061|852|1||SNV|HGNC|HGNC:37102||||||||||||||||||||||||||||||||||||||||||,A|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript||||||||||rs1441808061|113|1||SNV|HGNC|HGNC:37102|YES||1|||||||||||||||||||||||||||||||||||||||,A|intron_variant&non_coding_transcript_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene||10/10|ENST00000488147.1:n.1254-21C>T|||||||rs1441808061||-1||SNV|HGNC|HGNC:38034|YES|||||||||||||||||||||||||||||||||||||||||,A|downstream_gene_variant|MODIFIER|MIR6859-1|ENSG00000278267|Transcript|ENST00000619216|miRNA||||||||||rs1441808061|2847|-1||SNV|HGNC|HGNC:50039|YES|||||||||||||||||||||||||||||||||||||||||;Annotation=WASH7P;CADD=0.310                                                                                                    GT:AD:DP:GQ:PL               0/1:8,7:15:99:148,0,175
chr1    14542   .       A       G       228.64  VQSRTrancheSNP99.95to100.00     AC=1;AF=0.500;AN=2;AS_QD=16.36;BaseQRankSum=2.91;DP=14;ExcessHet=3.0103;FS=5.441;MLEAC=1;MLEAF=0.500;MQ=22.89;MQRankSum=1.11;QD=16.33;ReadPosRankSum=0.573;SOR=3.442;VQSLOD=-1.267e+01;culprit=MQ;CSQ=G|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene||||||||||rs1045833|872|1||SNV|HGNC|HGNC:37102||||||||||||||||||||||||||||||||||||||||||,G|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript||||||||||rs1045833|133|1||SNV|HGNC|HGNC:37102|YES||1|||||||||||||||||||||||||||||||||||||||,G|intron_variant&non_coding_transcript_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene||10/10|ENST00000488147.1:n.1254-41T>C|||||||rs1045833||-1||SNV|HGNC|HGNC:38034|YES|||||||||||||||||||||||||||||||||||||||||,G|downstream_gene_variant|MODIFIER|MIR6859-1|ENSG00000278267|Transcript|ENST00000619216|miRNA||||||||||rs1045833|2827|-1||SNV|HGNC|HGNC:50039|YES|||||||||||||||||||||||||||||||||||||||||;Annotation=WASH7P;CADD=2.989                                                                                                                    GT:AD:DP:GQ:PL               0/1:4,10:14:66:236,0,66
chr1    14574   .       A       G       199.64  VQSRTrancheSNP99.95to100.00     AC=1;AF=0.500;AN=2;AS_QD=18.18;BaseQRankSum=-2.362e+00;DP=15;ExcessHet=3.0103;FS=7.404;MLEAC=1;MLEAF=0.500;MQ=22.30;MQRankSum=-8.750e-01;QD=18.15;ReadPosRankSum=-7.780e-01;SOR=4.615;VQSLOD=-1.422e+01;culprit=MQ;CSQ=G|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene||||||||||rs28503599|904|1||SNV|HGNC|HGNC:37102||||||||||||||||||||||||||||||||||||||||||,G|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript||||||||||rs28503599|165|1||SNV|HGNC|HGNC:37102|YES||1|||||||||||||||||||||||||||||||||||||||,G|intron_variant&non_coding_transcript_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene||10/10|ENST00000488147.1:n.1254-73T>C|||||||rs28503599||-1||SNV|HGNC|HGNC:38034|YES|||||||||||||||||||||||||||||||||||||||||,G|downstream_gene_variant|MODIFIER|MIR6859-1|ENSG00000278267|Transcript|ENST00000619216|miRNA||||||||||rs28503599|2795|-1||SNV|HGNC|HGNC:50039|YES|||||||||||||||||||||||||||||||||||||||||;Annotation=WASH7P;CADD=5.831                                                                                               GT:AD:DP:GQ:PL               0/1:2,9:11:21:207,0,21
chr1    14590   .       G       A       95.64   VQSRTrancheSNP99.95to100.00     AC=1;AF=0.500;AN=2;AS_QD=9.60;BaseQRankSum=2.20;DP=10;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=21.82;MQRankSum=-2.510e-01;QD=9.56;ReadPosRankSum=0.921;SOR=1.911;VQSLOD=-1.108e+01;culprit=MQ;CSQ=A|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene||||||||||rs707679|920|1||SNV|HGNC|HGNC:37102||||||||||||||||||||||||||||||||||||||||||,A|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript||||||||||rs707679|181|1||SNV|HGNC|HGNC:37102|YES||1|||||||||||||||||||||||||||||||||||||||,A|intron_variant&non_coding_transcript_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene||10/10|ENST00000488147.1:n.1254-89C>T|||||||rs707679||-1||SNV|HGNC|HGNC:38034|YES|||||||||||||||||||||||||||||||||||||||||,A|downstream_gene_variant|MODIFIER|MIR6859-1|ENSG00000278267|Transcript|ENST00000619216|miRNA||||||||||rs707679|2779|-1||SNV|HGNC|HGNC:50039|YES|||||||||||||||||||||||||||||||||||||||||;Annotation=WASH7P;CADD=1.988                                                                                                                    GT:AD:DP:GQ:PL               0/1:5,5:10:99:103,0,103
chr1    14599   .       T       A       236.97  VQSRTrancheSNP99.95to100.00     AC=2;AF=1.00;AN=2;AS_QD=28.20;DP=6;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=21.17;QD=25.00;SOR=3.912;VQSLOD=-1.058e+01;culprit=MQ;CSQ=A|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000450305|transcribed_unprocessed_pseudogene||||||||||rs707680|929|1||SNV|HGNC|HGNC:37102||||||||||||||||0.1476|0.121|0.1758|0.0893|0.161|0.2096||||||||||||0.2096|SAS||||||||,A|downstream_gene_variant|MODIFIER|DDX11L1|ENSG00000223972|Transcript|ENST00000456328|processed_transcript||||||||||rs707680|190|1||SNV|HGNC|HGNC:37102|YES||1|||||||||||||0.1476|0.121|0.1758|0.0893|0.161|0.2096||||||||||||0.2096|SAS||||||||,A|intron_variant&non_coding_transcript_variant|MODIFIER|WASH7P|ENSG00000227232|Transcript|ENST00000488147|unprocessed_pseudogene||10/10|ENST00000488147.1:n.1254-98A>T|||||||rs707680||-1||SNV|HGNC|HGNC:38034|YES|||||||||||||||0.1476|0.121|0.1758|0.0893|0.161|0.2096||||||||||||0.2096|SAS||||||||,A|downstream_gene_variant|MODIFIER|MIR6859-1|ENSG00000278267|Transcript|ENST00000619216|miRNA||||||||||rs707680|2770|-1||SNV|HGNC|HGNC:50039|YES|||||||||||||||0.1476|0.121|0.1758|0.0893|0.161|0.2096||||||||||||0.2096|SAS||||||||;Annotation=WASH7P;CADD=4.933     GT:AD:DP:GQ:PGT:PID:PL:PS    1|1:0,6:6:18:1|1:14599_T_A:251,18,0:14599
[...]
```
