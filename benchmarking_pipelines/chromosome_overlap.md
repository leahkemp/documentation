# Chromosomes

Some chromosomes/chromosome patches are unique to the pipeline output vcf or the truth vcf (see). Variant calls on these regions could inflate the number of 'false-positive' and 'false-negative' calls.

This document describes the chromosomes for which variants were called during benchmarking for the pipeline output for ([human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) and [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline)) runs on the Genome In A Bottle (GIAB) sample [NA12878](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) (NIST7035 and NIST7086). This is compared against the chromosomes on which variants were called in the truth vcf.

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/bench1.0/
grep -v "#" ./vcf_annotation_pipeline/filtered/NIST7035_NIST_filtered.vcf | awk '{print $1}' | uniq
grep -v "#" ../known/project.NIST.hc.snps.indels.vcf | cut -f NIST7035 | awk '{print $1}' | uniq
```

- NIST7035

| pipeline output: NIST7035_NIST_filtered.vcf | truth vcf: project.NIST.hc.snps.indels.vcf.gz |
|---------------------------------------------|-----------------------------------------------|
| chrM                                        | chrM                                          |
| chr1                                        | chr1                                          |
| chr2                                        | chr2                                          |
| chr3                                        | chr3                                          |
| chr4                                        | chr4                                          |
| chr5                                        | chr5                                          |
| chr6                                        | chr6                                          |
| chr7                                        | chr7                                          |
| chr8                                        | chr8                                          |
| chr9                                        | chr9                                          |
| chr10                                       | chr10                                         |
| chr11                                       | chr11                                         |
| chr12                                       | chr12                                         |
| chr13                                       | chr13                                         |
| chr14                                       | chr14                                         |
| chr15                                       | chr15                                         |
| chr16                                       | chr16                                         |
| chr17                                       | chr17                                         |
| chr18                                       | chr18                                         |
| chr19                                       | chr19                                         |
| chr20                                       | chr20                                         |
| chr21                                       | chr21                                         |
| chr22                                       | chr22                                         |
| chrX                                        | chrX                                          |
| chrY                                        | chrY                                          |
| chr1_gl000191_random                        | chr1_gl000191_random                          |
| chr1_gl000192_random                        | chr1_gl000192_random                          |
| chr4_ctg9_hap1                              | chr4_ctg9_hap1                                |
| chr4_gl000193_random                        | chr4_gl000193_random                          |
| chr4_gl000194_random                        | chr4_gl000194_random                          |
| -                                           | chr6_apd_hap1                                 |
| chr6_cox_hap2                               | chr6_cox_hap2                                 |
| chr6_dbb_hap3                               | chr6_dbb_hap3                                 |
| chr6_mann_hap4                              | chr6_mann_hap4                                |
| chr6_mcf_hap5                               | chr6_mcf_hap5                                 |
| chr6_qbl_hap6                               | chr6_qbl_hap6                                 |
| chr6_ssto_hap7                              | chr6_ssto_hap7                                |
| chr7_gl000195_random                        | chr7_gl000195_random                          |
| chr8_gl000197_random                        | -                                             |
| chr9_gl000198_random                        | chr9_gl000198_random                          |
| chr9_gl000199_random                        | chr9_gl000199_random                          |
| chr17_ctg5_hap1                             | chr17_ctg5_hap1                               |
| chr17_gl000203_random                       | chr17_gl000203_random                         |
| chr17_gl000204_random                       | chr17_gl000204_random                         |
| chr17_gl000205_random                       | chr17_gl000205_random                         |
| chr17_gl000206_random                       | -                                             |
| chr18_gl000207_random                       | chr18_gl000207_random                         |
| chr19_gl000208_random                       | chr19_gl000208_random                         |
| chr19_gl000209_random                       | chr19_gl000209_random                         |
| chrUn_gl000211                              | chrUn_gl000211                                |
| chrUn_gl000212                              | chrUn_gl000212                                |
| chrUn_gl000213                              | chrUn_gl000213                                |
| chrUn_gl000214                              | chrUn_gl000214                                |
| chrUn_gl000215                              | chrUn_gl000215                                |
| chrUn_gl000216                              | chrUn_gl000216                                |
| chrUn_gl000217                              | chrUn_gl000217                                |
| chrUn_gl000218                              | chrUn_gl000218                                |
| chrUn_gl000219                              | chrUn_gl000219                                |
| chrUn_gl000220                              | chrUn_gl000220                                |
| chrUn_gl000221                              | chrUn_gl000221                                |
| chrUn_gl000222                              | chrUn_gl000222                                |
| chrUn_gl000223                              | chrUn_gl000223                                |
| chrUn_gl000224                              | chrUn_gl000224                                |
| chrUn_gl000225                              | chrUn_gl000225                                |
| chrUn_gl000226                              | chrUn_gl000226                                |
| chrUn_gl000227                              | -                                             |
| chrUn_gl000228                              | chrUn_gl000228                                |
| chrUn_gl000229                              | chrUn_gl000229                                |
| chrUn_gl000230                              | chrUn_gl000230                                |
| chrUn_gl000231                              | chrUn_gl000231                                |
| chrUn_gl000232                              | chrUn_gl000232                                |
| chrUn_gl000233                              | chrUn_gl000233                                |
| chrUn_gl000234                              | chrUn_gl000234                                |
| chrUn_gl000235                              | chrUn_gl000235                                |
| chrUn_gl000236                              | -                                             |
| chrUn_gl000237                              | chrUn_gl000237                                |
| chrUn_gl000238                              | chrUn_gl000238                                |
| chrUn_gl000239                              | chrUn_gl000239                                |
| chrUn_gl000240                              | chrUn_gl000240                                |
| chrUn_gl000241                              | chrUn_gl000241                                |
| chrUn_gl000242                              | chrUn_gl000242                                |
| chrUn_gl000243                              | chrUn_gl000243                                |
| chrUn_gl000244                              | -                                             |
| chrUn_gl000245                              | -                                             |
| chrUn_gl000246                              | -                                             |
| chrUn_gl000247                              | chrUn_gl000247                                |
| -                                           | chrUn_gl000248                                |

- NIST7086

```bash
cd /store/lkemp/exome_project/benchmarking/NA12878_exome/bench1.0/
grep -v "#" ./vcf_annotation_pipeline/filtered/NIST7086_NIST_filtered.vcf | awk '{print $1}' | uniq
grep -v "#" ../known/project.NIST.hc.snps.indels.vcf | cut -f NIST7086 | awk '{print $1}' | uniq
```
