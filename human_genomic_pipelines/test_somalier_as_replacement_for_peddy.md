# Containerising madeline2

Created: 2020/08/31 11:36:14
Last modified: 2020/08/31 16:00:41

- **Aim:** Try/test [somalier](https://github.com/brentp/somalier) and see if the outputs can be integrated with [scout](https://github.com/Clinical-Genomics/scout)
- **Prerequisite software:** [singularity 2.5.2](https://singularity.lbl.gov/)
- **OS:** Ubuntu 16.04 (Wintermute - research server)

In our whfole genome sequencing (WGS) workflow (fastq to vcf to preparing data to be ingested into scout) have been using [peddy](https://github.com/brentp/peddy) to evaluate/check the sex/familial relationships for the joint vcf's of our cohort samples. The files that peddy outputs are used and presented in [scout](https://github.com/Clinical-Genomics/scout) to show whether the genomic data supports the familial relationships and sex's defined in the pedigree file. We have had some inconsistency in the results of these checks with peddy depending on how the data is subset (eg. exome capture regions only and whether non-pass variants removed from the vcf). [somalier](https://github.com/brentp/somalier) has been released which ["is a more scalable, faster, replacement for peddy that uses some of the same methods as peddy along with some new ones"](https://github.com/brentp/peddy#fast-pedigreevcf-qc). It also allows you to do these checks based on bam's instead of vcf's. So we want to try it out and see if it's outputs can still be integrated into scout!
## Table of contents

- [Containerising madeline2](#containerising-madeline2)
  - [Table of contents](#table-of-contents)
  - [Try somalier](#try-somalier)

## Try somalier

A docker image for somalier is hosted [here](https://hub.docker.com/r/brentp/somalier)

```bash
singularity pull docker://brentp/somalier:latest

BIND_PATHS="/store/lkemp/manual_pipeline_run/"
CONTAINER="/store/lkemp/manual_pipeline_run/WGS_b37/test/leahkemp-madeline2_container-master-latest.simg"

singularity exec \
-B ${BIND_PATHS} \
${CONTAINER} \
somalier extract \
-d extracted/ \
--sites sites.vcf.gz \
-f /data/human/g1k_v37_decoy.fa \
.bam
```