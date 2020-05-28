# Benchmarking results

## Table of contents

- [Benchmarking results](#benchmarking-results)
  - [Table of contents](#table-of-contents)
  - [Known vcf](#known-vcf)
  - [intra_truth_comparison](#intra_truth_comparison)
    - [Run parameters/settings](#run-parameterssettings)
    - [Results](#results)
        - [Compared with bedops intersect](#compared-with-bedops-intersect)
        - [Compared with hap.py + RTG tools](#compared-with-happy--rtg-tools)
  - [bench1.0](#bench10)
    - [Run parameters/settings](#run-parameterssettings-1)
    - [Results](#results-1)
      - [human_genomics_pipeline + minimal vcf_annotation_pipeline](#human_genomics_pipeline--minimal-vcf_annotation_pipeline)
        - [Compared with bedops intersect](#compared-with-bedops-intersect-1)
        - [Compared with hap.py + RTG tools](#compared-with-happy--rtg-tools-1)
      - [parabricks germline pipeline](#parabricks-germline-pipeline)
        - [Compared with bedops intersect](#compared-with-bedops-intersect-2)
        - [Compared with hap.py + RTG tools](#compared-with-happy--rtg-tools-2)
  - [bench1.1](#bench11)
    - [Run parameters/settings](#run-parameterssettings-2)
    - [Results](#results-2)
      - [human_genomics_pipeline + minimal vcf_annotation_pipeline](#human_genomics_pipeline--minimal-vcf_annotation_pipeline-1)
        - [Compared with bedops intersect](#compared-with-bedops-intersect-3)
        - [Compared with hap.py + RTG tools](#compared-with-happy--rtg-tools-3)
  - [bench1.2](#bench12)
    - [Run parameters/settings](#run-parameterssettings-3)
    - [Results](#results-3)
      - [human_genomics_pipeline + minimal vcf_annotation_pipeline](#human_genomics_pipeline--minimal-vcf_annotation_pipeline-2)
        - [Compared with bedops intersect](#compared-with-bedops-intersect-4)
        - [Compared with hap.py + RTG tools](#compared-with-happy--rtg-tools-4)
      - [parabricks germline pipeline?](#parabricks-germline-pipeline-1)
        - [Compared with bedops intersect](#compared-with-bedops-intersect-5)
        - [Compared with hap.py + RTG tools](#compared-with-happy--rtg-tools-5)
  - [bench1.3](#bench13)
    - [Run parameters/settings](#run-parameterssettings-4)
    - [Results](#results-4)
      - [human_genomics_pipeline + minimal vcf_annotation_pipeline](#human_genomics_pipeline--minimal-vcf_annotation_pipeline-3)
        - [Compared with bedops intersect](#compared-with-bedops-intersect-6)
        - [Compared with hap.py + RTG tools](#compared-with-happy--rtg-tools-6)

## Known vcf

Parameters used to create the Genome In A Bottle (GIAB) sample [NA12878](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) (both NIST7035 and NIST7086) known vcf file:

- **Inputs:**
  - Reference genome: ucsc.hg19.fasta
  - dbSNP database: dbsnp_137.hg19.vcf

- **Description:**
  - NIST7035: (/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7035.vcf)
    - Total number of variants in NIST7035: 336,003
    - There are 99,141 variants found in NIST7035 that are not found in NIST7086
    - 29.51% of variants in NIST7035 are unique to NIST7035

  - NIST7086 (/store/lkemp/publicData/exomes/NA12878_exome/project.NIST.hc.snps.indels.NIST7086.vcf)
    - Total number of variants in NIST7086: 317,524
    - There are 80,661 variants found in NIST7086 that are not found in NIST7035
    - 25.40% of variants in NIST7086 are unique to NIST7086

## intra_truth_comparison

### Run parameters/settings

- **Aim:** Compare the two 'samples' of the Genome In A Bottle (GIAB) sample [NA12878](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) (NIST7035 and NIST7086)

### Results

Results dir: /store/lkemp/exome_project/benchmarking/NA12878_exome/intra_truth_comparison/ (wintermute)

##### Compared with bedops intersect

| file                                                                                   | count   | type                                               |
|----------------------------------------------------------------------------------------|---------|----------------------------------------------------|
| common_project.NIST.hc.snps.indels.NIST7035_v_project.NIST.hc.snps.indels.NIST7086.vcf | 237,338 | Common                                             |
| unique_project.NIST.hc.snps.indels.NIST7035.vcf                                        | 99,141  | Unique to project.NIST.hc.snps.indels.NIST7035.vcf |
| unique_project.NIST.hc.snps.indels.NIST7086.vcf                                        | 80,661  | Unique to project.NIST.hc.snps.indels.NIST7086.vcf |

##### Compared with hap.py + RTG tools

- NIST7035 ('baseline') compared to NIST7086 ('truth')

| Type  | Filter | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | QUERY.UNK | FP.gt | FP.al | METRIC.Recall | METRIC.Precision | METRIC.Frac_NA | METRIC.F1_Score | TRUTH.TOTAL.TiTv_ratio | QUERY.TOTAL.TiTv_ratio | TRUTH.TOTAL.het_hom_ratio | QUERY.TOTAL.het_hom_ratio |
|-------|--------|-------------|----------|----------|-------------|----------|-----------|-------|-------|---------------|------------------|----------------|-----------------|------------------------|------------------------|---------------------------|---------------------------|
| INDEL | ALL    | 29465       | 26677    | 2788     | 29855       | 2816     | 110       | 2772  | 44    | 0.905379      | 0.905329         | 0.003684       | 0.905354        |                        |                        | 0.93563464                | 1.004317906               |
| INDEL | PASS   | 29465       | 26677    | 2788     | 29855       | 2816     | 110       | 2772  | 44    | 0.905379      | 0.905329         | 0.003684       | 0.905354        |                        |                        | 0.93563464                | 1.004317906               |
| SNP   | ALL    | 207256      | 194644   | 12612    | 207307      | 12615    | 0         | 12603 | 12    | 0.939148      | 0.939148         | 0              | 0.939148        | 2.101499147            | 2.101166617            | 0.702357185               | 0.720801056               |
| SNP   | PASS   | 207256      | 194644   | 12612    | 207307      | 12615    | 0         | 12603 | 12    | 0.939148      | 0.939148         | 0              | 0.939148        | 2.101499147            | 2.101166617            | 0.702357185               | 0.720801056               |

- NIST7086 ('baseline') compared to NIST7035 ('truth')

| Type  | Filter | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | QUERY.UNK | FP.gt | FP.al | METRIC.Recall | METRIC.Precision | METRIC.Frac_NA | METRIC.F1_Score | TRUTH.TOTAL.TiTv_ratio | QUERY.TOTAL.TiTv_ratio | TRUTH.TOTAL.het_hom_ratio | QUERY.TOTAL.het_hom_ratio |
|-------|--------|-------------|----------|----------|-------------|----------|-----------|-------|-------|---------------|------------------|----------------|-----------------|------------------------|------------------------|---------------------------|---------------------------|
| INDEL | ALL    | 29467       | 26682    | 2785     | 29849       | 2806     | 109       | 2767  | 39    | 0.905487      | 0.905649         | 0.003652       | 0.905568        |                        |                        | 0.961824827               | 0.976642044               |
| INDEL | PASS   | 29467       | 26682    | 2785     | 29849       | 2806     | 109       | 2767  | 39    | 0.905487      | 0.905649         | 0.003652       | 0.905568        |                        |                        | 0.961824827               | 0.976642044               |
| SNP   | ALL    | 207255      | 194643   | 12612    | 207307      | 12616    | 0         | 12604 | 12    | 0.939147      | 0.939143         | 0              | 0.939145        | 2.101641405            | 2.101009497            | 0.720383593               | 0.702762213               |
| SNP   | PASS   | 207255      | 194643   | 12612    | 207307      | 12616    | 0         | 12604 | 12    | 0.939147      | 0.939143         | 0              | 0.939145        | 2.101641405            | 2.101009497            | 0.720383593               | 0.702762213               |

## bench1.0

### Run parameters/settings

- **Aim:** Benchmarking pipelines against the Genome In A Bottle (GIAB) sample [NA12878](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) (both NIST7035 and NIST7086) using the same input data were used to create the truth vcf and using hte default pipeline settings

- **Pipelines:**
  - [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) + minimal [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline) (no annotation) an
  - [parabricks germline pipeline](https://www.nvidia.com/en-us/docs/parabricks/germline/)

- **Inputs human_genomics_pipeline:**
  - Reference genome: ucsc.hg19.fasta
  - dbSNP database: dbsnp_138.hg19.vcf

- **Inputs vcf_annotation_pipeline:**
  - Reference genome: ucsc.hg19.fasta
  - dbSNP database: dbsnp_138.hg19.vcf
  - Hapmap: hapmap_3.3.hg19.sites.vcf.gz
  - Mills: Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
  
- **Inputs for parabricks germline pipeline:**
  - Reference genome: ucsc.hg19.fasta
  - dbSNP database: dbsnp_138.hg19.vcf

- **Other settings:**
  - No padding
  - No intervals
  - 2D model with pre-trained architecture (for rule gatk CNNScoreVariants)
  - SNP tranches: 99.95 (for rule gatk4_FilterVariantTranches)
  - Indel tranches: 99.40 (for rule gatk4_FilterVariantTranches)

*Note. dbsnp_137.hg19.vcf doesn't appear to be available anymore, therefore I used the closest version available*

(see run settings/output for [human_genomics_pipeline - bench1.0](https://github.com/ESR-NZ/human_genomics_pipeline/tree/bench1.0), [vcf_annotation_pipeline - bench1.0](https://github.com/ESR-NZ/vcf_annotation_pipeline/tree/bench1.0) and [parabricks - bench1.0]() for more detail)

### Results

Results dir: /store/lkemp/exome_project/benchmarking/NA12878_exome/bench1.0/ (wintermute)

#### human_genomics_pipeline + minimal vcf_annotation_pipeline

##### Compared with bedops intersect

- NIST7035

| file                                                                     | count   | type                                               |
|--------------------------------------------------------------------------|---------|----------------------------------------------------|
| common_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035.vcf | 188,127 | Common                                             |
| unique_NIST7035_NIST_filtered.vcf                                        | 280,797 | Unique to h_g_p + v_a_p output                     |
| unique_project.NIST.hc.snps.indels.NIST7035.vcf                          | 148,216 | Unique to known vcf                                |

- NIST7086

| file                                                                     | count   | type                                               |
|--------------------------------------------------------------------------|---------|----------------------------------------------------|
| common_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086.vcf | 190,599 | Common                                             |
| unique_NIST7086_NIST_filtered.vcf                                        | 311,418 | Unique to h_g_p + v_a_p output                     |
| unique_project.NIST.hc.snps.indels.NIST7086.vcf                          | 127,323 | Unique to known vcf                                |

##### Compared with hap.py + RTG tools

- NIST7035

| Type  | Filter | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | QUERY.UNK | FP.gt | FP.al | METRIC.Recall | METRIC.Precision | METRIC.Frac_NA | METRIC.F1_Score | TRUTH.TOTAL.TiTv_ratio | QUERY.TOTAL.TiTv_ratio | TRUTH.TOTAL.het_hom_ratio | QUERY.TOTAL.het_hom_ratio |
|-------|--------|-------------|----------|----------|-------------|----------|-----------|-------|-------|---------------|------------------|----------------|-----------------|------------------------|------------------------|---------------------------|---------------------------|
| INDEL | ALL    | 29473       | 20889    | 8584     | 57705       | 1684     | 34942     | 687   | 848   | 0.70875       | 0.92602          | 0.605528       | 0.802947        |                        |                        | 0.935843514               | 0.530099933               |
| INDEL | PASS   | 29473       | 20226    | 9247     | 56456       | 1597     | 34445     | 679   | 780   | 0.686255      | 0.927445         | 0.610121       | 0.788826        |                        |                        | 0.935843514               | 0.498892595               |
| SNP   | ALL    | 207268      | 163391   | 43877    | 411540      | 1985     | 246129    | 1112  | 372   | 0.788308      | 0.988            | 0.598068       | 0.876929        | 2.101415235            | 1.423485638            | 0.702447561               | 0.362133473               |
| SNP   | PASS   | 207268      | 161764   | 45504    | 408910      | 1908     | 245202    | 1082  | 347   | 0.780458      | 0.988345         | 0.599648       | 0.872185        | 2.101415235            | 1.422810364            | 0.702447561               | 0.35380422                |

- NIST7086

| Type  | Filter | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | QUERY.UNK | FP.gt | FP.al | METRIC.Recall | METRIC.Precision | METRIC.Frac_NA | METRIC.F1_Score | TRUTH.TOTAL.TiTv_ratio | QUERY.TOTAL.TiTv_ratio | TRUTH.TOTAL.het_hom_ratio | QUERY.TOTAL.het_hom_ratio |
|-------|--------|-------------|----------|----------|-------------|----------|-----------|-------|-------|---------------|------------------|----------------|-----------------|------------------------|------------------------|---------------------------|---------------------------|
| INDEL | ALL    | 29475       | 21439    | 8036     | 62207       | 1750     | 38823     | 744   | 883   | 0.727362      | 0.925163         | 0.624094       | 0.814424        |                        |                        | 0.961967923               | 0.537195213               |
| INDEL | PASS   | 29475       | 20764    | 8711     | 60928       | 1658     | 38312     | 735   | 812   | 0.704461      | 0.926689         | 0.628808       | 0.800437        |                        |                        | 0.961967923               | 0.507733625               |
| SNP   | ALL    | 207264      | 165322   | 41942    | 440207      | 1947     | 272906    | 1187  | 339   | 0.79764       | 0.988362         | 0.619949       | 0.882818        | 2.101558993            | 1.419227578            | 0.720450017               | 0.3631554                 |
| SNP   | PASS   | 207264      | 163448   | 43816    | 437144      | 1872     | 271792    | 1151  | 327   | 0.788598      | 0.988679         | 0.621745       | 0.877376        | 2.101558993            | 1.418276144            | 0.720450017               | 0.354108053               |

#### parabricks germline pipeline

##### Compared with bedops intersect

- NIST7035

| file                                                                     | count   | type                                               |
|--------------------------------------------------------------------------|---------|----------------------------------------------------|
| common_NIST7035_NIST_v_project.NIST.hc.snps.indels.NIST7035.vcf          | 188,161 | Common                                             |
| unique_NIST7035_NIST.vcf                                                 | 280,446 | Unique to parabricks output                        |
| unique_project.NIST.hc.snps.indels.NIST7035.vcf                          | 148,197 | Unique to known vcf                                |

- NIST7086

| file                                                                     | count   | type                                               |
|--------------------------------------------------------------------------|---------|----------------------------------------------------|
| common_NIST7086_NIST_v_project.NIST.hc.snps.indels.NIST7086.vcf          | 190,707 | Common                                             |
| unique_NIST7086_NIST.vcf                                                 | 311,598 | Unique to parabricks output                        |
| unique_project.NIST.hc.snps.indels.NIST7086.vcf                          | 127,260 | Unique to known vcf                                |

##### Compared with hap.py + RTG tools

- NIST7035

| Type  | Filter | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | QUERY.UNK | FP.gt | FP.al | METRIC.Recall | METRIC.Precision | METRIC.Frac_NA | METRIC.F1_Score | TRUTH.TOTAL.TiTv_ratio | QUERY.TOTAL.TiTv_ratio | TRUTH.TOTAL.het_hom_ratio | QUERY.TOTAL.het_hom_ratio |
|-------|--------|-------------|----------|----------|-------------|----------|-----------|-------|-------|---------------|------------------|----------------|-----------------|------------------------|------------------------|---------------------------|---------------------------|
| INDEL | ALL    | 29473       | 20616    | 8857     | 54367       | 1153     | 32444     | 536   | 491   | 0.699488      | 0.947407         | 0.596759       | 0.804787        |                        |                        | 0.935843514               | 0.450247227               |
| INDEL | PASS   | 29473       | 20616    | 8857     | 54367       | 1153     | 32444     | 536   | 491   | 0.699488      | 0.947407         | 0.596759       | 0.804787        |                        |                        | 0.935843514               | 0.450247227               |
| SNP   | ALL    | 207268      | 163947   | 43321    | 414135      | 2033     | 248128    | 1158  | 375   | 0.79099       | 0.987754         | 0.599148       | 0.878489        | 2.101415235            | 1.42376204             | 0.702447561               | 0.360109056               |
| SNP   | PASS   | 207268      | 163947   | 43321    | 414135      | 2033     | 248128    | 1158  | 375   | 0.79099       | 0.987754         | 0.599148       | 0.878489        | 2.101415235            | 1.42376204             | 0.702447561               | 0.360109056               |

- NIST7086

| Type  | Filter | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | QUERY.UNK | FP.gt | FP.al | METRIC.Recall | METRIC.Precision | METRIC.Frac_NA | METRIC.F1_Score | TRUTH.TOTAL.TiTv_ratio | QUERY.TOTAL.TiTv_ratio | TRUTH.TOTAL.het_hom_ratio | QUERY.TOTAL.het_hom_ratio |
|-------|--------|-------------|----------|----------|-------------|----------|-----------|-------|-------|---------------|------------------|----------------|-----------------|------------------------|------------------------|---------------------------|---------------------------|
| INDEL | ALL    | 29473       | 21269    | 8204     | 58506       | 1111     | 35951     | 511   | 495   | 0.721644      | 0.950743         | 0.614484       | 0.820501        |                        |                        | 0.961832578               | 0.458301072               |
| INDEL | PASS   | 29473       | 21269    | 8204     | 58506       | 1111     | 35951     | 511   | 495   | 0.721644      | 0.950743         | 0.614484       | 0.820501        |                        |                        | 0.961832578               | 0.458301072               |
| SNP   | ALL    | 207253      | 165941   | 41312    | 443676      | 1940     | 275768    | 1194  | 333   | 0.800669      | 0.988446         | 0.621553       | 0.884703        | 2.101533628            | 1.419960733            | 0.720358685               | 0.361391578               |
| SNP   | PASS   | 207253      | 165941   | 41312    | 443676      | 1940     | 275768    | 1194  | 333   | 0.800669      | 0.988446         | 0.621553       | 0.884703        | 2.101533628            | 1.419960733            | 0.720358685               | 0.361391578               |

## bench1.1

### Run parameters/settings

- **Aim:** See if the quality of variants (tp, fp, fn) are different when variants are filtered using the 1D model with pre-trained architecture (for rule gatk CNNScoreVariants) compared to the 2D model (see description of these models [here](https://gatk.broadinstitute.org/hc/en-us/articles/360037226672-CNNScoreVariants))

- **Pipelines:**
  - [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) + minimal [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline) (no annotation)

- **Inputs human_genomics_pipeline:**
  - Reference genome: ucsc.hg19.fasta
  - dbSNP database: dbsnp_138.hg19.vcf

- **Inputs vcf_annotation_pipeline:**
  - Reference genome: ucsc.hg19.fasta
  - dbSNP database: dbsnp_138.hg19.vcf
  - Hapmap: hapmap_3.3.hg19.sites.vcf.gz
  - Mills: Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
  
- **Inputs for parabricks germline pipeline:**
  - Reference genome: ucsc.hg19.fasta
  - dbSNP database: dbsnp_138.hg19.vcf

- **Other settings:**
  - No padding
  - No intervals
  - :star: 1D model with pre-trained architecture (for rule gatk CNNScoreVariants) :star:
  - SNP tranches: 99.95 (for rule gatk4_FilterVariantTranches)
  - Indel tranches: 99.40 (for rule gatk4_FilterVariantTranches)

*Note. because the settings in the first pipeline (human_genomics_pipeline) was not changed in this bench1.0 (just changes to settings in vcf_annotation_pipeline), I used the human_genomics_pipeline output from bench1.0*

(see run settings/output for [human_genomics_pipeline - bench1.1](https://github.com/ESR-NZ/human_genomics_pipeline/tree/bench1.1) and [vcf_annotation_pipeline - bench1.1](https://github.com/ESR-NZ/vcf_annotation_pipeline/tree/bench1.1) for more detail)

### Results

Results dir: /store/lkemp/exome_project/benchmarking/NA12878_exome/bench1.1/ (wintermute)

#### human_genomics_pipeline + minimal vcf_annotation_pipeline

##### Compared with bedops intersect

- NIST7035

| file                                                                     | count   | type                                               |
|--------------------------------------------------------------------------|---------|----------------------------------------------------|
| common_NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7035.vcf | 188,127 | Common                                             |
| unique_NIST7035_NIST_filtered.vcf                                        | 280,797 | Unique to NIST7035_NIST_filtered                   |
| unique_project.NIST.hc.snps.indels.NIST7035.vcf                          | 148,216 | Unique to project.NIST.hc.snps.indels.NIST7035.vcf |

- NIST7086

| file                                                                     | count   | type                                               |
|--------------------------------------------------------------------------|---------|----------------------------------------------------|
| common_NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels.NIST7086.vcf | 190,599 | Common                                             |
| unique_NIST7086_NIST_filtered.vcf                                        | 311,418 | Unique to NIST7086_NIST_filtered                   |
| unique_project.NIST.hc.snps.indels.NIST7086.vcf                          | 127,323 | Unique to project.NIST.hc.snps.indels.NIST7086.vcf |

##### Compared with hap.py + RTG tools

- NIST7035

| Type  | Filter | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | QUERY.UNK | FP.gt | FP.al | METRIC.Recall | METRIC.Precision | METRIC.Frac_NA | METRIC.F1_Score | TRUTH.TOTAL.TiTv_ratio | QUERY.TOTAL.TiTv_ratio | TRUTH.TOTAL.het_hom_ratio | QUERY.TOTAL.het_hom_ratio |
|-------|--------|-------------|----------|----------|-------------|----------|-----------|-------|-------|---------------|------------------|----------------|-----------------|------------------------|------------------------|---------------------------|---------------------------|
| INDEL | ALL    | 29473       | 20889    | 8584     | 57705       | 1684     | 34942     | 687   | 848   | 0.70875       | 0.92602          | 0.605528       | 0.802947        |                        |                        | 0.935843514               | 0.530099933               |
| INDEL | PASS   | 29473       | 20420    | 9053     | 56725       | 1617     | 34499     | 682   | 793   | 0.692838      | 0.927247         | 0.60818        | 0.793084        |                        |                        | 0.935843514               | 0.504986667               |
| SNP   | ALL    | 207268      | 163391   | 43877    | 411540      | 1985     | 246129    | 1112  | 372   | 0.788308      | 0.988            | 0.598068       | 0.876929        | 2.101415235            | 1.423485638            | 0.702447561               | 0.362133473               |
| SNP   | PASS   | 207268      | 161444   | 45824    | 408308      | 1905     | 244923    | 1084  | 347   | 0.778914      | 0.98834          | 0.599849       | 0.871218        | 2.101415235            | 1.421803843            | 0.702447561               | 0.351807141               |

- NIST7086

| Type  | Filter | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | QUERY.UNK | FP.gt | FP.al | METRIC.Recall | METRIC.Precision | METRIC.Frac_NA | METRIC.F1_Score | TRUTH.TOTAL.TiTv_ratio | QUERY.TOTAL.TiTv_ratio | TRUTH.TOTAL.het_hom_ratio | QUERY.TOTAL.het_hom_ratio |
|-------|--------|-------------|----------|----------|-------------|----------|-----------|-------|-------|---------------|------------------|----------------|-----------------|------------------------|------------------------|---------------------------|---------------------------|
| INDEL | ALL    | 29475       | 21439    | 8036     | 62207       | 1750     | 38823     | 744   | 883   | 0.727362      | 0.925163         | 0.624094       | 0.814424        |                        |                        | 0.961967923               | 0.537195213               |
| INDEL | PASS   | 29475       | 20904    | 8571     | 61136       | 1680     | 38358     | 735   | 832   | 0.709211      | 0.926245         | 0.627421       | 0.803327        |                        |                        | 0.961967923               | 0.511903579               |
| SNP   | ALL    | 207264      | 165322   | 41942    | 440207      | 1947     | 272906    | 1187  | 339   | 0.79764       | 0.988362         | 0.619949       | 0.882818        | 2.101558993            | 1.419227578            | 0.720450017               | 0.3631554                 |
| SNP   | PASS   | 207264      | 163334   | 43930    | 436902      | 1876     | 271660    | 1158  | 327   | 0.788048      | 0.988647         | 0.621787       | 0.877023        | 2.101558993            | 1.417673127            | 0.720450017               | 0.353316271               |

## bench1.2

### Run parameters/settings

- **Aim:** Improve the accuracy of variant filtering in the gatk4_FilterVariantTranches rule (see description of this step [here](https://gatk.broadinstitute.org/hc/en-us/articles/360037227632-FilterVariantTranches)). I'll use a number of difference tranches to get a gauge of what kind of sensitivity will be needed for the filtering of snps and indels.

- **Pipelines:**
  - [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) + minimal [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline) (no annotation)
  - [parabricks germline pipeline](https://www.nvidia.com/en-us/docs/parabricks/germline/)

- **Inputs human_genomics_pipeline:**
  - Reference genome: ucsc.hg19.fasta
  - dbSNP database: dbsnp_138.hg19.vcf

- **Inputs vcf_annotation_pipeline:**
  - Reference genome: ucsc.hg19.fasta
  - dbSNP database: dbsnp_138.hg19.vcf
  - Hapmap: hapmap_3.3.hg19.sites.vcf.gz
  - Mills: Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
  
- **Inputs for parabricks germline pipeline:**
  - Reference genome: ucsc.hg19.fasta
  - dbSNP database: dbsnp_138.hg19.vcf

- **Other settings:**
  - No padding
  - No intervals
  - 2D model with pre-trained architecture (for rule gatk CNNScoreVariants)
  - :star: SNP tranches: 80.00 81.00 82.00 83.00 84.00 85.00 86.00 87.00 88.00 89.00 90.00 91.00 92.00 93.00 94.00 95.00 96.00 97.00 98.00 99.00 (for rule gatk4_FilterVariantTranches) :star:
  - :star: Indel tranches: 80.00 81.00 82.00 83.00 84.00 85.00 86.00 87.00 88.00 89.00 90.00 91.00 92.00 93.00 94.00 95.00 96.00 97.00 98.00 99.00 (for rule gatk4_FilterVariantTranches) :star:

*Note. because the settings in the first pipeline (human_genomics_pipeline) was not changed in this bench1.0 (just changes to settings in vcf_annotation_pipeline), I used the human_genomics_pipeline output from bench1.0*

(see run settings/output for [human_genomics_pipeline - bench1.2](https://github.com/ESR-NZ/human_genomics_pipeline/tree/bench1.2), [vcf_annotation_pipeline - bench1.2](https://github.com/ESR-NZ/vcf_annotation_pipeline/tree/bench1.2) and [parabricks - bench1.2]() for more detail)

### Results

Results dir: /store/lkemp/exome_project/benchmarking/NA12878_exome/bench1.2/ (Wintermute) and /home/lkemp/benchmarking_resources/ (Methead)

Variants from tranche 82.00 or less were extracted for the following comparisons. 'PASS' indicates a variant in a tranche below 80.00.

#### human_genomics_pipeline + minimal vcf_annotation_pipeline

##### Compared with bedops intersect

- NIST7035

| file                                                                                     | count   | type                                               |
|------------------------------------------------------------------------------------------|---------|----------------------------------------------------|
| common_NIST7035_NIST_filtered_less_than_82.00_v_project.NIST.hc.snps.indels.NIST7035.vcf | 139,328 | Common                                             |
| unique_NIST7035_NIST_filtered_less_than_82.00.vcf                                        | 172,013 | Unique to h_g_p + v_a_p output                     |
| unique_project.NIST.hc.snps.indels.NIST7035.vcf                                          | 196,693 | Unique to known vcf                                |

- NIST7086

| file                                                                                     | count   | type                                               |
|------------------------------------------------------------------------------------------|---------|----------------------------------------------------|
| common_NIST7086_NIST_filtered_less_than_82.00_v_project.NIST.hc.snps.indels.NIST7086.vcf | 140,929 | Common                                             |
| unique_NIST7086_NIST_filtered_less_than_82.00.vcf                                        | 192,300 | Unique to h_g_p + v_a_p output                     |
| unique_project.NIST.hc.snps.indels.NIST7086.vcf                                          | 176,602 | Unique to known vcf                                |

##### Compared with hap.py + RTG tools

- NIST7035

| Type  | Filter | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | QUERY.UNK | FP.gt | FP.al | METRIC.Recall | METRIC.Precision | METRIC.Frac_NA | METRIC.F1_Score | TRUTH.TOTAL.TiTv_ratio | QUERY.TOTAL.TiTv_ratio | TRUTH.TOTAL.het_hom_ratio | QUERY.TOTAL.het_hom_ratio |
|-------|--------|-------------|----------|----------|-------------|----------|-----------|-------|-------|---------------|------------------|----------------|-----------------|------------------------|------------------------|---------------------------|---------------------------|
| INDEL | ALL    | 29470       | 14545    | 14925    | 40762       | 783      | 25280     | 421   | 323   | 0.493553      | 0.949425         | 0.620185       | 0.649478        |                        |                        | 0.935772466               | 0.344969335               |
| INDEL | PASS   | 29470       | 14138    | 15332    | 39434       | 750      | 24395     | 411   | 307   | 0.479742      | 0.95013          | 0.618629       | 0.637564        |                        |                        | 0.935772466               | 0.339235756               |
| SNP   | ALL    | 207268      | 123729   | 83539    | 271067      | 422      | 146881    | 315   | 51    | 0.596952      | 0.996602         | 0.541862       | 0.746662        | 2.101415235            | 1.717486092            | 0.702447561               | 0.252478245               |
| SNP   | PASS   | 207268      | 121830   | 85438    | 261115      | 405      | 138845    | 306   | 50    | 0.58779       | 0.996688         | 0.531739       | 0.739478        | 2.101415235            | 1.763760279            | 0.702447561               | 0.258217343               |

- NIST7086

| Type  | Filter | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | QUERY.UNK | FP.gt | FP.al | METRIC.Recall | METRIC.Precision | METRIC.Frac_NA | METRIC.F1_Score | TRUTH.TOTAL.TiTv_ratio | QUERY.TOTAL.TiTv_ratio | TRUTH.TOTAL.het_hom_ratio | QUERY.TOTAL.het_hom_ratio |
|-------|--------|-------------|----------|----------|-------------|----------|-----------|-------|-------|---------------|------------------|----------------|-----------------|------------------------|------------------------|---------------------------|---------------------------|
| INDEL | ALL    | 29471       | 14531    | 14940    | 43800       | 871      | 28247     | 457   | 373   | 0.493061      | 0.943998         | 0.644909       | 0.647779        |                        |                        | 0.962095573               | 0.350436614               |
| INDEL | PASS   | 29471       | 14152    | 15319    | 42375       | 840      | 27235     | 443   | 360   | 0.480201      | 0.944518         | 0.642714       | 0.636699        |                        |                        | 0.962095573               | 0.345535714               |
| SNP   | ALL    | 207267      | 125305   | 81962    | 289995      | 403      | 164254    | 311   | 39    | 0.604558      | 0.996795         | 0.566403       | 0.752639        | 2.101557474            | 1.727640365            | 0.720474925               | 0.248334919               |
| SNP   | PASS   | 207267      | 123441   | 83826    | 279425      | 386      | 155566    | 302   | 33    | 0.595565      | 0.996884         | 0.556736       | 0.745656        | 2.101557474            | 1.776386902            | 0.720474925               | 0.253609857               |

#### parabricks germline pipeline?

##### Compared with bedops intersect

- NIST7035

- NIST7086

##### Compared with hap.py + RTG tools

- NIST7035

- NIST7086

## bench1.3

### Run parameters/settings

- **Aim:** Improve the accuracy of variant filtering in the gatk4_FilterVariantTranches rule (see description of this step [here](https://gatk.broadinstitute.org/hc/en-us/articles/360037227632-FilterVariantTranches)). There are still a number of false positive being called even when the tranches (for both snps and indels) are reduced to from the defaults of around 99.00 to 82.00. I will try another (lower) set of tranches to get a gauge of what kind of sensitivity will be needed for the filtering of snps and indels.

- **Pipelines:**
  - [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) + minimal [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline) (no annotation)
  - [parabricks germline pipeline](https://www.nvidia.com/en-us/docs/parabricks/germline/)

- **Inputs human_genomics_pipeline:**
  - Reference genome: ucsc.hg19.fasta
  - dbSNP database: dbsnp_138.hg19.vcf

- **Inputs vcf_annotation_pipeline:**
  - Reference genome: ucsc.hg19.fasta
  - dbSNP database: dbsnp_138.hg19.vcf
  - Hapmap: hapmap_3.3.hg19.sites.vcf.gz
  - Mills: Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
  
- **Inputs for parabricks germline pipeline:**
  - Reference genome: ucsc.hg19.fasta
  - dbSNP database: dbsnp_138.hg19.vcf

- **Other settings:**
  - No padding
  - No intervals
  - 2D model with pre-trained architecture (for rule gatk CNNScoreVariants)
  - :star: SNP tranches: 70.00 71.00 72.00 73.00 74.00 75.00 76.00 77.00 78.00 79.00 80.00 81.00 82.00 (for rule gatk4_FilterVariantTranches) :star:
  - :star: Indel tranches: 70.00 71.00 72.00 73.00 74.00 75.00 76.00 77.00 78.00 79.00 80.00 81.00 82.00 (for rule gatk4_FilterVariantTranches) :star:

*Note. because the settings in the first pipeline (human_genomics_pipeline) was not changed in this bench1.0 (just changes to settings in vcf_annotation_pipeline), I used the human_genomics_pipeline output from bench1.0*

(see run settings/output for [human_genomics_pipeline - bench1.3](https://github.com/ESR-NZ/human_genomics_pipeline/tree/bench1.3), [vcf_annotation_pipeline - bench1.3](https://github.com/ESR-NZ/vcf_annotation_pipeline/tree/bench1.3) and [parabricks - bench1.3]() for more detail)

### Results

Results dir: /store/lkemp/exome_project/benchmarking/NA12878_exome/bench1.3/ (Wintermute) and /home/lkemp/benchmarking_resources/ (Methead)

Variants from tranche 82.00 or less were extracted for the following comparisons. 'PASS' indicates a variant in a tranche below 70.00.

#### human_genomics_pipeline + minimal vcf_annotation_pipeline

##### Compared with bedops intersect

- NIST7035

| file                                                                                     | count   | type                                               |
|------------------------------------------------------------------------------------------|---------|----------------------------------------------------|
| common_NIST7035_NIST_filtered_less_than_82.00_v_project.NIST.hc.snps.indels.NIST7035.vcf | 139,328 | Common                                             |
| unique_NIST7035_NIST_filtered_less_than_82.00.vcf                                        | 172,013 | Unique to h_g_p + v_a_p output                     |
| unique_project.NIST.hc.snps.indels.NIST7035.vcf                                          | 196,693 | Unique to known vcf                                |

- NIST7086

| file                                                                                     | count   | type                                               |
|------------------------------------------------------------------------------------------|---------|----------------------------------------------------|
| common_NIST7086_NIST_filtered_less_than_82.00_v_project.NIST.hc.snps.indels.NIST7086.vcf | 140,929 | Common                                             |
| unique_NIST7086_NIST_filtered_less_than_82.00.vcf                                        | 192,300 | Unique to h_g_p + v_a_p output                     |
| unique_project.NIST.hc.snps.indels.NIST7086.vcf                                          | 176,602 | Unique to known vcf                                |

##### Compared with hap.py + RTG tools

- NIST7035

| Type  | Filter | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | QUERY.UNK | FP.gt | FP.al | METRIC.Recall | METRIC.Precision | METRIC.Frac_NA | METRIC.F1_Score | TRUTH.TOTAL.TiTv_ratio | QUERY.TOTAL.TiTv_ratio | TRUTH.TOTAL.het_hom_ratio | QUERY.TOTAL.het_hom_ratio |
|-------|--------|-------------|----------|----------|-------------|----------|-----------|-------|-------|---------------|------------------|----------------|-----------------|------------------------|------------------------|---------------------------|---------------------------|
| INDEL | ALL    | 29470       | 14545    | 14925    | 40762       | 783      | 25280     | 421   | 323   | 0.493553      | 0.949425         | 0.620185       | 0.649478        |                        |                        | 0.935772466               | 0.344969335               |
| INDEL | PASS   | 29470       | 12356    | 17114    | 33299       | 628      | 20182     | 348   | 256   | 0.419274      | 0.952123         | 0.606084       | 0.582181        |                        |                        | 0.935772466               | 0.314825097               |
| SNP   | ALL    | 207268      | 123729   | 83539    | 271067      | 422      | 146881    | 315   | 51    | 0.596952      | 0.996602         | 0.541862       | 0.746662        | 2.101415235            | 1.717486092            | 0.702447561               | 0.252478245               |
| SNP   | PASS   | 207268      | 113396   | 93872    | 214503      | 345      | 100727    | 270   | 43    | 0.547098      | 0.996968         | 0.469583       | 0.706498        | 2.101415235            | 2.07447586             | 0.702447561               | 0.309019661               |

- NIST7086

| Type  | Filter | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | QUERY.UNK | FP.gt | FP.al | METRIC.Recall | METRIC.Precision | METRIC.Frac_NA | METRIC.F1_Score | TRUTH.TOTAL.TiTv_ratio | QUERY.TOTAL.TiTv_ratio | TRUTH.TOTAL.het_hom_ratio | QUERY.TOTAL.het_hom_ratio |
|-------|--------|-------------|----------|----------|-------------|----------|-----------|-------|-------|---------------|------------------|----------------|-----------------|------------------------|------------------------|---------------------------|---------------------------|
| INDEL | ALL    | 29471       | 14531    | 14940    | 43800       | 871      | 28247     | 457   | 373   | 0.493061      | 0.943998         | 0.644909       | 0.647779        |                        |                        | 0.962095573               | 0.350436614               |
| INDEL | PASS   | 29471       | 12517    | 16954    | 35865       | 706      | 22509     | 386   | 291   | 0.424723      | 0.94714          | 0.627604       | 0.586461        |                        |                        | 0.962095573               | 0.325171582               |
| SNP   | ALL    | 207267      | 125305   | 81962    | 289995      | 403      | 164254    | 311   | 39    | 0.604558      | 0.996795         | 0.566403       | 0.752639        | 2.101557474            | 1.727640365            | 0.720474925               | 0.248334919               |
| SNP   | PASS   | 207267      | 115544   | 91723    | 230523      | 347      | 114600    | 280   | 28    | 0.557465      | 0.997007         | 0.49713        | 0.715093        | 2.101557474            | 2.08224489             | 0.720474925               | 0.29987085                |