# Benchmarking results

## Table of contents

- [Benchmarking results](#benchmarking-results)
  - [Table of contents](#table-of-contents)
  - [Known vcf](#known-vcf)
  - [bench1.0](#bench10)
    - [Run parameters/settings](#run-parameterssettings)
    - [Results](#results)
      - [NIST7035](#nist7035)
      - [NIST7086](#nist7086)
  - [bench1.1](#bench11)
    - [Run parameters/settings](#run-parameterssettings-1)
    - [Results](#results-1)
      - [NIST7035](#nist7035-1)
      - [NIST7086](#nist7086-1)

## Known vcf

Parameters used to create the Genome In A Bottle (GIAB) sample [NA12878](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) (both NIST7035 and NIST7086) known vcf file:

- **Inputs:**
  - Reference genome: ucsc.hg19.fasta
  - dbSNP database: dbsnp_137.hg19.vcf

- **Description:**
  - Number of variants: 416,689

## bench1.0

### Run parameters/settings

- **Aim:** Benchmarking against against the Genome In A Bottle (GIAB) sample [NA12878](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) (both NIST7035 and NIST7086)
- **Pipeline:** [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) and minimal [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline) (no annotation)
- **Inputs human_genomics_pipeline:**
  - Reference genome: ucsc.hg19.fasta
  - dbSNP database: dbsnp_138.hg19.vcf

- **Inputs vcf_annotation_pipeline:**
  - Reference genome: ucsc.hg19.fasta
  - dbSNP database: dbsnp_138.hg19.vcf
  - Hapmap: hapmap_3.3.hg19.sites.vcf.gz
  - Mills: Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz

- **Other settings:**
  - No padding
  - No intervals
  - 2D model with pre-trained architecture (for rule gatk CNNScoreVariants)

*Note. dbsnp_137.hg19.vcf doesn't appear to be available anymore, therefore I used the closest version available*

(see run settings/output for [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline/tree/bench1.0) and [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline/tree/bench1.0) for more detail)

### Results

#### NIST7035

- Compared with bcftool isec

| file                                                          | count   | count (chromosome overlap adjusted) | type      |
|---------------------------------------------------------------|---------|-------------------------------------|-----------|
| NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels/0000.vcf | 215,888 | 215,875                             | False-pos |
| NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels/0001.vcf | 164,030 | 164,027                             | False-neg |
| NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels/0002.vcf | 252,659 | 252,659                             | True-pos  |
| NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels/0003.vcf | 252,659 | 252,659                             | True-pos  |

- Compared with rtg vcfeval

|Threshold | True-pos-baseline   | True-pos-call  | False-pos | False-neg | Precision | Sensitivity | F-measure |
|----------|---------------------|----------------|-----------|-----------|-----------|-------------|-----------|
| 1.000    | 52,666              | 52,667         | 2957      | 4030      | 0.9468    | 0.9289      | 0.9378    |
| None     | 52,668              | 52,669         | 2961      | 4028      | 0.9468    | 0.9290      | 0.9378    |

- Compared with hap.py

| Type  | Filter | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | QUERY.UNK | FP.gt | FP.al | METRIC.Recall | METRIC.Precision | METRIC.Frac_NA | METRIC.F1_Score | TRUTH.TOTAL.TiTv_ratio | QUERY.TOTAL.TiTv_ratio | TRUTH.TOTAL.het_hom_ratio | QUERY.TOTAL.het_hom_ratio |
|-------|--------|-------------|----------|----------|-------------|----------|-----------|-------|-------|---------------|------------------|----------------|-----------------|------------------------|------------------------|---------------------------|---------------------------|
| INDEL | ALL    | 33933       | 32145    | 1788     | 57706       | 1687     | 23671     | 390   | 1102  | 0.947308      | 0.950433         | 0.4102         | 0.948868        |                        |                        | 0.575677711               | 0.530085807               |
| INDEL | PASS   | 33933       | 31461    | 2472     | 56457       | 1597     | 23195     | 384   | 1044  | 0.927151      | 0.951987         | 0.410844       | 0.939405        |                        |                        | 0.575677711               | 0.498879283               |
| SNP   | ALL    | 219952      | 216836   | 3116     | 411540      | 190697   | 3974      | 1365  | 2846  | 0.985833      | 0.532108         | 0.009656       | 0.691159        | 2.032894465            | 1.423485638            | 0.577210523               | 0.362133473               |
| SNP   | PASS   | 219952      | 215188   | 4764     | 408910      | 189819   | 3870      | 1336  | 2760  | 0.978341      | 0.531357         | 0.009464       | 0.688679        | 2.032894465            | 1.422810364            | 0.577210523               | 0.35380422                |

#### NIST7086

- Compared with bcftools isec

| file                                                          | count   | count (chromosome overlap adjusted) | type      |
|---------------------------------------------------------------|---------|-------------------------------------|-----------|
| NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels/0000.vcf | 230,065 | 230,053                             | False-pos |
| NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels/0001.vcf | 145,127 | 145,121                             | False-neg |
| NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels/0002.vcf | 271,562 | 271,562                             | True-pos  |
| NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels/0003.vcf | 271,562 | 271,562                             | True-pos  |

- Compared with rtg vcfeval

|Threshold | True-pos-baseline  | True-pos-call  | False-pos | False-neg | Precision | Sensitivity | F-measure |
|----------|--------------------|----------------|-----------|-----------|-----------|-------------|-----------|
| 4.000    | 53,077             | 53,078         | 2872      | 3853      | 0.9487    | 0.9323      | 0.9404    |
| None     | 53,102             | 53,103         | 2912      | 3828      | 0.9480    | 0.9328      | 0.9403    |

## bench1.1

### Run parameters/settings

- **Aim:** Benchmarking against against the Genome In A Bottle (GIAB) sample [NA12878](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) (both NIST7035 and NIST7086)
- **Pipeline:** [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline) and minimal [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline) (no annotation)
- **Inputs human_genomics_pipeline:**
  - Reference genome: ucsc.hg19.fasta
  - dbSNP database: dbsnp_138.hg19.vcf

- **Inputs vcf_annotation_pipeline:**
  - Reference genome: ucsc.hg19.fasta
  - dbSNP database: dbsnp_138.hg19.vcf
  - Hapmap: hapmap_3.3.hg19.sites.vcf.gz
  - Mills: Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz

- **Other settings:**
  - No padding
  - No intervals
  - 1D model with pre-trained architecture (for rule gatk CNNScoreVariants)

*Note. dbsnp_137.hg19.vcf doesn't appear to be available anymore, therefore I used the closest version available*

(see run settings/output for [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline/tree/bench1.0) and [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline/tree/bench1.0) for more detail)

### Results

#### NIST7035

- Compared with bcftool isec

| file                                                          | count   | count (chromosome overlap adjusted) | type      |
|---------------------------------------------------------------|---------|-------------------------------------|-----------|
| NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels/0000.vcf | 215,888 | 215,875                             | False-pos |
| NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels/0001.vcf | 164,030 | 164,027                             | False-neg |
| NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels/0002.vcf | 252,659 | 252,659                             | True-pos  |
| NIST7035_NIST_filtered_v_project.NIST.hc.snps.indels/0003.vcf | 252,659 | 252,659                             | True-pos  |

- Compared with rtg vcfeval

|Threshold | True-pos-baseline  | True-pos-call  | False-pos | False-neg | Precision | Sensitivity | F-measure |
|----------|--------------------|----------------|-----------|-----------|-----------|-------------|-----------|
| 1.000    | 52,717             | 52,718         | 2944      | 3979      | 0.9471    | 0.9298      | 0.9384    |
| None     | 52,719             | 52,720         | 2948      | 3977      | 0.9470    | 0.9299      | 0.9384    |

#### NIST7086

- Compared with bcftools isec

| file                                                          | count   | count (chromosome overlap adjusted) | type      |
|---------------------------------------------------------------|---------|-------------------------------------|-----------|
| NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels/0000.vcf | 230,065 | 230,053                             | False-pos |
| NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels/0001.vcf | 145,127 | 145,121                             | False-neg |
| NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels/0002.vcf | 271,562 | 271,562                             | True-pos  |
| NIST7086_NIST_filtered_v_project.NIST.hc.snps.indels/0003.vcf | 271,562 | 271,562                             | True-pos  |

- Compared with rtg vcfeval

|Threshold | True-pos-baseline  | True-pos-call  | False-pos | False-neg | Precision | Sensitivity | F-measure |
|----------|--------------------|----------------|-----------|-----------|-----------|-------------|-----------|
| 4.000    | 53,197             | 53,198         | 2941      | 3733      | 0.9476    | 0.9344      | 0.9410    |
| None     | 53,222             | 53,223         | 2981      | 3708      | 0.9470    | 0.9349      | 0.9409    |
