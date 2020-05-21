# Benchmarking results

## Table of contents

- [Benchmarking results](#benchmarking-results)
  - [Table of contents](#table-of-contents)
  - [Known vcf](#known-vcf)
  - [bench1.0](#bench10)
    - [Run parameters/settings](#run-parameterssettings)
    - [Results](#results)
      - [NIST7035](#nist7035)
        - [human_genomics_pipeline + minimal vcf_annotation_pipeline](#humangenomicspipeline--minimal-vcfannotationpipeline)
        - [parabricks germline pipeline](#parabricks-germline-pipeline)
      - [NIST7086](#nist7086)
        - [human_genomics_pipeline + minimal vcf_annotation_pipeline](#humangenomicspipeline--minimal-vcfannotationpipeline-1)
        - [parabricks germline pipeline](#parabricks-germline-pipeline-1)
  - [bench1.1](#bench11)
    - [Run parameters/settings](#run-parameterssettings-1)
    - [Results](#results-1)
      - [NIST7035](#nist7035-1)
        - [human_genomics_pipeline + minimal vcf_annotation_pipeline](#humangenomicspipeline--minimal-vcfannotationpipeline-2)
        - [parabricks germline pipeline](#parabricks-germline-pipeline-2)
      - [NIST7086](#nist7086-1)
        - [human_genomics_pipeline + minimal vcf_annotation_pipeline](#humangenomicspipeline--minimal-vcfannotationpipeline-3)
        - [parabricks germline pipeline](#parabricks-germline-pipeline-3)

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

*Note. dbsnp_137.hg19.vcf doesn't appear to be available anymore, therefore I used the closest version available*

(see run settings/output for [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline/tree/bench1.0) and [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline/tree/bench1.0) for more detail)

### Results

#### NIST7035

##### human_genomics_pipeline + minimal vcf_annotation_pipeline

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
| INDEL | ALL    | 33935       | 32265    | 1670     | 57706       | 1595     | 23655     | 733   | 694   | 0.950788      | 0.953158         | 0.409923       | 0.951972        |                        |                        | 0.575724774               | 0.530085807               |
| INDEL | PASS   | 33935       | 31582    | 2353     | 56457       | 1505     | 23181     | 725   | 628   | 0.930662      | 0.954772         | 0.410596       | 0.942563        |                        |                        | 0.575724774               | 0.498879283               |
| SNP   | ALL    | 219952      | 216841   | 3111     | 411540      | 190689   | 3974      | 1429  | 6845  | 0.985856      | 0.532127         | 0.009656       | 0.691181        | 2.032894465            | 1.423485638            | 0.577210523               | 0.362133473               |
| SNP   | PASS   | 219952      | 215192   | 4760     | 408910      | 189812   | 3870      | 1399  | 6523  | 0.978359      | 0.531375         | 0.009464       | 0.688698        | 2.032894465            | 1.422810364            | 0.577210523               | 0.35380422                |

##### parabricks germline pipeline

- Compared with bcftool isec

- Compared with hap.py

| Type  | Filter | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | QUERY.UNK | FP.gt | FP.al | METRIC.Recall | METRIC.Precision | METRIC.Frac_NA | METRIC.F1_Score | TRUTH.TOTAL.TiTv_ratio | QUERY.TOTAL.TiTv_ratio | TRUTH.TOTAL.het_hom_ratio | QUERY.TOTAL.het_hom_ratio |
|-------|--------|-------------|----------|----------|-------------|----------|-----------|-------|-------|---------------|------------------|----------------|-----------------|------------------------|------------------------|---------------------------|---------------------------|
| INDEL | ALL    | 33935       | 32805    | 1130     | 54368       | 1236     | 20169     | 598   | 467   | 0.966701      | 0.963859         | 0.370972       | 0.965278        |                        |                        | 0.575724774               | 0.450235194               |
| INDEL | PASS   | 33935       | 32805    | 1130     | 54368       | 1236     | 20169     | 598   | 467   | 0.966701      | 0.963859         | 0.370972       | 0.965278        |                        |                        | 0.575724774               | 0.450235194               |
| SNP   | ALL    | 219952      | 218293   | 1659     | 414135      | 195815   | 0         | 1550  | 7445  | 0.992457      | 0.527171         | 0              | 0.688583        | 2.032894465            | 1.42376204             | 0.577210523               | 0.360109056               |
| SNP   | PASS   | 219952      | 218293   | 1659     | 414135      | 195815   | 0         | 1550  | 7445  | 0.992457      | 0.527171         | 0              | 0.688583        | 2.032894465            | 1.42376204             | 0.577210523               | 0.360109056               |

#### NIST7086

##### human_genomics_pipeline + minimal vcf_annotation_pipeline

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

- Compared with hap.py

##### parabricks germline pipeline

## bench1.1

### Run parameters/settings

- **Aim:** Benchmarking against against the Genome In A Bottle (GIAB) sample [NA12878](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/) (both NIST7035 and NIST7086)
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
  - 1D model with pre-trained architecture (for rule gatk CNNScoreVariants)

*Note. dbsnp_137.hg19.vcf doesn't appear to be available anymore, therefore I used the closest version available*

(see run settings/output for [human_genomics_pipeline](https://github.com/ESR-NZ/human_genomics_pipeline/tree/bench1.0) and [vcf_annotation_pipeline](https://github.com/ESR-NZ/vcf_annotation_pipeline/tree/bench1.0) for more detail)

### Results

#### NIST7035

##### human_genomics_pipeline + minimal vcf_annotation_pipeline

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

- Compared with hap.py

##### parabricks germline pipeline

#### NIST7086

##### human_genomics_pipeline + minimal vcf_annotation_pipeline

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

- Compared with hap.py

##### parabricks germline pipeline

- Compared with bcftools isec
  
- Compared with hap.py
