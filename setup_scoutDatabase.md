# Create a scout database

Created: 2020/03/11 11:25:43
Last modified: 2020/03/11 11:28:02

- **Aim:** Create a [scout](http://www.clinicalgenomics.se/scout/) database
- **Prerequisite software:** [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html)
- **OS:** Ubuntu 16.04 (Wintermute - research server)

These notes are adapted/extended from documentation written by Miles Benton

## Table of contents

- [Create a scout database](#create-a-scout-database)
  - [Table of contents](#table-of-contents)
  - [Setup](#setup)
  - [Create database](#create-database)
  - [Populate database](#populate-database)
  - [Index the database](#index-the-database)
  - [Create a case](#create-a-case)
    - [Load cases from a config file](#load-cases-from-a-config-file)
    - [Load cases from the command line](#load-cases-from-the-command-line)
  - [Testing](#testing)
  - [Create madeline pedigree drawing](#create-madeline-pedigree-drawing)
    - [Peddy test](#peddy-test)
  - [OMIM access](#omim-access)

## Setup

Make sure:

- you have activated your conda environment you installed scout and all it's dependencies
- a mongod process is running

(see documentation: installation_scout.md)

## Create database

Create a mongoDB database

```bash
scout -db test-database setup
```

Load a new institute :confused:

```bash
scout load institute -i ESR001 --display-name ESR
```

Add users

```bash
scout -db test-database load user -i ESR001 -u "Leah Kemp" -m leahmhkemp@gmail.com --admin
scout -db test-database load user -i ESR001 -u "Miles Benton" -m miles.benton84@gmail.com --admin
scout -db test-database load user -i ESR001 -u "Joep de ligt" -m joepio@gmail.com --admin
```

To view the underlying mongod databases and delete them

```bash
mongo # Open a mongo shell
show dbs
use scout-demo
db.dropDatabase()
quit() # Exit mongo shell
```

See the institute and users

```bash
scout -db test-database view institutes
scout -db test-database view users
```

## Populate database

```bash
scout -db test-database update genes --build 37 --api-key EU_PX7n2SDyx5CnQII5xyQ
scout -db test-database load exons --build 37
scout -db test-database update diseases --api-key EU_PX7n2SDyx5CnQII5xyQ
scout -db test-database update hpo # Didn't seem to work
# Be careful how often you run 'update genes' and 'update diseases' since you can only 
# download these a limited number of times
scout -db test-database load panel --panel-app --institute ESR001
```

Note. extra steps are required to download genes and exons for the 38 build of the reference genome

## Index the database

```bash
scout index # Didn't seem to change anything
```

## Create a case

For example I created this yaml file saved as 'trio_4001.config.yaml'

```yaml
---

owner: ESR001

family: 'internal_id_4'
family: 'Trio_4001'
samples:
  - analysis_type: wes
    sample_id: WDHB4001
    father: WDHB4003
    mother: WDHB4004
    sample_name: WDHB4002
    phenotype: affected
    sex: male
    tissue_type: blood
    expected_coverage: 30
    bam_path: /store/lkemp/exome_project/data/mapped/4001_bwa_recal.cram

  - analysis_type: wes
    sample_id: WDHB4003
    father: 0
    mother: 0
    sample_name: WDHB4003
    phenotype: unaffected
    sex: female
    tissue_type: blood
    expected_coverage: 30
    bam_path: /store/lkemp/exome_project/data/mapped/4003_bwa_recal.cram

  - analysis_type: wes
    sample_id: WDHB4004
    father: 0
    mother: 0
    sample_name: WDHB4004
    phenotype: unaffected
    sex: male
    tissue_type: blood
    expected_coverage: 30
    bam_path: /store/lkemp/exome_project/data/mapped/4004_bwa_recal.cram

#vcf_snv: /data/WCHP-Clinical-Genetics/Wellington/vcf/annot/vcf/cohort_WDHB3988.inherit.filtered.recalibrated_dbNSFP_VEP_clean.vcf.gz
vcf_snv: /store/lkemp/exome_project/data/vcf/annot/4001/trio_WDHB4001.inherit.annot.filtered.recalibrated.clean.dbnsfp.VEP.mod.probandfilter.CADD.vcf

peddy_ped: /store/lkemp/exome_project/data/pedigrees/madeline/trio_4001/Trio_4001.peddy.ped
peddy_ped_check: /store/lkemp/exome_project/data/pedigrees/madeline/trio_4001/Trio_4001.ped_check.csv
peddy_sex_check: /store/lkemp/exome_project/data/pedigrees/madeline/trio_4001/Trio_4001.sex_check.csv

madeline: /store/lkemp/exome_project/data/pedigrees/madeline/trio_4001/Trio_4001_pedigree.xml

# default_gene_panels: [558aa423bb5a16630e15b63c,55b605f722c1fc05fd2345af]
# gene_panels: [558aa423bb5a16630e15b63c,55b605f722c1fc05fd2345af]

# this allows you to attach reports and notes etc
#delivery_report: /data/WCHP-Clinical-Genetics/Wellington/scout_test/delivery_report.html

# meta data
analysis_date: 2019-11-04 14:35:20
```

Run scout!

```bash
scout -db test-database serve
```

### Load cases from a config file

```bash
scout -db test-database load case /store/lkemp/exome_project/scout/trio_4001.config.yaml
```

### Load cases from the command line

```bash
scout load case --owner ESR001 \
--ped /data/WCHP-Clinical-Genetics/Wellington/pedigrees/madeline/trio_1_3988/pedigree_3988_tab.ped \
--vcf /data/WCHP-Clinical-Genetics/Wellington/vcf/annot/vcf/cohort_WDHB3988.inherit.filtered.recalibrated_dbNSFP_VEP_clean.vcf.gz \
--peddy-ped /data/WCHP-Clinical-Genetics/Wellington/pedigrees/madeline/trio_1_3988/Trio_1_3988.peddy.ped \
--peddy-sex /data/WCHP-Clinical-Genetics/Wellington/pedigrees/madeline/trio_1_3988/Trio_1_3988.sex_check.csv \
--peddy-check /data/WCHP-Clinical-Genetics/Wellington/pedigrees/madeline/trio_1_3988/Trio_1_3988.ped_check.csv

# scout load case --owner ESR001 -u \
# --ped /data/WCHP-Clinical-Genetics/Wellington/pedigrees/madeline/trio_1_3988/pedigree_3988_tab.ped \
# --peddy-ped /data/WCHP-Clinical-Genetics/Wellington/pedigrees/madeline/trio_1_3988/Trio_1_3988.peddy.ped \
# --peddy-sex /data/WCHP-Clinical-Genetics/Wellington/pedigrees/madeline/trio_1_3988/Trio_1_3988.sex_check.csv \
# --peddy-check /data/WCHP-Clinical-Genetics/Wellington/pedigrees/madeline/trio_1_3988/Trio_1_3988.ped_check.csv
```

## Testing

Deleting 'old' cases

```bash
scout delete case -i ESR001 -c Trio_1
```

## Create madeline pedigree drawing  

```bash
#3988
ped_parser ./pedigrees/pedigree_3988.ped --to_madeline --outfile ./pedigrees/madeline/trio_1_3988/WDHB_trio1_3988.ped.data
madeline2 -L "IndividualId Affected" ./pedigrees/madeline/trio_1_3988/WDHB_trio1_3988.ped.data --outputext xml --color
#4002
ped_parser ./pedigrees/trio_4002.ped --to_madeline --outfile ./pedigrees/madeline/trio_4002/WDHB_trio_4002.ped.data
madeline2 -L "IndividualId Affected" ./pedigrees/madeline/trio_4002/WDHB_trio_4002.ped.data --outputext xml --color
peddy -p 8 --plot --prefix Trio_4002 ../../../vcf/annot/4002/trio_WDHB4002.inherit.annot.filtered.recalibrated.clean.vcf.gz ./trio_4002_tab.ped
#4001
ped_parser ./pedigrees/trio_4001.ped --to_madeline --outfile ./pedigrees/madeline/trio_4001/WDHB_trio_4001.ped.data
madeline2 -L "IndividualId Affected" ./pedigrees/madeline/trio_4001/WDHB_trio_4001.ped.data --outputext xml --color
peddy -p 8 --plot --prefix Trio_4001 ../../../vcf/annot/4002/trio_WDHB4001.inherit.annot.filtered.recalibrated.clean.vcf.gz ./trio_4001_tab.ped
```

### Peddy test

Run in /data/WCHP-Clinical-Genetics/Wellington/pedigrees/madeline/trio_1_3988

```bash
peddy -p 8 --plot --prefix Trio_1_3988 ../../../vcf/annot/vcf/cohort_WDHB3988.inherit.filtered.recalibrated_dbNSFP_VEP_clean.vcf.gz ./pedigree_3988_tab.ped
```

## OMIM access

You can get each file you have access to using the URLs listed below:
https://omim.org/static/omim/data/mim2gene.txt
https://data.omim.org/downloads/EU_PX7n2SDyx5CnQII5xyQ/mimTitles.txt
https://data.omim.org/downloads/EU_PX7n2SDyx5CnQII5xyQ/genemap2.txt
https://data.omim.org/downloads/EU_PX7n2SDyx5CnQII5xyQ/morbidmap.txt

You can access the API using the following API key:
[API Key](EU_PX7n2SDyx5CnQII5xyQ)
[API Host](api.omim.org)
[API Base URL](https://api.omim.org/api)
[API Web Interface](https://api.omim.org/api/html/index.html)
[API Documentation](https://omim.org/help/api(scout_env))