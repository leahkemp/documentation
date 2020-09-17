# Test somalier as replacement for peddy

Created: 2020/08/31 11:36:14
Last modified: 2020/09/17 17:20:54

- **Aim:** Try/test [somalier](https://github.com/brentp/somalier) and see if the outputs can be integrated with [scout](https://github.com/Clinical-Genomics/scout)
- **Prerequisite software:** [singularity 2.5.2](https://singularity.lbl.gov/)
- **OS:** Ubuntu 16.04 (Wintermute - research server)

In our whole genome sequencing (WGS) workflow (fastq to vcf to preparing data to be ingested into scout) have been using [peddy](https://github.com/brentp/peddy) to evaluate/check the sex/familial relationships for the joint vcf's of our cohort samples. The files that peddy outputs are used and presented in [scout](https://github.com/Clinical-Genomics/scout) to show whether the genomic data supports the familial relationships and sex's defined in the pedigree file. We have had some inconsistency in the results of these checks with peddy depending on how the data is subset (eg. exome capture regions only and whether non-pass variants removed from the vcf). [somalier](https://github.com/brentp/somalier) has been released which ["is a more scalable, faster, replacement for peddy that uses some of the same methods as peddy along with some new ones"](https://github.com/brentp/peddy#fast-pedigreevcf-qc). It also allows you to do these checks based on bam's instead of vcf's. So we want to try it out and see if it's outputs can still be integrated into scout, and therefore included in our pipeline/workflow!

## Table of contents

- [Test somalier as replacement for peddy](#test-somalier-as-replacement-for-peddy)
  - [Table of contents](#table-of-contents)
  - [Try somalier](#try-somalier)
    - [Based on vcf](#based-on-vcf)
    - [Based on bam](#based-on-bam)
  - [Try incorporate somalier output into scout](#try-incorporate-somalier-output-into-scout)
    - [Scout script that parses a peddy.ped file](#scout-script-that-parses-a-peddyped-file)
    - [Scout script that parses a .ped_check.csv file](#scout-script-that-parses-a-ped_checkcsv-file)
    - [Scout script that parses a .ped_check.csv file](#scout-script-that-parses-a-ped_checkcsv-file-1)
  - [Notes](#notes)

## Try somalier

Get data to play with from manual pipeline runs

```bash
cd /store/lkemp/test_somalier/
mkdir based_on_vcf
mkdir based_on_bam
cp /store/lkemp/manual_pipeline_run/WGS_b37/DHB4205_GRCh37_WGS_trio.combined.gpu.genotyped.vqsr.withPosteriors.denovo_ann.dbnsfp.vep.cadd.models.sorted.dbsnp.onlyproband.dualalleles.rankscore.exomes.clean.vcf.gz* .
cp /store/lkemp/manual_pipeline_run/WGS_b37/DHB4205_pedigree.ped .
```

A docker image for somalier is hosted [here](https://hub.docker.com/r/brentp/somalier)

```bash
singularity pull docker://brentp/somalier:latest
```

Get sites files from [here](https://github.com/brentp/somalier/releases/tag/v0.2.12)

```bash
wget https://github.com/brentp/somalier/files/3412453/sites.hg19.vcf.gz
wget https://github.com/brentp/somalier/files/3412454/sites.hg38.nochr.vcf.gz
wget https://github.com/brentp/somalier/files/3412455/sites.GRCh37.vcf.gz
wget https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz
```

Get ancestry files

```bash
wget https://raw.githubusercontent.com/brentp/somalier/master/scripts/ancestry-labels-1kg.tsv
wget https://zenodo.org/record/3479773/files/1kg.somalier.tar.gz?download=1
tar -xf 1kg.somalier.tar.gz?download=1
```

### Based on vcf

Extract sites from vcf

```bash
cd /store/lkemp/test_somalier/based_on_vcf/

BIND_PATHS="/store/"
CONTAINER="/store/lkemp/test_somalier/somalier-latest.simg"

# DHB4205
singularity exec \
-B ${BIND_PATHS} \
${CONTAINER} \
somalier extract \
-d DHB4205_extracted/ \
--sites ../sites.GRCh37.vcf.gz \
-f /store/lkemp/publicData/b37/human_g1k_v37_decoy.fasta \
/store/lkemp/manual_pipeline_run/WGS_b37/DHB4205_GRCh37_WGS_trio.combined.gpu.genotyped.vqsr.withPosteriors.denovo_ann.dbnsfp.vep.cadd.models.sorted.dbsnp.onlyproband.dualalleles.rankscore.vcf.gz

# DHB4177
singularity exec \
-B ${BIND_PATHS} \
${CONTAINER} \
somalier extract \
-d DHB4177_extracted/ \
--sites ../sites.GRCh37.vcf.gz \
-f /store/lkemp/publicData/b37/human_g1k_v37_decoy.fasta \
/store/lkemp/manual_pipeline_run/WGS_b37/DHB4177_GRCh37_WGS_trio.combined.gpu.genotyped.hardfiltered.withPosteriors.denovo_ann.dbnsfp.vep.cadd.models.sorted.dbsnp.onlyproband.dualalleles.rankscore.vcf.gz
```

Calculate relatedness on the extracted data

```bash
# DHB4205
singularity exec \
-B ${BIND_PATHS} \
${CONTAINER} \
somalier relate \
--ped /store/lkemp/manual_pipeline_run/WGS_b37/DHB4205_pedigree.ped \
DHB4205_extracted/*.somalier \
--output-prefix=DHB4205_sex_relate

# DHB4177
singularity exec \
-B ${BIND_PATHS} \
${CONTAINER} \
somalier relate \
--ped /store/lkemp/manual_pipeline_run/WGS_b37/DHB4177_pedigree.ped \
DHB4177_extracted/*.somalier \
--output-prefix=DHB4177_sex_relate
```

Calculate ancestry

```bash
# DHB4205
singularity exec \
-B ${BIND_PATHS} \
${CONTAINER} \
somalier ancestry \
--labels ../ancestry-labels-1kg.tsv \
../1kg-somalier/*.somalier ++ DHB4205_extracted/*.somalier \
--output-prefix=DHB4205

# DHB4177
singularity exec \
-B ${BIND_PATHS} \
${CONTAINER} \
somalier ancestry \
--labels ../ancestry-labels-1kg.tsv \
../1kg-somalier/*.somalier ++ DHB4177_extracted/*.somalier \
--output-prefix=DHB4177
```

### Based on bam

Extract sites from the three bams in the cohort

```bash
cd /store/lkemp/test_somalier/based_on_bam/

# DHB4205
for f in /store/mbenton/WGS_b37/bams/DHB420*.bam; do
    singularity exec \
    -B ${BIND_PATHS} \
    ${CONTAINER} \
    somalier extract \
    -d DHB4205_extracted/ \
    --sites ../sites.GRCh37.vcf.gz \
    -f /store/lkemp/publicData/b37/human_g1k_v37_decoy.fasta \
    $f
done

# DHB4177
for f in /store/mbenton/WGS_b37/bams/DHB4177_GRCh37_WGS.bam /store/mbenton/WGS_b37/bams/DHB4178_GRCh37_WGS.bam /store/mbenton/WGS_b37/bams/DHB4179_GRCh37_WGS.bam ; do
    singularity exec \
    -B ${BIND_PATHS} \
    ${CONTAINER} \
    somalier extract \
    -d DHB4177_extracted/ \
    --sites ../sites.GRCh37.vcf.gz \
    -f /store/lkemp/publicData/b37/human_g1k_v37_decoy.fasta \
    $f
done
```

Calculate relatedness on the extracted data

```bash
# DHB4205
singularity exec \
-B ${BIND_PATHS} \
${CONTAINER} \
somalier relate \
--ped /store/lkemp/manual_pipeline_run/WGS_b37/DHB4205_pedigree.ped \
DHB4205_extracted/*.somalier \
--output-prefix=DHB4205_sex_relate

# DHB4177
singularity exec \
-B ${BIND_PATHS} \
${CONTAINER} \
somalier relate \
--ped /store/lkemp/manual_pipeline_run/WGS_b37/DHB4177_pedigree.ped \
DHB4177_extracted/*.somalier \
--output-prefix=DHB4177_sex_relate
```

Calculate ancestry

```bash
# DHB4205
singularity exec \
-B ${BIND_PATHS} \
${CONTAINER} \
somalier ancestry \
--labels ../ancestry-labels-1kg.tsv \
../1kg-somalier/*.somalier ++ DHB4205_extracted/*.somalier \
--output-prefix=DHB4205

# DHB4177
singularity exec \
-B ${BIND_PATHS} \
${CONTAINER} \
somalier ancestry \
--labels ../ancestry-labels-1kg.tsv \
../1kg-somalier/*.somalier ++ DHB4177_extracted/*.somalier \
--output-prefix=DHB4177
```

## Try incorporate somalier output into scout

Edit the config file which will pass

```bash
cd /store/lkemp/test_somalier/
cp /store/lkemp/manual_pipeline_run/WGS_b37/DHB4205_config.yaml .
cp /store/lkemp/manual_pipeline_run/WGS_b37/DHB4177_config.yaml .
```

Try loading into an existing demo database

```bash
scout -db demo-database load case /store/lkemp/test_somalier/DHB4205_config.yaml
scout -db demo-database load case /store/lkemp/test_somalier/DHB4177_config.yaml
```

Error when trying to load *.pairs.tsv (in place of *.ped_check.csv) by passing it to `peddy_check:`

```bash
(scout_env) lkemp@Wintermute:/store/lkemp/exome_project/scout$ scout -db demo-database load case /store/lkemp/test_somalier/DHB4205_config.yaml
sh: 0: getcwd() failed: No such file or directory
2020-09-17 13:03:44 Wintermute scout.commands.base[12468] INFO Running scout version 4.19
2020-09-17 13:03:44 Wintermute scout.commands.base[12468] DEBUG Debug logging enabled.
2020-09-17 13:03:44 Wintermute scout.commands.load.case[12468] ERROR KEYERROR 'hets_a' missing when loading '/store/lkemp/test_somalier/DHB4205_config.yaml'
2020-09-17 13:03:44 Wintermute scout.commands.load.case[12468] DEBUG Stack trace: Traceback (most recent call last):
  File "/home/lkemp/gunicorn_test/scout/scout/commands/load/case.py", line 95, in case
    peddy_check=peddy_check,
  File "/home/lkemp/gunicorn_test/scout/scout/parse/case.py", line 119, in parse_case_data
    add_peddy_information(config_data)
  File "/home/lkemp/gunicorn_test/scout/scout/parse/case.py", line 237, in add_peddy_information
    for pair_info in parse_peddy_ped_check(file_handle):
  File "/home/lkemp/gunicorn_test/scout/scout/parse/peddy.py", line 65, in parse_peddy_ped_check
    pair_info["hets_a"] = convert_number(pair_info["hets_a"])
KeyError: 'hets_a'

Aborted!
```

Although the column 'hets_a' does exist in *.pairs.tsv


### Scout script that parses a peddy.ped file

```python
def parse_peddy_ped(lines):
    """Parse a peddy.ped file
    
    Args:
        lines(iterable(str))
    
    Returns:
        peddy_ped(list(dict))
    """
    peddy_ped = []
    header = []
    for i, line in enumerate(lines):
        line = line.rstrip()
        if i == 0:
            # Header line
            header = line.lstrip("#").split("\t")
        else:
            ind_info = dict(zip(header, line.split("\t")))

            # PC1/PC2/PC3/PC4: the first 4 values after this sample was
            # projected onto the thousand genomes principle components.
            ind_info["PC1"] = convert_number(ind_info["PC1"])
            ind_info["PC2"] = convert_number(ind_info["PC2"])
            ind_info["PC3"] = convert_number(ind_info["PC3"])
            # ancestry-prediction one of AFR AMR EAS EUR SAS UNKNOWN

            ind_info["het_call_rate"] = convert_number(ind_info["het_call_rate"])

            # idr_baf: inter-decile range (90th percentile - 10th percentile)
            # of b-allele frequency. We make a distribution of all sites of
            # alts / (ref + alts) and then report the difference between the
            # 90th and the 10th percentile.
            # Large values indicated likely sample contamination.
            ind_info["het_idr_baf"] = convert_number(ind_info["het_idr_baf"])

            ind_info["het_mean_depth"] = convert_number(ind_info["het_mean_depth"])

            peddy_ped.append(ind_info)
    return peddy_ped
```

### Scout script that parses a .ped_check.csv file

```python
def parse_peddy_ped_check(lines):
    """Parse a .ped_check.csv file
    
    Args:
        lines(iterable(str))
    
    Returns:
        ped_check(list(dict))
    """
    ped_check = []
    header = []
    for i, line in enumerate(lines):
        line = line.rstrip()
        if i == 0:
            # Header line
            header = line.lstrip("#").split(",")
        else:
            pair_info = dict(zip(header, line.split(",")))

            # the number of sites at which sample_a was heterozygous
            pair_info["hets_a"] = convert_number(pair_info["hets_a"])

            # the number of sites at which sample_b was heterozygous
            pair_info["hets_b"] = convert_number(pair_info["hets_b"])

            # the number of sites at which the 2 samples shared no alleles
            # (should approach 0 for parent-child pairs).
            pair_info["ibs0"] = convert_number(pair_info["ibs0"])

            # the number of sites and which the 2 samples where both
            # hom-ref, both het, or both hom-alt.
            pair_info["ibs2"] = convert_number(pair_info["ibs2"])

            # the number of sites that was used to predict the relatedness.
            pair_info["n"] = convert_number(pair_info["n"])

            # the relatedness reported in the ped file.
            pair_info["rel"] = convert_number(pair_info["rel"])

            # the relatedness reported in the ped file.
            pair_info["pedigree_relatedness"] = convert_number(pair_info["pedigree_relatedness"])

            # difference between the preceding 2 colummns.
            pair_info["rel_difference"] = convert_number(pair_info["rel_difference"])

            # the number of sites at which both samples were hets.
            pair_info["shared_hets"] = convert_number(pair_info["shared_hets"])

            # boolean indicating that this pair is a parent-child pair
            # according to the ped file.
            pair_info["pedigree_parents"] = make_bool(pair_info.get("pedigree_parents"))

            # boolean indicating that this pair is expected to be a parent-child
            # pair according to the ibs0 (< 0.012) calculated from the genotypes.
            pair_info["predicted_parents"] = make_bool(pair_info.get("predicted_parents"))

            # boolean indicating that the preceding 2 columns do not match
            pair_info["parent_error"] = make_bool(pair_info.get("parent_error"))

            #  boolean indicating that rel > 0.75 and ibs0 < 0.012
            pair_info["sample_duplication_error"] = make_bool(
                pair_info.get("sample_duplication_error")
            )

            ped_check.append(pair_info)

    return ped_check
```

### Scout script that parses a .ped_check.csv file

```python

def parse_peddy_sex_check(lines):
    """Parse a .ped_check.csv file
    
    Args:
        lines(iterable(str))
    
    Returns:
        sex_check(list(dict))
    """
    sex_check = []
    header = []
    for i, line in enumerate(lines):
        line = line.rstrip()
        if i == 0:
            # Header line
            header = line.lstrip("#").split(",")
        else:
            ind_info = dict(zip(header, line.split(",")))

            # boolean indicating wether there is a mismatch between X
            # genotypes and ped sex.
            ind_info["error"] = make_bool(ind_info.get("error"))

            # number of homozygous-alternate calls
            ind_info["hom_alt_count"] = convert_number(ind_info["hom_alt_count"])
            # number of homozygous-reference calls
            ind_info["hom_ref_count"] = convert_number(ind_info["hom_ref_count"])
            # number of heterozygote calls
            ind_info["het_count"] = convert_number(ind_info["het_count"])

            # ratio of het_count / hom_alt_count. Low for males, high for females
            ind_info["het_ratio"] = convert_number(ind_info["het_ratio"])

            sex_check.append(ind_info)

    return sex_check
```

## Notes

- Everything (calculating relatedness and ancestry) seems to work well
- Takes longer to extract sites on a bam than vcf (makes sense given the much larger size of bams)
- Calculating ancestry is a little slow and takes a fair amount of resources, probably because of the size of the 1000 genome dataset it's comparing to (but we will be able to control the threading with `export OMP_NUM_THREADS=8`, see [here](https://github.com/brentp/somalier/issues/45))
- It looks like it will be fairly straightforward to include ethnic representation beyond what is included in the [1000 genomes project](https://www.internationalgenome.org/) for the ancestry estimates. It looks as though all we would need to do is run somalier on samples from ethnicities we want to include in the dataset and add them to the ancestry-labels-1kg.tsv file passed to `--labels`
