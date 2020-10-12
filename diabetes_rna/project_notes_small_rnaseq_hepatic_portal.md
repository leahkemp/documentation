# Project notes - small RNAseq - hepatic portal

Created: 2020/10/13 12:21:30
Last modified: 2020/10/13 13:06:37

- **Aim:** Evaluate the current pipelines available for processing RNA-seq data. This will help us decide if there is a pipeline currently available for our use, one we could adapt, or if we will need to create an RNA-seq pipeline from scratch

## Table of contents

- [Project notes - small RNAseq - hepatic portal](#project-notes---small-rnaseq---hepatic-portal)
  - [Table of contents](#table-of-contents)
  - [Data](#data)
  - [Research questions](#research-questions)
  - [Other notes](#other-notes)
  - [Analysis](#analysis)
  - [Wider context](#wider-context)

## Data

- 60 samples
- 30 type 2 diabetic and 30 *not* type 2 diabetic
- HiSeq
- 100 bp
- Single end
- All female
- Mostly European
- Hepatic portal serum (drains from gut/lymph before going to liver)
- Wakefield hospital - private hospital

## Research questions

- Look for any differences between type 2 diabetic and *not* type 2 diabetic
- Explore any associations of type 2 diabetes with the other measurements taken (eg. preOp_G_Fasting)
- Initial exploratory analysis to inform which further sequencing/analyses/research to undertake

## Other notes

- There are likely a number of routes to diabetes
- More info at http://leviathan/human_genomics/ (on ESR's network)

## Analysis

- Looking at using [smrnaseq](https://github.com/nf-core/smrnaseq) for the initial data analysis/qc
- Will look at/play around with the other tools available for analyses to follow on from [smrnaseq](https://github.com/nf-core/smrnaseq)
- Likely use [limma](https://www.rdocumentation.org/packages/limma/versions/3.28.14) and [Deseq2](https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3) (in R) (use both packages since they have different strengths)
- Interested in analysing isomers (difference in RNA sequences between individuals)
  - Note. some software calls differences in length as different isomers, however differences in length can be an artifact of sequencing and/or trimming
  - It's a good idea to only look at isomers that are different in the sequences
- Mature vs. pre-curser (hairpin)
- Keep unmapped reads - it's likely that it will map to bacteria as gut microbiota could be included in the sample (software to map unmapped reads to bacteria: [kraken](https://bio.tools/kraken), [centrifuge](https://bio.tools/Centrifuge))

- QC and sanity checks
  - Good idea to check the lengths of the RNA's between individuals as a sanity check
  - Plot read counts (should be ~12 million reads per sample, 2 million would be very low and data we might not want to analyse)
  - Mapping efficiency (% mapped and % not mapped)

- Trimming:
  - Trimming is quite important for small RNA data because you tend to get alot more adapters sequenced since sequencing won't drop out before the end of the sequence (like it would with longer fragments)
  - [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) and trim galore are generally robust trimming tools
  - I remember the RNA Carpentries mentioning they preferred cutadapt - will have to look into what it's benefit was
  - Automatic trimming often works well, but there are many times where it is insufficient
  - Good idea to apply a hard cap (minimum length) of around 17, this is because smaller smRNA's will be more likely to map to many locations
  - Alarm bells should ring if there is no trimming, or over trimming (often comes up as a problem)

- Other
  - Don't be alarmed if there are alot of sequences not mapped, you are extracting only the small non-coding RNA's out of the whole pool of small RNA's
  - [tximport](https://www.rdocumentation.org/packages/tximport/versions/1.0.3) is a good package for importing the bams (output by [smrnaseq](https://github.com/nf-core/smrnaseq)) into R for further analysis

## Wider context

- Have access to an extensive tissue bank that can be used for further research
- 400 samples collected during gastric bypass surgeries, eg.
  - Abdominal adipose tissue
  - Subcutaneous adipose tissues
  - Blood
  - Jejunum
  - Pancreas
  - Individuals with and without type 2 diabetes

- Option for pre-post surgery analysis
  