# Create a configuration files for scout

Created: 2020-05-08 13:00:23
Last modified: 2020/05/08 13:31:28

This document describes how to create config files for loading cases into [scout](http://www.clinicalgenomics.se/scout/)

## Table of contents

- [Create a configuration files for scout](#create-a-configuration-files-for-scout)
  - [Table of contents](#table-of-contents)
  - [Minimal example](#minimal-example)
    - [Single sample](#single-sample)
    - [Cohort](#cohort)
  - [Notes](#notes)

## Minimal example

Both single samples and cohorts can be put into scout

### Single sample

```yaml
---

owner: cust001

family: 'internal_id_2'
samples:
  - analysis_type: wes
    sample_id: ADM1059A1
    father: 0
    mother: 0
    sample_name: NA12882
    phenotype: unknown
    sex: unknown
    expected_coverage: 30
    bam_path: /scout/demo/reduced_mt.bam

vcf_snv: /scout/demo/643595.clinical.vcf.gz
```

### Cohort

```yaml
---

owner: cust001

family: 'internal_id_2'
samples:
  - analysis_type: wes
    sample_id: ADM1059A1
    father: 0
    mother: 0
    sample_name: NA12882
    phenotype: unknown
    sex: unknown
    expected_coverage: 30
    bam_path: /scout/demo/reduced_mt.bam

  - analysis_type: wes
    sample_id: ADM1059A2
    father: 0
    mother: 0
    sample_name: NA12882
    phenotype: unknown
    sex: unknown
    expected_coverage: 30
    bam_path: /scout/demo/reduced_mt.bam

  - analysis_type: wes
    sample_id: ADM1059A3
    father: 0
    mother: 0
    sample_name: NA12882
    phenotype: unknown
    sex: unknown
    expected_coverage: 30
    bam_path: /scout/demo/reduced_mt.bam

vcf_snv: /scout/demo/643595.clinical.vcf.gz
```

## Notes

- The 'sample_id' value needs to match the sample id/s in your vcf file

These can be found in the columns of your vcf file, for example for an example cohort vcf file provided by scout (643595.clinical.vcf.gz, found in the demo folder), there are three samples that were jointly called:

```bash
#CHROM  POS     ID      REF     ALT     QUAL    FILTER         INFO    FORMAT  ADM1059A1       ADM1059A2       ADM1059A3
1       7847367 .       T       TC      14      IndelGap        AC=3;AF=0.5;AN=6;ANN=TC|sequence_feature|LOW|PER3|ENSG00000049246|domain:PAS_1|ENST00000377532|protein_coding|3/21|c.390+471_390+472insC||||||,TC|sequence_feature|LOW|PER3|ENSG00000049246|domain:PAS_1|ENST00000361923|protein_coding|3/21|c.390+471_390+472insC||||||,TC|upstream_gene_variant|MODIFIER|PER3|ENSG00000049246|transcript|ENST00000602883|retained_intron||n.-15_-14insC|||||14|,TC|intron_variant|MODIFIER|PER3|ENSG00000049246|transcript|ENST00000377532|protein_coding|3/20|c.390+471_390+472insC||||||,TC|intron_variant|MODIFIER|PER3|ENSG00000049246|transcript|ENST00000377541|protein_coding|4/9|c.390+471_390+472insC||||||,TC|intron_variant|MODIFIER|PER3|ENSG00000049246|transcript|ENST00000361923|protein_coding|3/20|c.390+471_390+472insC||||||,TC|intron_variant|MODIFIER|PER3|ENSG00000049246|transcript|ENST00000473653|retained_intron|3/3|n.229+61_229+62insC||||||;CSQ=C|intron_variant|MODIFIER|PER3|ENSG00000049246|Transcript|ENST00000361923|protein_coding||3/20|ENST00000361923.2:c.390+471_390+472insC|||||||||1||HGNC|8847|YES|||CCDS89.1|ENSP00000355031|P56645|Q8TAR6&B4DR65&A2I2N5|UPI0000167B1D|||||||||0.978||||,C|intron_variant|MODIFIER|PER3|ENSG00000049246|Transcript|ENST00000377532|protein_coding||3/20|ENST00000377532.3:c.390+471_390+472insC|||||||||1||HGNC|8847|||||ENSP00000366755|P56645||UPI00003664CA|||||||||0.978||||,C|intron_variant|MODIFIER|PER3|ENSG00000049246|Transcript|ENST00000377541|protein_coding||4/9|ENST00000377541.1:c.390+471_390+472insC|||||||||1||HGNC|8847|||||ENSP00000366764||Q8TAR6|UPI0000070D9E|||||||||0.978||||,C|intron_variant&non_coding_transcript_variant|MODIFIER|PER3|ENSG00000049246|Transcript|ENST00000473653|retained_intron||3/3|ENST00000473653.1:n.229+61_229+62insC|||||||||1||HGNC|8847|||||||||||||||||0.978||||,C|upstream_gene_variant|MODIFIER|PER3|ENSG00000049246|Transcript|ENST00000602883|retained_intron|||||||||||14|1||HGNC|8847|||||||||||||||||0.978||||:ENST00000361923:intron_variant:ENST00000377532:intron_variant:ENST00000377541:intron_variant:ENST00000473653:intron_variant&non_coding_transcript_variant:ENST00000602883:upstream_gene_variant;DP=124;DP4=2,1,0,4;HOB=0.166667;ICB=1;IDV=1;IMF=0.0212766;INDEL;MQ=49;MQ0F=0;MQSB=0.8;RankResult=0|0|-12|0|1|0|0|4|0|0;RankScore=internal_id_2:-7;SGB=-0.426121;VDB=0.388649;most_severe_consequence=8847:C|non_coding_transcript_variant;set=FilteredInAll      GT:AD:DV:PL      0/0:2,0:0:0,6,31        0/1:1,2:2:20,0,8        1/1:0,2:2:31,6,0
1       7879627 rs10462018      C       G       3380.38 PASS    AB=0;ABP=0;AC=0;AF=0;AN=6;ANN=G|intron_variant|MODIFIER|PER3|ENSG00000049246|transcript|ENST00000377532|protein_coding|13/20|c.1658+147C>G||||||,G|intron_variant|MODIFIER|RP3-467L1.4|ENSG00000236266|transcript|ENST00000451646|antisense|1/2|n.239+7537G>C||||||,G|intron_variant|MODIFIER|PER3|ENSG00000049246|transcript|ENST00000361923|protein_coding|13/20|c.1634+147C>G||||||;AO=0;BQB=0.43525;BaseQRankSum=4.42;CIGAR=1X;CSQ=G|intron_variant|MODIFIER|PER3|ENSG00000049246|Transcript|ENST00000361923|protein_coding||13/20|ENST00000361923.2:c.1634+147C>G|||||||||1||HGNC|8847|YES|||CCDS89.1|ENSP00000355031|P56645|Q8TAR6&B4DR65&A2I2N5|UPI0000167B1D|||||||||0.978||||,G|intron_variant|MODIFIER|PER3|ENSG00000049246|Transcript|ENST00000377532|protein_coding||13/20|ENST00000377532.3:c.1658+147C>G|||||||||1||HGNC|8847|||||ENSP00000366755|P56645||UPI00003664CA|||||||||0.978||||:ENST00000361923:intron_variant:ENST00000377532:intron_variant;ClippingRankSum=0;DB;DP=276;DP4=18,35,4,9;DPB=69;DPRA=0;EPP=0;EPPR=13.3047;ExcessHet=4.449;FS=9.052;GTI=0;HOB=0.0555556;ICB=0.128205;InbreedingCoeff=-0.1017;LEN=0;MEANALT=0;MQ0=0;MQ0F=0;MQB=0.000396205;MQM=0;MQMR=60;MQRankSum=0.252;MQSB=0.724672;NEGATIVE_TRAIN_SITE;NS=3;NUMALT=1;ODDS=19.4275;OLD_MULTIALLELIC=1:7879627:C/T/G;Obs=6;PAIRED=0;PAIREDR=1;PAO=0;POSITIVE_TRAIN_SITE;PQA=0;PQR=0;PRO=0;QA=0;QD=15.58;QR=1950;RO=54;RPB=0.915887;RPL=0;RPP=0;RPPR=95.6598;RPR=0;RUN=0;RankResult=0|3|-12|0|1|0|0|2|0|0;RankScore=internal_id_2:-6;ReadPosRankSum=0.665;SAF=0;SAP=0;SAR=0;SGB=12.2736;SOR=1.754;SPIDEX=-0.035;SRF=18;SRP=16.0391;SRR=36;SWEREFAC_Hemi=0;SWEREFAC_Het=6;SWEREFAC_Hom=0;SWEREFAF=0.003;TYPE=snp;VDB=0.0388906;VQSLOD=1.5;culprit=MQ;most_severe_consequence=8847:G|intron_variant;set=gatk-filterInsamtools-freebayes;technology.ILLUMINA=0  GT:DP:GQ:PP     0/.:22:99:499,249       0/0:23:48:0,648        0/0:25:78:0,1027
1       7889972 rs57875989      GAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGCCTCCCATGA G       33161   LOW_VQSLOD      AC=2;AF=0.333;AN=6;ANN=G|inframe_deletion|MODERATE|PER3|ENSG00000049246|transcript|ENST00000377532|protein_coding|18/21|c.3046_3099delGCTCTGTCCACAGGATCGCCTCCCATGAAGAATCCATCCCATCCTACTGCCAGC|p.Ala1016_Ser1033del|3270/6279|3046/3633|1016/1210||INFO_REALIGN_3_PRIME,G|inframe_deletion|MODERATE|PER3|ENSG00000049246|transcript|ENST00000361923|protein_coding|18/21|c.3019_3072delGCTCTGTCCACAGGATCGCCTCCCATGAAGAATCCATCCCATCCTACTGCCAGC|p.Ala1007_Ser1024del|3194/6203|3019/3606|1007/1201||INFO_REALIGN_3_PRIME,G|sequence_feature|LOW|PER3|ENSG00000049246|compositionally_biased_region:Ser-rich|ENST00000377532|protein_coding|19/21|c.2966_3019delAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGCCTCCCATGA||||||,G|sequence_feature|LOW|PER3|ENSG00000049246|repeat:1|ENST00000377532|protein_coding|18/21|c.2966_3019delAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGCCTCCCATGA||||||,G|sequence_feature|LOW|PER3|ENSG00000049246|region_of_interest:5_X_18_AA_tandem_repeats_of_S-[HP]-[AP]-T-[AT]-[GST]-[ATV]-L-S-[MT]-G-[LS]-P-P-[MRS]-[EKR]-[NST]-P|ENST00000377532|protein_coding|18/21|c.2966_3019delAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGCCTCCCATGA||||||,G|sequence_feature|LOW|PER3|ENSG00000049246|repeat:2|ENST00000377532|protein_coding|18/21|c.2966_3019delAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGCCTCCCATGA||||||,G|sequence_feature|LOW|PER3|ENSG00000049246|compositionally_biased_region:Ser-rich|ENST00000361923|protein_coding|18/21|c.2939_2992delAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGCCTCCCATGA||||||,G|sequence_feature|LOW|PER3|ENSG00000049246|repeat:1|ENST00000361923|protein_coding|18/21|c.2939_2992delAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGCCTCCCATGA||||||,G|sequence_feature|LOW|PER3|ENSG00000049246|region_of_interest:5_X_18_AA_tandem_repeats_of_S-[HP]-[AP]-T-[AT]-[GST]-[ATV]-L-S-[MT]-G-[LS]-P-P-[MRS]-[EKR]-[NST]-P|ENST00000361923|protein_coding|18/21|c.2939_2992delAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGCCTCCCATGA||||||,G|sequence_feature|LOW|PER3|ENSG00000049246|repeat:2|ENST00000361923|protein_coding|18/21|c.2939_2992delAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGCCTCCCATGA||||||,G|upstream_gene_variant|MODIFIER|RP3-467L1.4|ENSG00000236266|transcript|ENST00000451646|antisense||n.-2624_-2571delTCATGGGAGGCGATCCTGTGGACAGAGCGCTGGCAGTAGGATGGGATGGATTCT|||||2571|;BaseQRankSum=1.43;CADD=2.444;CSQ=-|inframe_deletion|MODERATE|PER3|ENSG00000049246|Transcript|ENST00000361923|protein_coding|18/21||ENST00000361923.2:c.3019_3072delGCTCTGTCCACAGGATCGCCTCCCATGAAGAATCCATCCCATCCTACTGCCAGC|ENSP00000355031.2:p.Ala1007_Ser1024del|3114-3167|2939-2992|980-998|ENPSHPTASALSTGSPPMK/E|gAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGCCTCCCATGAag/gag|||1||HGNC|8847|YES|||CCDS89.1|ENSP00000355031|P56645|Q8TAR6&B4DR65&A2I2N5|UPI0000167B1D|||Pfam_domain:PF12114&hmmpanther:PTHR11269&hmmpanther:PTHR11269:SF13|80|||||0.978||||,-|inframe_deletion|MODERATE|PER3|ENSG00000049246|Transcript|ENST00000377532|protein_coding|18/21||ENST00000377532.3:c.3046_3099delGCTCTGTCCACAGGATCGCCTCCCATGAAGAATCCATCCCATCCTACTGCCAGC|ENSP00000366755.3:p.Ala1016_Ser1033del|3190-3243|2966-3019|989-1007|ENPSHPTASALSTGSPPMK/E|gAGAATCCATCCCATCCTACTGCCAGCGCTCTGTCCACAGGATCGCCTCCCATGAag/gag|||1||HGNC|8847|||||ENSP00000366755|P56645||UPI00003664CA|||hmmpanther:PTHR11269&hmmpanther:PTHR11269:SF13&Pfam_domain:PF12114|80|||||0.978||||:ENST00000361923:inframe_deletion:ENST00000377532:inframe_deletion;ClippingRankSum=0;DB;DP=43;EXACAF=0.017;EXACMAX_AF=0.0160871302957634;ExcessHet=0;FS=0;InbreedingCoeff=0.9515;MQ=7.38;MQ0=0;MQRankSum=0.427;NEGATIVE_TRAIN_SITE;Obs=114;PG=0,0,0;QD=30.91;RankResult=0|0|-12|0|5|0|0|-12|0|0;RankScore=internal_id_2:-19;ReadPosRankSum=-0.524;SOR=0.551;VQSLOD=-3.837;culprit=MQ;most_severe_consequence=8847:-|inframe_deletion;set=FilteredInAll GT:AD:DP:GQ:JL:JP:PL:PP 0/0:30,0:30:75:5:40:0,76,1174:0,75,1173        0/1:9,2:11:42:5:40:43,0,362:42,0,422    0/1:2,0:2:47:5:40:4,0,77:47,0,77
1       7890024 .       T       G       15.99   LowQual AC=1;AF=0.167;AN=6;ANN=G|missense_variant|MODERATE|PER3|ENSG00000049246|transcript|ENST00000377532|protein_coding|18/21|c.3017T>G|p.Met1006Arg|3241/6279|3017/3633|1006/1210||,G|missense_variant|MODERATE|PER3|ENSG00000049246|transcript|ENST00000361923|protein_coding|18/21|c.2990T>G|p.Met997Arg|3165/6203|2990/3606|997/1201||,G|sequence_feature|LOW|PER3|ENSG00000049246|compositionally_biased_region:Ser-rich|ENST00000377532|protein_coding|19/21|c.3017T>G||||||,G|sequence_feature|LOW|PER3|ENSG00000049246|region_of_interest:5_X_18_AA_tandem_repeats_of_S-[HP]-[AP]-T-[AT]-[GST]-[ATV]-L-S-[MT]-G-[LS]-P-P-[MRS]-[EKR]-[NST]-P|ENST00000377532|protein_coding|18/21|c.3017T>G||||||,G|sequence_feature|LOW|PER3|ENSG00000049246|repeat:2|ENST00000377532|protein_coding|18/21|c.3017T>G||||||,G|sequence_feature|LOW|PER3|ENSG00000049246|compositionally_biased_region:Ser-rich|ENST00000361923|protein_coding|18/21|c.2990T>G||||||,G|sequence_feature|LOW|PER3|ENSG00000049246|region_of_interest:5_X_18_AA_tandem_repeats_of_S-[HP]-[AP]-T-[AT]-[GST]-[ATV]-L-S-[MT]-G-[LS]-P-P-[MRS]-[EKR]-[NST]-P|ENST00000361923|protein_coding|18/21|c.2990T>G||||||,G|sequence_feature|LOW|PER3|ENSG00000049246|repeat:2|ENST00000361923|protein_coding|18/21|c.2990T>G||||||,G|upstream_gene_variant|MODIFIER|RP3-467L1.4|ENSG00000236266|transcript|ENST00000451646|antisense||n.-2622A>C|||||2622|;BQB=0.76338;CADD=0.320;CSQ=G|missense_variant|MODERATE|PER3|ENSG00000049246|Transcript|ENST00000361923|protein_coding|18/21||ENST00000361923.2:c.2990T>G|ENSP00000355031.2:p.Met997Arg|3165|2990|997|M/R|aTg/aGg|||1||HGNC|8847|YES|||CCDS89.1|ENSP00000355031|P56645|Q8TAR6&B4DR65&A2I2N5|UPI0000167B1D|tolerated_low_confidence|unknown|Pfam_domain:PF12114&hmmpanther:PTHR11269&hmmpanther:PTHR11269:SF13||||||0.978||||,G|missense_variant|MODERATE|PER3|ENSG00000049246|Transcript|ENST00000377532|protein_coding|18/21||ENST00000377532.3:c.3017T>G|ENSP00000366755.3:p.Met1006Arg|3241|3017|1006|M/R|aTg/aGg|||1||HGNC|8847|||||ENSP00000366755|P56645||UPI00003664CA|tolerated_low_confidence|unknown|hmmpanther:PTHR11269&hmmpanther:PTHR11269:SF13&Pfam_domain:PF12114||||||0.978||||:ENST00000361923:missense_variant:ENST00000377532:missense_variant;DP=54;DP4=18,22,0,4;EXACAF=8.622e-06;EXACMAX_AF=6.19655471557814e-05;HOB=0.0555556;ICB=0.128205;MQ=45;MQ0F=0;MQB=0.319553;MQSB=0.969078;Obs=47;RPB=0.0168512;RankResult=5|0|-12|0|5|0|0|-12|0|0;RankScore=internal_id_2:-14;SGB=-2.50655;SPIDEX=-2.44;VDB=0.0221621;dbNSFP_phastCons100way_vertebrate=0.004000;dbNSFP_phastCons46way_primate=0.013000;dbNSFP_phyloP100way_vertebrate=-1.043000;dbNSFP_phyloP46way_primate=-1.459000;most_severe_consequence=8847:G|missense_variant;phastCons100way_vertebrate_prediction_term=NotConserved;phyloP100way_vertebrate_prediction_term=NotConserved;set=FilteredInAll        GT:AD:DV:PL     0/0:14,0:0:0,42,255     0/0:19,1:1:0,30,255     0/1:7,3:3:51,0,160
[...]
```

Their sample id's are:

- ADM1059A1
- ADM1059A2
- ADM1059A3

*Note: a single sample analysis will return a vcf file with a single column/sample id*