#!/bin/bash

dir="/NGS/scratch/KSCBIOM/HumanGenomics/exome_project/fastq/"

# Create dummy single data
cp $dir/CH_13BL2450_S1_R1.fastq.gz $dir/dummy01_1.fastq.gz
cp $dir/CH_13BL2450_S1_R2.fastq.gz $dir/dummy01_2.fastq.gz

cp $dir/GO_16BL0892_S4_R1.fastq.gz $dir/dummy02_1.fastq.gz
cp $dir/GO_16BL0892_S4_R2.fastq.gz $dir/dummy02_2.fastq.gz

cp $dir/JC_16BL0361_S6_R1.fastq.gz $dir/dummy03_1.fastq.gz
cp $dir/JC_16BL0361_S6_R2.fastq.gz $dir/dummy03_2.fastq.gz

cp $dir/LJ_17BL1030_S8_R1.fastq.gz $dir/dummy04_1.fastq.gz
cp $dir/LJ_17BL1030_S8_R2.fastq.gz $dir/dummy04_2.fastq.gz

cp $dir/MS_16BL0795_S5_R1.fastq.gz $dir/dummy05_1.fastq.gz
cp $dir/MS_16BL0795_S5_R2.fastq.gz $dir/dummy05_2.fastq.gz

cp $dir/NH_15BL1218_S7_R1.fastq.gz $dir/dummy06_1.fastq.gz
cp $dir/NH_15BL1218_S7_R2.fastq.gz $dir/dummy06_2.fastq.gz

cp $dir/PH_18BL0824_S2_R1.fastq.gz $dir/dummy07_1.fastq.gz
cp $dir/PH_18BL0824_S2_R2.fastq.gz $dir/dummy07_2.fastq.gz

cp $dir/SJ_18BL1034_S3_R1.fastq.gz $dir/dummy08_1.fastq.gz
cp $dir/SJ_18BL1034_S3_R1.fastq.gz $dir/dummy08_2.fastq.gz

# Create dummy cohort data
cp $dir/CH_13BL2450_S1_R1.fastq.gz $dir/dummy_mother01_1.fastq.gz
cp $dir/CH_13BL2450_S1_R2.fastq.gz $dir/dummy_mother01_2.fastq.gz
cp $dir/GO_16BL0892_S4_R1.fastq.gz $dir/dummy_father01_1.fastq.gz
cp $dir/GO_16BL0892_S4_R2.fastq.gz $dir/dummy_father01_2.fastq.gz
cp $dir/JC_16BL0361_S6_R1.fastq.gz $dir/dummy_proband01_1.fastq.gz
cp $dir/JC_16BL0361_S6_R2.fastq.gz $dir/dummy_proband01_2.fastq.gz

cp $dir/LJ_17BL1030_S8_R1.fastq.gz $dir/dummy_mother02_1.fastq.gz
cp $dir/LJ_17BL1030_S8_R2.fastq.gz $dir/dummy_mother02_2.fastq.gz
cp $dir/MS_16BL0795_S5_R1.fastq.gz $dir/dummy_father02_1.fastq.gz
cp $dir/MS_16BL0795_S5_R2.fastq.gz $dir/dummy_father02_2.fastq.gz
cp $dir/NH_15BL1218_S7_R1.fastq.gz $dir/dummy_proband02_1.fastq.gz
cp $dir/NH_15BL1218_S7_R2.fastq.gz $dir/dummy_proband02_2.fastq.gz

cp $dir/PH_18BL0824_S2_R1.fastq.gz $dir/dummy_mother03_1.fastq.gz
cp $dir/PH_18BL0824_S2_R2.fastq.gz $dir/dummy_mother03_2.fastq.gz
cp $dir/SJ_18BL1034_S3_R1.fastq.gz $dir/dummy_father03_1.fastq.gz
cp $dir/SJ_18BL1034_S3_R1.fastq.gz $dir/dummy_father03_2.fastq.gz
cp $dir/CH_13BL2450_S1_R1.fastq.gz $dir/dummy_proband03_1.fastq.gz
cp $dir/CH_13BL2450_S1_R2.fastq.gz $dir/dummy_proband03_2.fastq.gz

cp $dir/GO_16BL0892_S4_R1.fastq.gz $dir/dummy_mother04_1.fastq.gz
cp $dir/GO_16BL0892_S4_R2.fastq.gz $dir/dummy_mother04_2.fastq.gz
cp $dir/JC_16BL0361_S6_R1.fastq.gz $dir/dummy_father04_1.fastq.gz
cp $dir/JC_16BL0361_S6_R2.fastq.gz $dir/dummy_father04_2.fastq.gz
cp $dir/LJ_17BL1030_S8_R1.fastq.gz $dir/dummy_proband04_1.fastq.gz
cp $dir/LJ_17BL1030_S8_R2.fastq.gz $dir/dummy_proband04_2.fastq.gz

cp $dir/MS_16BL0795_S5_R1.fastq.gz $dir/dummy_mother05_1.fastq.gz
cp $dir/MS_16BL0795_S5_R2.fastq.gz $dir/dummy_mother05_2.fastq.gz
cp $dir/NH_15BL1218_S7_R1.fastq.gz $dir/dummy_father05_1.fastq.gz
cp $dir/NH_15BL1218_S7_R2.fastq.gz $dir/dummy_father05_2.fastq.gz
cp $dir/PH_18BL0824_S2_R1.fastq.gz $dir/dummy_proband05_1.fastq.gz
cp $dir/PH_18BL0824_S2_R2.fastq.gz $dir/dummy_proband05_2.fastq.gz

cp $dir/SJ_18BL1034_S3_R1.fastq.gz $dir/dummy_mother06_1.fastq.gz
cp $dir/SJ_18BL1034_S3_R1.fastq.gz $dir/dummy_mother06_2.fastq.gz
cp $dir/CH_13BL2450_S1_R1.fastq.gz $dir/dummy_father06_1.fastq.gz
cp $dir/CH_13BL2450_S1_R2.fastq.gz $dir/dummy_father06_2.fastq.gz
cp $dir/GO_16BL0892_S4_R1.fastq.gz $dir/dummy_proband06_1.fastq.gz
cp $dir/GO_16BL0892_S4_R2.fastq.gz $dir/dummy_proband06_2.fastq.gz

cp $dir/JC_16BL0361_S6_R1.fastq.gz $dir/dummy_mother07_1.fastq.gz
cp $dir/JC_16BL0361_S6_R2.fastq.gz $dir/dummy_mother07_2.fastq.gz
cp $dir/LJ_17BL1030_S8_R1.fastq.gz $dir/dummy_father07_1.fastq.gz
cp $dir/LJ_17BL1030_S8_R2.fastq.gz $dir/dummy_father07_2.fastq.gz
cp $dir/MS_16BL0795_S5_R1.fastq.gz $dir/dummy_proband07_1.fastq.gz
cp $dir/MS_16BL0795_S5_R2.fastq.gz $dir/dummy_proband07_2.fastq.gz

cp $dir/NH_15BL1218_S7_R1.fastq.gz $dir/dummy_mother08_1.fastq.gz
cp $dir/NH_15BL1218_S7_R2.fastq.gz $dir/dummy_mother08_2.fastq.gz
cp $dir/PH_18BL0824_S2_R1.fastq.gz $dir/dummy_father08_1.fastq.gz
cp $dir/PH_18BL0824_S2_R2.fastq.gz $dir/dummy_father08_2.fastq.gz
cp $dir/SJ_18BL1034_S3_R1.fastq.gz $dir/dummy_proband08_1.fastq.gz
cp $dir/SJ_18BL1034_S3_R1.fastq.gz $dir/dummy_proband08_2.fastq.gz
