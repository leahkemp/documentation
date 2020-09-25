#!/bin/bash

# Set directory to directory where NA12878 exome was downloaded
dir="/NGS/scratch/KSCBIOM/HumanGenomics/exome_project/fastq"

##### Single runs #####

# For runs on 1 sample - single runs
for fastq_dir in exome/*/01_sample/single_run/*/fastq/; do
    cp $dir/dummy01_1.fastq.gz $fastq_dir
    cp $dir/dummy01_2.fastq.gz $fastq_dir
done

# For runs on 2 samples - single runs
for fastq_dir in exome/*/02_sample/single_run/*/fastq/; do
    cp $dir/dummy01_1.fastq.gz $fastq_dir
    cp $dir/dummy01_2.fastq.gz $fastq_dir
    cp $dir/dummy02_1.fastq.gz $fastq_dir
    cp $dir/dummy02_2.fastq.gz $fastq_dir
done

# For runs on 4 samples - single runs
for fastq_dir in exome/*/04_sample/single_run/*/fastq/; do
    cp $dir/dummy01_1.fastq.gz $fastq_dir
    cp $dir/dummy01_2.fastq.gz $fastq_dir
    cp $dir/dummy02_1.fastq.gz $fastq_dir
    cp $dir/dummy02_2.fastq.gz $fastq_dir
    cp $dir/dummy03_1.fastq.gz $fastq_dir
    cp $dir/dummy03_2.fastq.gz $fastq_dir
    cp $dir/dummy04_1.fastq.gz $fastq_dir
    cp $dir/dummy04_2.fastq.gz $fastq_dir
done

# For runs on 8 samples - single runs
for fastq_dir in exome/*/08_sample/single_run/*/fastq/; do
    cp $dir/dummy01_1.fastq.gz $fastq_dir
    cp $dir/dummy01_2.fastq.gz $fastq_dir
    cp $dir/dummy02_1.fastq.gz $fastq_dir
    cp $dir/dummy02_2.fastq.gz $fastq_dir
    cp $dir/dummy03_1.fastq.gz $fastq_dir
    cp $dir/dummy03_2.fastq.gz $fastq_dir
    cp $dir/dummy04_1.fastq.gz $fastq_dir
    cp $dir/dummy04_2.fastq.gz $fastq_dir
    cp $dir/dummy05_1.fastq.gz $fastq_dir
    cp $dir/dummy05_2.fastq.gz $fastq_dir
    cp $dir/dummy06_1.fastq.gz $fastq_dir
    cp $dir/dummy06_2.fastq.gz $fastq_dir
    cp $dir/dummy07_1.fastq.gz $fastq_dir
    cp $dir/dummy07_2.fastq.gz $fastq_dir
    cp $dir/dummy08_1.fastq.gz $fastq_dir
    cp $dir/dummy08_2.fastq.gz $fastq_dir
done

##### Cohort runs #####

# For runs on 1 family - cohort runs
for fastq_dir in exome/*/01_sample/cohort_run/*/fastq/; do
    cp $dir/dummy_proband01_1.fastq.gz $fastq_dir
    cp $dir/dummy_proband01_2.fastq.gz $fastq_dir
    cp $dir/dummy_father01_1.fastq.gz $fastq_dir
    cp $dir/dummy_father01_2.fastq.gz $fastq_dir
    cp $dir/dummy_mother01_1.fastq.gz $fastq_dir
    cp $dir/dummy_mother01_2.fastq.gz $fastq_dir
    cp $dir/dummy_proband01_pedigree.ped $fastq_dir../pedigrees/
done

# For runs on 2 families - cohort runs
for fastq_dir in exome/*/02_sample/cohort_run/*/fastq/; do
    cp $dir/dummy_proband01_1.fastq.gz $fastq_dir
    cp $dir/dummy_proband01_2.fastq.gz $fastq_dir
    cp $dir/dummy_father01_1.fastq.gz $fastq_dir
    cp $dir/dummy_father01_2.fastq.gz $fastq_dir
    cp $dir/dummy_mother01_1.fastq.gz $fastq_dir
    cp $dir/dummy_mother01_2.fastq.gz $fastq_dir
    cp $dir/dummy_proband01_pedigree.ped $fastq_dir../pedigrees/
    cp $dir/dummy_proband02_1.fastq.gz $fastq_dir
    cp $dir/dummy_proband02_2.fastq.gz $fastq_dir
    cp $dir/dummy_father02_1.fastq.gz $fastq_dir
    cp $dir/dummy_father02_2.fastq.gz $fastq_dir
    cp $dir/dummy_mother02_1.fastq.gz $fastq_dir
    cp $dir/dummy_mother02_2.fastq.gz $fastq_dir
    cp $dir/dummy_proband02_pedigree.ped $fastq_dir../pedigrees/
done

# For runs on 4 families - cohort runs
for fastq_dir in exome/*/04_sample/cohort_run/*/fastq/; do
    cp $dir/dummy_proband01_1.fastq.gz $fastq_dir
    cp $dir/dummy_proband01_2.fastq.gz $fastq_dir
    cp $dir/dummy_father01_1.fastq.gz $fastq_dir
    cp $dir/dummy_father01_2.fastq.gz $fastq_dir
    cp $dir/dummy_mother01_1.fastq.gz $fastq_dir
    cp $dir/dummy_mother01_2.fastq.gz $fastq_dir
    cp $dir/dummy_proband01_pedigree.ped $fastq_dir../pedigrees/
    cp $dir/dummy_proband02_1.fastq.gz $fastq_dir
    cp $dir/dummy_proband02_2.fastq.gz $fastq_dir
    cp $dir/dummy_father02_1.fastq.gz $fastq_dir
    cp $dir/dummy_father02_2.fastq.gz $fastq_dir
    cp $dir/dummy_mother02_1.fastq.gz $fastq_dir
    cp $dir/dummy_mother02_2.fastq.gz $fastq_dir
    cp $dir/dummy_proband02_pedigree.ped $fastq_dir../pedigrees/
    cp $dir/dummy_proband03_1.fastq.gz $fastq_dir
    cp $dir/dummy_proband03_2.fastq.gz $fastq_dir
    cp $dir/dummy_father03_1.fastq.gz $fastq_dir
    cp $dir/dummy_father03_2.fastq.gz $fastq_dir
    cp $dir/dummy_mother03_1.fastq.gz $fastq_dir
    cp $dir/dummy_mother03_2.fastq.gz $fastq_dir
    cp $dir/dummy_proband03_pedigree.ped $fastq_dir../pedigrees/
    cp $dir/dummy_proband04_1.fastq.gz $fastq_dir
    cp $dir/dummy_proband04_2.fastq.gz $fastq_dir
    cp $dir/dummy_father04_1.fastq.gz $fastq_dir
    cp $dir/dummy_father04_2.fastq.gz $fastq_dir
    cp $dir/dummy_mother04_1.fastq.gz $fastq_dir
    cp $dir/dummy_mother04_2.fastq.gz $fastq_dir
    cp $dir/dummy_proband04_pedigree.ped $fastq_dir../pedigrees/
done

# For runs on 8 families - cohort runs
for fastq_dir in exome/*/08_sample/cohort_run/*/fastq/; do
    cp $dir/dummy_proband01_1.fastq.gz $fastq_dir
    cp $dir/dummy_proband01_2.fastq.gz $fastq_dir
    cp $dir/dummy_father01_1.fastq.gz $fastq_dir
    cp $dir/dummy_father01_2.fastq.gz $fastq_dir
    cp $dir/dummy_mother01_1.fastq.gz $fastq_dir
    cp $dir/dummy_mother01_2.fastq.gz $fastq_dir
    cp $dir/dummy_proband01_pedigree.ped $fastq_dir../pedigrees/
    cp $dir/dummy_proband02_1.fastq.gz $fastq_dir
    cp $dir/dummy_proband02_2.fastq.gz $fastq_dir
    cp $dir/dummy_father02_1.fastq.gz $fastq_dir
    cp $dir/dummy_father02_2.fastq.gz $fastq_dir
    cp $dir/dummy_mother02_1.fastq.gz $fastq_dir
    cp $dir/dummy_mother02_2.fastq.gz $fastq_dir
    cp $dir/dummy_proband02_pedigree.ped $fastq_dir../pedigrees/
    cp $dir/dummy_proband03_1.fastq.gz $fastq_dir
    cp $dir/dummy_proband03_2.fastq.gz $fastq_dir
    cp $dir/dummy_father03_1.fastq.gz $fastq_dir
    cp $dir/dummy_father03_2.fastq.gz $fastq_dir
    cp $dir/dummy_mother03_1.fastq.gz $fastq_dir
    cp $dir/dummy_mother03_2.fastq.gz $fastq_dir
    cp $dir/dummy_proband03_pedigree.ped $fastq_dir../pedigrees/
    cp $dir/dummy_proband04_1.fastq.gz $fastq_dir
    cp $dir/dummy_proband04_2.fastq.gz $fastq_dir
    cp $dir/dummy_father04_1.fastq.gz $fastq_dir
    cp $dir/dummy_father04_2.fastq.gz $fastq_dir
    cp $dir/dummy_mother04_1.fastq.gz $fastq_dir
    cp $dir/dummy_mother04_2.fastq.gz $fastq_dir
    cp $dir/dummy_proband04_pedigree.ped $fastq_dir../pedigrees/
    cp $dir/dummy_proband05_1.fastq.gz $fastq_dir
    cp $dir/dummy_proband05_2.fastq.gz $fastq_dir
    cp $dir/dummy_father05_1.fastq.gz $fastq_dir
    cp $dir/dummy_father05_2.fastq.gz $fastq_dir
    cp $dir/dummy_mother05_1.fastq.gz $fastq_dir
    cp $dir/dummy_mother05_2.fastq.gz $fastq_dir
    cp $dir/dummy_proband05_pedigree.ped $fastq_dir../pedigrees/
    cp $dir/dummy_proband06_1.fastq.gz $fastq_dir
    cp $dir/dummy_proband06_2.fastq.gz $fastq_dir
    cp $dir/dummy_father06_1.fastq.gz $fastq_dir
    cp $dir/dummy_father06_2.fastq.gz $fastq_dir
    cp $dir/dummy_mother06_1.fastq.gz $fastq_dir
    cp $dir/dummy_mother06_2.fastq.gz $fastq_dir
    cp $dir/dummy_proband06_pedigree.ped $fastq_dir../pedigrees/
    cp $dir/dummy_proband07_1.fastq.gz $fastq_dir
    cp $dir/dummy_proband07_2.fastq.gz $fastq_dir
    cp $dir/dummy_father07_1.fastq.gz $fastq_dir
    cp $dir/dummy_father07_2.fastq.gz $fastq_dir
    cp $dir/dummy_mother07_1.fastq.gz $fastq_dir
    cp $dir/dummy_mother07_2.fastq.gz $fastq_dir
    cp $dir/dummy_proband07_pedigree.ped $fastq_dir../pedigrees/
    cp $dir/dummy_proband08_1.fastq.gz $fastq_dir
    cp $dir/dummy_proband08_2.fastq.gz $fastq_dir
    cp $dir/dummy_father08_1.fastq.gz $fastq_dir
    cp $dir/dummy_father08_2.fastq.gz $fastq_dir
    cp $dir/dummy_mother08_1.fastq.gz $fastq_dir
    cp $dir/dummy_mother08_2.fastq.gz $fastq_dir
    cp $dir/dummy_proband08_pedigree.ped $fastq_dir../pedigrees/
done
