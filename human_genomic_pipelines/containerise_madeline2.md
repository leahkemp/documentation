# Containerising madeline2

Created: 2020-07-16 11:03:03
Last modified: 2020/08/18 10:27:59

- **Aim:** Create a containerized version of [madeline2](https://madeline.med.umich.edu/madeline/) ([singularity container](https://singularity.lbl.gov/))
- **Prerequisite software:** [singularity 2.5.2](https://singularity.lbl.gov/)
- **OS:** Ubuntu 16.04 (Wintermute - research server)

This document describes how I created a containerised version of [madeline2](https://madeline.med.umich.edu/madeline/) that could be used within our genomic pipelines (snakemake) to create pedigree information when processing cohort data. We want to use either a conda install of Madeline2 or a containerised version of Madeline2 in our genomic pipelines to ensure reproducibility and deployability on the production cluster. A conda environment for Madeline2 is not currently available on Anaconda cloud (https://anaconda.org/) nor is a containerised version of it available. Therefore the aim here is to containerise Madeline2 and deploy it to [singularity-hub](https://singularity-hub.org/)

## Table of contents

- [Containerising madeline2](#containerising-madeline2)
  - [Table of contents](#table-of-contents)
  - [Create a recipe file](#create-a-recipe-file)
  - [Build the recipe file](#build-the-recipe-file)
  - [Test](#test)
  - [Back up recipe file on github](#back-up-recipe-file-on-github)
  - [Host container on singularity-hub](#host-container-on-singularity-hub)
  - [Test container](#test-container)

## Create a recipe file

Create file

```bash
cd /home/lkemp/madeline2_container/
touch madeline2
```

Write recipe file. Contents of recipe file:

```txt
Bootstrap: docker
From: ubuntu:18.04

%help
This is a containerised version of the Madeline 2.0 Pedigree Drawing Engine (https://madeline.med.umich.edu/madeline/)

%labels
    MAINTAINER Leah Kemp (leah.kemp@esr.cri.nz)
    VERSION 0.1

%post
    apt-get -y update
    apt-get -y install git cmake build-essential libcurl4-openssl-dev libssl-dev libxml2-dev gettext
    git clone https://github.com/piratical/Madeline_2.0_PDE.git
    cd Madeline_2.0_PDE/
    ./configure --with-include-gettext
    make
    make install
```

*Notes:*

- Use `-y` flag to prompt the install to continue when it prompt for a y/n response

## Build the recipe file

Use the recipe file to build/create the container. This will create a container titled madeline2.simg

```bash
sudo singularity build madeline2.simg madeline2
```

## Test

```bash
cd /store/lkemp/exome_project/manual_pipeline_run/original_dbNSFP_database/

singularity exec \
-B /store/lkemp/exome_project/manual_pipeline_run/: \
/home/lkemp/madeline2_container/madeline2.simg \
madeline2 -L "IndividualId Affected" trio_4205.ped.data --outputext xml --color
```

It works!

## Back up recipe file on github

## Host container on [singularity-hub](https://singularity-hub.org/)

See [these instructions](https://singularityhub.github.io/singularityhub-docs/docs/getting-started/). Create/log into your user account and create container build/deployment from the website.

## Test container

The madeline2 container is now hosted [here](https://singularity-hub.org/collections/4639). Pull and interact with the container:

```bash
BIND_PATHS="/store/lkemp/manual_pipeline_run/"
CONTAINER="/home/lkemp/leahkemp-madeline2_container-master-latest.simg"

singularity exec \
-B ${BIND_PATHS} \
${CONTAINER} \
madeline2 \
-L "IndividualId Affected" DHB4177.ped.data \
--outputext xml \
--color
```

It works!
