# Install singularity and create a singularity container

Created: 2020/03/11 11:25:43
Last modified: 2020/03/11 11:27:21

- **Aim:** Install [singularity](https://singularity.lbl.gov/install-linux) and create a singularity container
- **Prerequisite software:** [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html)
- **OS:** Ubuntu 16.04 (Wintermute - research server)

## Table of contents

- [Install singularity and create a singularity container](#install-singularity-and-create-a-singularity-container)
  - [Table of contents](#table-of-contents)
  - [Install singularity](#install-singularity)
  - [Pull an image from docker](#pull-an-image-from-docker)
  - [Build a container from scratch](#build-a-container-from-scratch)

## Install singularity

Requires root access

```bash
VERSION=2.5.2
wget https://github.com/singularityware/singularity/releases/download/$VERSION/singularity-$VERSION.tar.gz
tar xvf singularity-$VERSION.tar.gz
cd singularity-$VERSION
./configure --prefix=/usr/local
```

I then received an error stating that I needed to install the [libarchive-devel](https://packages.debian.org/sid/libarchive-dev) package:

```bash
configure: error: Unable to find the libarchive headers, need package libarchive-devel (libarchive-dev on Debian/Ubuntu)
```

Install libarchive-devel package

```bash
sudo apt-get update -y && sudo apt-get install -y libarchive-dev
```

Then continue with singularity install

```bash
./configure --prefix=/usr/local
make
sudo make install
```

squashfs-tools is not required for configuring singularity, however it is required for full functionality. To install squashfs-tools:

```bash
sudo apt-get update && sudo apt-get install squashfs-tools
```

## Pull an image from docker

To pull a docker container ([bcftools](https://biocontainers.pro/#/tools/bcftools)) from the [BioContainers](https://biocontainers.pro/#/) community. This container can be used for genomic variant calling and manipulation of vcf/bcf files.

```bash
sudo singularity build bcftools.simg docker://quay.io/biocontainers/bcftools:1.7--0
```

## Build a container from scratch

Make a recipe file, for example create this file and save as 'Singularity'

```bash
Bootstrap: docker
From: ubuntu:16.04

%post
    apt-get -y update
    apt-get -y install fortune cowsay lolcat

%environment
    export LC_ALL=C
    export PATH=/usr/games:$PATH

%runscript
    fortune | cowsay | lolcat
```

Use this recipe file to build/create a container. This will create the file 'scout.simg'.

```bash
sudo singularity build scout.simg Singularity
```
