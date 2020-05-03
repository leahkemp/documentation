# Testing HUPAN (HUman Pan-genome ANalysis)

Created: 2020-04-30 11:38:42
Last modified: 2020/04/30 13:57:02

- **Aim:** Test the utility of [HUPAN (HUman Pan-genome ANalysis)](https://github.com/SJTU-CGM/HUPAN) for use in creating human pan genomes/graph genomes
- **Prerequisite software:** [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html), [R 3.4.3](https://www.r-project.org/)
- **OS:** Ubuntu 16.04 (Wintermute - research server)

## Table of contents

- [Testing HUPAN (HUman Pan-genome ANalysis)](#testing-hupan-human-pan-genome-analysis)
  - [Table of contents](#table-of-contents)
  - [Download HUPAN](#download-hupan)

## Download HUPAN

Create a conda environment

```bash
conda create -n hupan_env python=3.7.6
conda activate hupan_env
```

Clone HUPAN

```bash
git clone git@github.com:SJTU-CGM/HUPAN.git
```

Install R dependencies

```bash
cd HUPAN
Rscript installRPac
```

Compile tools

```bash
make
```

Add bin/ to PATH and add lib/ to LD_LIBRARY_PATH. To do this, add the following text to ~/.bash_profile:

```bash
export PATH=$PATH:/home/lkemp/HUPAN/bin:
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/lkemp/HUPAN/lib/:
export PERL5LIB=$PERL5LIB:/home/lkemp/HUPAN/lib/
```

```bash
source ~/.profile
```

Didn't install properly

