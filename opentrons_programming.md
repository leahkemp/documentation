# Programming opentrons robot

Created: 2020-04-20 11:35:06
Last modified: 2020/04/28 09:24:00

- **Aim:** Program a portion of the new [opentrons pipetting robot](https://opentrons.com/) using the [opentrons api (OT-2 API V2)](https://docs.opentrons.com/v2/index.html). More specifically, automating the dilution step (concentration normalisation) after pcr amplification and a picogreen assay in the arctic protocol
- **Prerequisite software:** [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html), [pip](https://pypi.org/project/pip/)
- **OS:** Ubuntu 16.04 (Wintermute - research server)

## Table of contents

- [Programming opentrons robot](#programming-opentrons-robot)
  - [Table of contents](#table-of-contents)
  - [Setup](#setup)
  - [Testing script](#testing-script)

## Setup

Create and activate a conda environment to get the minimum python version required for a [non-jupyter opentrons installation](https://docs.opentrons.com/v2/writing.html#non-jupyter-installation)

```bash
conda create -n opentrons_env python=3.7.6
conda activate opentrons_env
```

Install opentrons within conda environment

```bash
pip install opentrons
```

## Testing script

Simulate the protocol [from the command line](https://docs.opentrons.com/v2/writing.html#from-the-command-line)

```bash
opentrons_simulate dilutions_opentrons.py
```
