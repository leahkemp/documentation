# Programming opentrons robot

Created: 2020-04-20 11:35:06
Last modified: 2020/04/21 14:36:53

- **Aim:** Program a portion of the new [opentrons pipetting robot](https://opentrons.com/) using the [opentrons api (OT-2 API V2)](https://docs.opentrons.com/v2/index.html). More specifically, automating the dilution step (concentration normalisation) after pcr amplification and a picogreen assay in the arctic protocol
- **Prerequisite software:** [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html), [pip](https://pypi.org/project/pip/)
- **OS:** Ubuntu 16.04 (Wintermute - research server)

## Table of contents

- [Programming opentrons robot](#programming-opentrons-robot)
  - [Table of contents](#table-of-contents)
  - [Setup](#setup)
  - [Writing python script](#writing-python-script)
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

## Writing python script

Install packages to read in xlsx file

```bash
conda install pandas
conda install xlrd
```

Adapted the [normalization protocol provided by opentrons](https://protocols.opentrons.com/protocol/normalization)

Notes:

- Want to use the lowest supported OT-2 Python Protocol API as possible since this keeps your protocol portable between API versions

To considered for later modification

- Any speed modifications

Inputs:

- What labware we have (Opentrons labware vs. other software, pipettes, tip racks, reagents)
- Hardware modules

## Testing script

Simulate the protocol [from the command line](https://docs.opentrons.com/v2/writing.html#from-the-command-line)

```bash
opentrons_simulate dilutions_opentrons.py
```
