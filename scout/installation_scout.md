# Install and run scout

Created: 2020/03/11 11:25:43
Last modified: 2020/05/07 18:13:11

- **Aim:** Install [scout](http://www.clinicalgenomics.se/scout/) and run a demo instance
- **Prerequisite software:** [Conda 4.8.2](https://docs.conda.io/projects/conda/en/latest/index.html)
- **OS:** Ubuntu 16.04 (Wintermute - research server)

## Table of contents

- [Install and run scout](#install-and-run-scout)
  - [Table of contents](#table-of-contents)
  - [Install scout](#install-scout)
    - [Create a conda environment](#create-a-conda-environment)
    - [Install/run MongoDB](#installrun-mongodb)
    - [Install scout](#install-scout-1)
  - [Setup/run an example scout database](#setuprun-an-example-scout-database)

## Install scout

Their [website](http://www.clinicalgenomics.se/scout/install/) describes the overall installation process. The easiest way I found to install scout was to create a conda environment in which to install scout and it's dependant software is installed. This saves you from some python dependency issues.

### Create a conda environment

Create and activate a conda environment

```bash
conda create --name scout_env python=3.7
conda activate scout_env
```

Now install software in this environment

### Install/run MongoDB

Install [MongoDB](https://docs.mongodb.com/manual/tutorial/install-mongodb-on-ubuntu/) (I used db version v4.2.3)

```bash
wget -qO - https://www.mongodb.org/static/pgp/server-4.2.asc | sudo apt-key add -
echo "deb [ arch=amd64,arm64 ] https://repo.mongodb.org/apt/ubuntu xenial/mongodb-org/4.2 multiverse" | sudo tee /etc/apt/sources.list.d/mongodb-org-4.2.list
sudo apt-get update
sudo apt-get install -y mongodb-org
```

Start a mongod process and check it is running

```bash
sudo systemctl start mongod
sudo systemctl status mongod
```

If you need to stop a mongod process

```bash
sudo systemctl stop mongod
```

### Install scout

Clone the scout repository and install it's dependencies (I used Scout version 4.12.3, this code will download the latest stable release)

```bash
git clone https://github.com/Clinical-Genomics/scout
cd scout
pip install --requirement requirements.txt --editable .
```

## Setup/run an example scout database

Initialize a working instance with all genes, diseases etc. This will setup the database with a curated collection of gene definitions with links to OMIM along with HPO phenotype terms.

```bash
scout setup database
```

for more info

```bash
scout --help
```

Setup a fully working Scout demo. This will setup an instance of scout with a database called `scout-demo`

```bash
scout setup demo
```

Now run the demo

```bash
scout --demo serve
```

Either directly use Wintermute or create an SSH tunnel so that we can use it's GUI to explore scout.

Connect to the [connection](http://localhost:5000/) through a browser to open the scout GUI. Use the user email clark.kent@mail.com to get access.

See the scout [admin guide](http://www.clinicalgenomics.se/scout/admin-guide/) for more info on how to set up scout as an administrator.
