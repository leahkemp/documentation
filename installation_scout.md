# Install and run scout

**Aim:** Install [scout](http://www.clinicalgenomics.se/scout/) and run a demo instance
**Prerequisite software:** Conda 4.8.2
**OS:** Ubuntu 16.04 (Wintermute - research server)

## Install scout

Their [website](http://www.clinicalgenomics.se/scout/install/) describes the overall installation process. The easiest way I found to install scout was to create a conda environment in which to install scout and it's dependant software is installed. This saves you from *some* python dependency issues.

### Create a conda environment with python3.7

Create the conda environment and list all the environments

```bash
conda create --name scout_env python=3.7
conda env list
```

Activate the environment

```bash
conda activate scout_env
```

Now install software in this environment

### MongoDB

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

### Scout

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

Physically move to use Wintermute so that we can use it's GUI

From the command line, access the conda environment we created for scout

```bash
conda activate scout_env
```

Now run the demo

```bash
scout --demo serve
```

Connect to the [connection](http://localhost:5000/) through a browser to open the scout GUI. Use the user email clark.kent@mail.com to get access.

See the scout [admin guide](http://www.clinicalgenomics.se/scout/admin-guide/) for more info on how to set up scout as an administrator
