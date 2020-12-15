# Update scout in production

Created: 2020/12/16 09:45:00
Last modified: 2020/12/16 10:15:14

- **Aim:** Document how to update [scout](http://www.clinicalgenomics.se/scout/) while scout is being actively used in production. This assumes you have scout already have installed in a conda environment (see [scout installation docs](./installation_scout.md)).
- **Prerequisite software:** [Conda 4.9.0](https://docs.conda.io/projects/conda/en/latest/index.html), [git 2.7.4](https://git-scm.com/), [pip 20.1.1](https://pypi.org/project/pip/)
- **OS:** Ubuntu 16.04 (Wintermute - research server)

## Update scout

- Check where scout is installed

```bash
which scout
```

```bash
/home/lkemp/anaconda3/envs/scout_env/bin/scout
```

- Check current scout version

```bash
scout --version
```

```bash
scout, version 4.23
```

- Make sure I'm in the conda env in which I installed scout, in my case `scout_env`

- Move to where I installed cloned the scout github repo, in my case, it is at `/store/lkemp/GA_clinical_genomics/scout/`

```bash
cd /store/lkemp/GA_clinical_genomics/scout/
```

- Check the status of the repo and pull from the remote repositiory

```bash
# Chcek status
git status

# Update my local github repo with the main repo
git pull origin master
```

- Check out the release to update, in my case, version 4.27

```bash
# Get all the tags from the main repo
git fetch --all --tags

# Checkout the version to use
git checkout tags/v4.27
```

- Re-install scout

```bash
pip install --requirement requirements.txt --editable .
```

- Check the current verison of scout

```bash
scout --version
```

```bash
scout, version 4.27
```

- Re-serve the scout database, in my case:

```bash
scout -db dhb-database serve
```
