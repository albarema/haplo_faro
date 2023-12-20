# Workflow guidelines
## Snakemake
Snakemake is a **text-based workflow system using Python interpreter**. Please, follow the instructions to download [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). 

Snakemake is very well-documented. A few useful links:
- Snakemake tutorial: https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html
- Snakemake slides: https://slides.com/johanneskoester/snakemake-tutorial 
- Snakemake advanced: https://f1000research.com/articles/10-33/v1
  
### How to run Snakemake

Pre-running-pipeline step:
Values key set-up using a configuration file in YAML format (do not modify the key name unless you change it accordingly in all Snakemake files. If you change the config filename, you will need to add this flag --configfile when running Snakemake).

The workflow is executed as follows:

```bash
# simple run 
snakemake --snakefile haplonet.smk -j5
# if the filename is Snakefile and is located in the working directory you don't have to provide the name 
snakemake -j5
# dry-run. The -np flags specify that job execution will be simulated (-n) and the individual rule commands printed (-p)
snakemake -np

```
Specify the path using the flag ```--snakefile``` if the file is not located in your working directory. I usually have all my snakefiles in a separate directory named  ```rules```. This applied as well to the config.yaml file ```--configfile```. The location can be indicated in the Snakefiles. 


## Conda environments

I recommend creating your conda environment for each project (which will require different software and potentially different versions). More information on [conda environments](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html). 

Some useful commands to get started:

```bash
# list all the current channels
conda config --show channels
# add new channels to the front of the priority list (you only need to do this once)
conda config --prepend channels bioconda
conda config --prepend channels conda-forge
# list the available environments
conda env list
# If an environmental.yaml file is given (versions specified), all dependencies and packages can be installed in a new env as follow: 
conda env create --name XXX --file environment.yaml
# create an empty new environment
conda env create --name XXX
# activate the environment to use it
conda activate XXX
# find latest version of snakemake (or whatever other package you are interested)
conda search snakemake
# install it - one package at a time
conda install snakemake=5.15.0
# export the environment file (you will get the updated dependencies - if you had install new ones after the creation of the env)
- A. only built from user
conda env export --no-builds --from-history > environment.yaml
- B. all
conda env export > environment.yaml
# deactivate an environemnt
conda deactivate
# delete environment
conda env remove --name XXX
```
