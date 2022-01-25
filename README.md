# GAARD_work
Work and outputs linked to the GAARD project

## Conda environment
The folder conda_yml contains yaml files to build what will hopefully be a largely shared conda environment. Doesn't mean everything has to be run in this shared environment, but if there is some work that would benefit from being run on different machines, this will facilitate that (I hope). 

The environments contain both python and R packages. It uses python 3.8 and R 4.1. 

Yaml files have name gaard_\#.yml or gaard_basic_\#.yml. Look for the largest value of \# for the latest file. Try the gaard_\#.yml version first. If you run into incompatibility issues, try gaard_basic_\#.yml.

example:

conda env create -f gaard_01.yml

conda activate gaard

