# EpiPI
[![Operating systems](https://github.com/SABS-R3-Epidemiology/EpiPI/actions/workflows/os_versions.yml/badge.svg)](https://github.com/SABS-R3-Epidemiology/EpiPI/actions/workflows/os_versions.yml)
[![Python package](https://github.com/SABS-R3-Epidemiology/EpiPI/actions/workflows/python_versions.yml/badge.svg)](https://github.com/SABS-R3-Epidemiology/EpiPI/actions/workflows/python_versions.yml)
[![Style tests (flake8)](https://github.com/SABS-R3-Epidemiology/EpiPI/actions/workflows/style.yml/badge.svg)](https://github.com/SABS-R3-Epidemiology/EpiPI/actions/workflows/style.yml)
[![Documentation Status](https://readthedocs.org/projects/epipi/badge/?version=latest)](https://epipi.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/SABS-R3-Epidemiology/EpiPI/graph/badge.svg?token=X36987JKP0)](https://codecov.io/gh/SABS-R3-Epidemiology/EpiPI)

## General Information
This package consists of scripts to run inference of incidecnes on each day based on the prevalence data generated from the sampling. This repository has a Jupyter Notebook file summarising work can be done with this package, which is in the file `work_summary.ipynb`.

## Installation
Currently, EpiPI is only tested on MacOS and Linux. Windows are not supported at the moment.

Installation of EpiPI depends on the CPU of the machine that you are using.

### For x86 CPU:
EpiPI is not yet available on [PyPI](https://pypi.org/), but the module can be pip installed locally. The directory should first be downloaded to your local machine, and can then be installed using the command:

```console
pip install -e .
```

We also recommend you to install the [EpiABM](https://github.com/SABS-R3-Epidemiology/epiabm) model and [EpiOS](https://github.com/SABS-R3-Epidemiology/EpiOS) to generate the data of infection simulation and do sampling.

### For other CPU:
For Apple Silicon users, you should firstly use `conda` to install `pystan`. After this, you should manually uninstall the following packages: `httpstan` and `pysimdjson`. Next, you should go to the github repository of `httpstan` to download the `4.10.1` version and install it according to the documentation of `httpstan` for `non-x86 CPU`. Finally, use `pip` to install `pysimdjson==5.0.2`. Now you will have the correct `pystan` working on your Apple Silicon machine.

## Documentation

 Documentations can be accessed via the above `docs` badge.

## Example Usage
Model description, examples and work summary can be found in the file `work_summary.ipynb`, with codes to use EpiPI to predict $R_t$. There are also definitions, description and reference for data used in this model.