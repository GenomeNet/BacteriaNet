# BacteriaNet

BacteriaNet is a tool that allows to generate phenotypic prediction of a bacteria genome based on nucleotide data.

[![Build and Test Conda Package](https://github.com/GenomeNet/BacteriaNet/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/GenomeNet/BacteriaNet/actions/workflows/python-package-conda.yml) [![Anaconda-Server Badge](https://anaconda.org/genomenet/BacteriaNet/badges/version.svg)](https://anaconda.org/genomenet/BacteriaNet) [![Anaconda-Server Badge](https://anaconda.org/genomenet/BacteriaNet/badges/latest_release_relative_date.svg)](https://anaconda.org/genomenet/BacteriaNet) [![Anaconda-Server Badge](https://anaconda.org/genomenet/BacteriaNet/badges/downloads.svg)](https://anaconda.org/genomenet/BacteriaNet)

## Install

Make sure you have conda installed. If you don't have it installed, you can download and install Miniconda from the official website: https://docs.conda.io/en/latest/miniconda.html

```bash
conda create -n BacteriaNet python=3.11 -y
conda activate BacteriaNet
conda install -c anaconda -c conda-forge -c genomenet BacteriaNet -y
```

## Installation using Mamba (recommended)

```bash
conda install mamba -c conda-forge -y
mamba create -n BacteriaNet python=3.11 -y
mamba activate BacteriaNet
mamba install -c genomenet -c anaconda -c conda-forge genomenet::BacteriaNet -y
```

## Usage

Download the models

```
BacteriaNet download
```

Make predictions


```
BacteriaNet predict --input test/ecoli.fasta --output example_results
```

This output will be written to the screen

```
rgStart_len: -0.0233813473411525
rgEnd_len: -0.0892123775556684
rgStart_wid: 0.331785323719184
rgEnd_wid: 0.114853305121263
cellshape: rod.shaped
flagellum: monotrichous_polar
gram: negative
is_motile: TRUE
biosafety: 3
pathogenicity_human: FALSE
pathogenicity_animal: FALSE
pathogenicity_plant: FALSE
oxygen_growth: aerobe
oxygen_facultative: facultative.aerobe
oxygen_obligate: obligate.aerobe
microaerophile: FALSE
ability_spore: FALSE
```

## Development

### Building this conda package

```
mamba create -n build-env python=3.11 boa anaconda-client -y
micromamba activate build-env
boa build . --prefix-lengt 30
anaconda upload --user genomenet --channel genomenet path_to_BacteriaNet.tar.bz2
```
