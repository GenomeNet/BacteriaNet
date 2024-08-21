# BacteriaNet

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
BacteriaNet predict --mode binary --input test/ecoli.fasta --output results_binary
```

This output will be written to the screen

```
TODO
```

```
BacteriaNet predict --mode genus --input test/ecoli.fasta --output results_genus
```

Which will produce this output

```
TODO
```

## Development

### Building this conda package

```
mamba create -n build-env python=3.11 boa anaconda-client -y
micromamba activate build-env
boa build .
anaconda upload --user genomenet --channel genomenet path_to_BacteriaNet.tar.bz2
```
