# VirusNet

[![Build and Test Conda Package](https://github.com/GenomeNet/VirusNet/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/GenomeNet/VirusNet/actions/workflows/python-package-conda.yml) [![Anaconda-Server Badge](https://anaconda.org/genomenet/virusnet/badges/version.svg)](https://anaconda.org/genomenet/virusnet) [![Anaconda-Server Badge](https://anaconda.org/genomenet/virusnet/badges/latest_release_relative_date.svg)](https://anaconda.org/genomenet/virusnet) [![Anaconda-Server Badge](https://anaconda.org/genomenet/virusnet/badges/downloads.svg)](https://anaconda.org/genomenet/virusnet)

## Install

Make sure you have conda installed. If you don't have it installed, you can download and install Miniconda from the official website: https://docs.conda.io/en/latest/miniconda.html

```bash
conda create -n virusnet python=3.11 -y
conda activate virusnet
conda install -c anaconda -c conda-forge -c genomenet virusnet -y
```

## Installation using Mamba (recommended)

```bash
conda install mamba -c conda-forge -y
mamba create -n virusnet python=3.11 -y
mamba activate virusnet
mamba install -c genomenet -c anaconda -c conda-forge genomenet::virusnet -y
```

## Usage

Download the models

```
virusnet download
```

Make predictions


```
virusnet predict --mode binary --input test/covid.fasta --output results_binary
```

This output will be written to the screen

```
[INFO] Checking input file 
[INFO] Number of FASTA entries in the file: 1 
[INFO] Loading binary model 
[INFO] Performing predictions 
[INFO] Using non-metagenomic mode 
[INFO] Summarized results written to: results_binary/binary_results_summarized.csv 
[INFO] Non-summarized results written to: results_binary/binary_results.csv 
[INFO] Number of contigs classified as viral: 1 
[INFO] Number of contigs classified as non-viral: 0 
[INFO] FASTA data written to: results_binary/viral_contigs.fasta 
[WARN] No contigs found 
```

```
virusnet predict --mode genus --input test/covid.fasta --output results_genus
```

Which will produce this output

```
[INFO] Checking input file 
[INFO] Number of FASTA entries in the file: 1 
[INFO] Loading genus model 
[INFO] Performing genus predictions 
Contig Summary:
Contig 1: Virus - Alphacoronavirus, Probability - 37%
```

## Development

### Building this conda package

```
mamba create -n build-env python=3.11 boa anaconda-client -y
micromamba activate build-env
boa build .
anaconda upload --user genomenet --channel genomenet path_to_virusnet.tar.bz2
```
