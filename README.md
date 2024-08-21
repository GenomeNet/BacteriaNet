# BacteriaNet

BacteriaNet is a tool that allows for the generation of phenotypic predictions of a bacterial genome based on nucleotide data. This tool leverages deep learning models to predict various bacterial phenotypes directly from genomic sequences, providing a comprehensive and efficient method for phenotype prediction without the need for manual curation or preprocessing of input data.


[![Build and Test Conda Package](https://github.com/GenomeNet/BacteriaNet/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/GenomeNet/BacteriaNet/actions/workflows/python-package-conda.yml) [![Anaconda-Server Badge](https://anaconda.org/genomenet/BacteriaNet/badges/version.svg)](https://anaconda.org/genomenet/BacteriaNet) [![Anaconda-Server Badge](https://anaconda.org/genomenet/BacteriaNet/badges/latest_release_relative_date.svg)](https://anaconda.org/genomenet/BacteriaNet) [![Anaconda-Server Badge](https://anaconda.org/genomenet/BacteriaNet/badges/downloads.svg)](https://anaconda.org/genomenet/BacteriaNet)

## Description

BacteriaNet was developed to address the need for high-quality phenotypic data, which is often time-consuming to generate manually. By training deep learning models on data from the DSMZ BacDive database, BacteriaNet provides researchers with access to predicted phenotypes for missing entries in BacDive, enhancing the database's comprehensiveness and utility.

The models were trained using deepG, a deep learning library (https://github.com/GenomeNet/deepG/). Unlike other methods, BacteriaNet does not require manual curation of input data or preprocessing and is not limited to coding regions, making it a powerful tool for phenotype prediction.

### Key Features:
- Predicts 14 different bacterial phenotypes from genomic sequences.
- Utilizes deep learning models trained on a comprehensive dataset linking genomic and phenotypic data.
- Outperforms traditional phenotype prediction methods like Traitar for many phenotypes.
- Provides predictions for phenotypes such as spore forming, oxygen growth, gram staining, motility, and biosafety level.

### Model Performance:
The following table summarizes the performance of the BacteriaNet models for various phenotypes:

| Task                  | Type       | Target                                                                 | Score | Metric              |
|-----------------------|------------|------------------------------------------------------------------------|--------|---------------------|
| Spore forming         | Binary     | [TRUE, FALSE]                                                          | 91%   | Balanced accuracy   |
| Oxygen growth         | Binary     | [aerobe, anaerobe]                                                     | 87%   | Balanced accuracy   |
| Oxygen (obligate)     | Binary     | [aerobe, anaerobe]                                                     | 85%   | Balanced accuracy   |
| Oxygen (microaerophile)| Binary    | [aerobe, anaerobe]                                                     | 82%   | Balanced accuracy   |
| Oxygen (facultative)  | Binary     | [aerobe, anaerobe]                                                     | 72%   | Balanced accuracy   |
| Gram staining         | Multi-class| [gram stain negative, gram stain positive, gram stain variable]        | 85%   | Balanced accuracy   |
| Motility              | Binary     | [TRUE, FALSE]                                                          | 74%   | Balanced accuracy   |
| Biosafety level (BSL) | Multi-Class| [biosafety level 1, biosafety level 2, biosafety level 3]              | 40%   | Balanced accuracy   |
| Flagellum             | Multi-Class| [monotrichous, monotrichous_polar, polar, peritrichous, lophotrichous, gliding] | 29%    | Balanced accuracy   |
| Cell shape            | Multi-Class| [rod shaped, coccus shaped, vibrio shaped, filament shaped, sphere shaped, ovoid shaped, pleomorphic shaped, spiral shaped, curved shaped, oval shaped, other] | 20% | Balanced accuracy   |
| Cell size (length)    | Regression | Numeric                                                                | 0.15   | Correlation         |
| Cell size (width)     | Regression | Numeric                                                                | 0.23   | Correlation         |
| Pathogenicity (human) | Binary     | [TRUE, FALSE]                                                          | 0.68   | AUC                 |
| Pathogenicity (animal)| Binary     | [TRUE, FALSE]                                                          | 0.72   | AUC                 |
| Pathogenicity (plant) | Binary     | [TRUE, FALSE]                                                          | 0.68   | AUC                 |


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
