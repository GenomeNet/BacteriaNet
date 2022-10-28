
<img src="https://github.com/GenomeNet/BacteriaNet/blob/main/logo.png?raw=true" width="328" height="49">

<hr />

*BacteriaNet* is a deep neuronal network developed with the *deepG* package that predics viral taxonomy on the genus level. A typical viral sequence can be analyzed in ~3 seconds. 

## Dependencies
1. Docker

## Installation

```
docker pull genomenet/bacterianet:alpha
```

## Usage
1. Input and output are defined as stdin/stdout as done in the following command

```
cat example.fasta | docker run --rm -i genomenet/bacterianet:alpha morphology.r > prediction.csv
```

2. it will report the predicted morphologies (probabilities) to the screen
3. it will write all predictions to stdout (here `prediction.csv`). 
