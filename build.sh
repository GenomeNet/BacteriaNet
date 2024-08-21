#!/bin/bash

set -e  # Exit immediately if a command exits with a non-zero status
set -x  # Print commands and their arguments as they are executed

mkdir -p $PREFIX/bin

# Install the Python script
echo "Installing the scripts..."
cp $SRC_DIR/bacterianet.py $PREFIX/bin/bacterianet
cp $SRC_DIR/utils.r $PREFIX/bin/utils.r
cp $SRC_DIR/setup_logger.r $PREFIX/bin/setup_logger.r
cp $SRC_DIR/models.json $PREFIX/bin/models.json
cp $SRC_DIR/predict_phenotypes.r $PREFIX/bin/predict_phenotypes.r

# Set the R library directory
mkdir -p "$PREFIX/lib/R/library"
chmod -R u+w "$PREFIX/lib/R/library"
R_LIBS_SITE="$PREFIX/lib/R/library"
"$BUILD_PREFIX/bin/R" -e "remotes::install_github('GenomeNet/deepg', lib='$PREFIX/lib/R/library')" 