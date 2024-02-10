#!/bin/bash

set -eo pipefail

echo "Install Miniconda .." >> artifacts/test_artifact.log

export PATH=/opt/miniconda3/bin:$PATH

wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/miniconda3 && \
    rm ~/miniconda.sh && \

conda init bash

conda update -n base -c defaults conda

conda install mamba -n base -c conda-forge

echo 'alias conda=mamba' >> ~/.bash_aliases
