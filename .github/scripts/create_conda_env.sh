#!/bin/bash

set -eo pipefail

export PATH=/opt/miniconda3/bin:$PATH

mamba env create --prefix ${HOME}/.conda/envs/ncov-qc-86edc7b3dd9d65d41e4044034891dd43 --file environments/ncov-tools-1.9.yml

nextflow pull BCCDC-PHL/ncov2019-artic-nf -r v1.3.4

mamba env create --prefix ${HOME}/.conda/envs/artic-ncov2019-illumina-b69e5bb8f90dfe5ddd3f76cc24acf6d3 --file ${HOME}/.nextflow/assets/BCCDC-PHL/ncov2019-artic-nf/environments/illumina/environment.yml
