#!/bin/bash

set -eo pipefail

export PATH=/opt/miniconda3/bin:$PATH

mamba env create --prefix /home/runner/.conda/envs/ncov-qc-86edc7b3dd9d65d41e4044034891dd43 --file /home/runner/work/ncov-tools-nf/ncov-tools-nf/environments/ncov-tools-1.9.yml
