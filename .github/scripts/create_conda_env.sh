#!/bin/bash

set -eo pipefail

export PATH=/opt/miniconda3/bin:$PATH

mamba env create --prefix /home/runner/.conda/envs/ncov-qc-724a4524bfccdf191bf8550202dac824 --file /home/runner/work/ncov-tools-nf/ncov-tools-nf/environments/ncov-tools-1.5.yml
