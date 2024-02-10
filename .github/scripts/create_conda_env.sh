#!/bin/bash

set -eo pipefail

export PATH=/opt/miniconda3/bin:$PATH

mamba env create --prefix /home/runner/.conda/envs/ncov-qc-40e1fefdb1b312663435a35c6d6cb9da --file /home/runner/work/ncov-tools-nf/ncov-tools-nf/environments/ncov-tools-1.5.yml
