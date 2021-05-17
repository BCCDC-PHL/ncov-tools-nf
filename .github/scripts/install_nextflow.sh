#!/bin/bash
set -eo pipefail

echo Install Nextflow .. >> artifacts/test_artifact.log

wget -qO- https://get.nextflow.io | bash

mkdir -p /opt/nextflow/bin

mv nextflow /opt/nextflow/bin

echo "export PATH=/opt/nextflow/bin:$PATH" >> ~/.bashrc

NXF_VER=20.10.0 /opt/nextflow/bin/nextflow -quiet run hello
