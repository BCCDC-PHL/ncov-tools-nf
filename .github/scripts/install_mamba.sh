#!/bin/bash
set -eo pipefail

echo "Install Micromamba.." >> artifacts/test_artifact.log

mkdir micromamba && cd micromamba
wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
cd ..  & mv micromamba /opt

/opt/micromamba/bin/micromamba shell init -s bash -p ~/micromamba

echo 'alias conda=micromamba' >> ~/.bash_aliases
