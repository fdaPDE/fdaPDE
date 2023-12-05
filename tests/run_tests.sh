#!/bin/bash

IMAGE=rhub/debian-gcc-release

docker pull $IMAGE 

docker run --name=tests-debian-gcc -v $(pwd)/../:/root/fdaPDE/ --rm -ti $IMAGE bin/bash -c '

export RGL_USE_NULL=TRUE
export DISPLAY=0

apt-get install -y libgl1-mesa-dev libglu1-mesa-dev
apt-get install -y libxml2-dev
apt-get install -y libssl-dev

cd /root/fdaPDE/tests/

Rscript install_dependencies.R
Rscript build_data.R
Rscript run_tests.R

exit
'
