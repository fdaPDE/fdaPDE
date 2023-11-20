#!/bin/bash

IMAGE=rhub/fedora-gcc-devel

docker pull $IMAGE 

docker run --name=build-fedora-gcc -v $(pwd)/../../:/root/fdaPDE --rm -ti $IMAGE bin/bash -c '

export RGL_USE_NULL=TRUE
export DISPLAY=99.0
cd root/
sudo ln -s /opt/R-devel/bin/R /usr/local/bin/R
sudo ln -s /opt/R-devel/bin/Rscript /usr/local/bin/Rscript
sudo dnf install -y mesa-libGL-devel mesa-libGLU-devel libpng-devel

Rscript fdaPDE/tests/install_dependencies.R 
Rscript fdaPDE/tests/building_check/check.R 

exit
'
