#!/bin/bash

IMAGE=archlinux

docker pull $IMAGE 

docker run --name=build-archlinux-gcc -v $(pwd)/../../:/root/fdaPDE --rm -ti $IMAGE bin/bash -c "

pacman --noconfirm -Sy gcc
pacman --noconfirm -Sy icu
pacman --noconfirm -Sy make
pacman --noconfirm -Sy libxml2
pacman --noconfirm -Sy freetype2
pacman --noconfirm -Sy harfbuzz
pacman --noconfirm -Sy fribidi
pacman --noconfirm -Sy libtiff
pacman --noconfirm -Sy r-base

export RGL_USE_NULL=TRUE
export DISPLAY=99.0

cd 
mkdir .R; touch .R/Makevars
echo CPPFLAGS+=-I'/usr/include/freetype2/' > .R/Makevars
echo CPPFLAGS+=-I'/usr/include/harfbuzz/' >> .R/Makevars
echo CPPFLAGS+=-I'/usr/include/fribidi/' >> .R/Makevars

Rscript fdaPDE/tests/install_dependencies.R
Rscript fdaPDE/tests/building_check/check.R
exit
"

