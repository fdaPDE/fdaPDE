#!/bin/bash

IMAGE=archlinux

docker pull $IMAGE 

docker run --name=build-archlinux-gcc -v $(pwd)/../../:/root/fdaPDE --rm -ti $IMAGE bin/bash -c "

pacman --noconfirm -Sy gcc
pacman --noconfirm -Sy icu
pacman --noconfirm -Sy r-base
pacman --noconfirm -Sy make
pacman --noconfirm -Sy libxml2
pacman --noconfirm -Sy freetype2
pacman --noconfirm -Sy harfbuzz
pacman --noconfirm -Sy fribidi
pacman --noconfirm -Sy libtiff

export RGL_USE_NULL=TRUE
export DISPLAY=99.0

cd 
mkdir .R
cd .R
touch Makevars
echo CPPFLAGS+=-I'/usr/include/freetype2/' > Makevars
echo CPPFLAGS+=-I'/usr/include/harfbuzz/' >> Makevars
echo CPPFLAGS+=-I'/usr/include/fribidi/' >> Makevars
cd ..

Rscript fdaPDE/tests/install_dependencies.R
Rscript fdaPDE/tests/building_check/check.R
exit
"

