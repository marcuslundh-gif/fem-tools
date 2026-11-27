#!/bin/bash
# Bygger CalculiX 2.21 med Intel MKL Pardiso (2–3 gånger snabbare än apt-versionens SPOOLES) som /usr/local/bin/ccx_pardiso.
# Ubuntu 24.04 (moln eller WSL), körs som root eller med sudo. MKL tas från pip-paketen mkl och mkl-include
# (installeras i en venv om pip saknar systemrättigheter). Tar cirka 5 minuter på 4 kärnor.
set -e
W=${1:-/tmp/ccx_build}; mkdir -p $W && cd $W
SUDO=$([ "$(id -u)" = 0 ] || echo sudo)
$SUDO sed -i 's/^Types: deb$/Types: deb deb-src/' /etc/apt/sources.list.d/ubuntu.sources
$SUDO apt-get update -qq
$SUDO apt-get install -y -qq libspooles-dev libarpack2-dev gfortran liblapack-dev libblas-dev dpkg-dev python3-venv
apt-get source -qq calculix-ccx
python3 -m venv $W/venv && $W/venv/bin/pip install -q mkl mkl-include
INC=$(dirname $(find $W/venv -name mkl_service.h | head -1))
MKL=$(find $W/venv -name "libmkl_rt.so*" | sort | head -1)
cd calculix-ccx-2.21/ccx_2.21/src
make -j"$(nproc)" CFLAGS="-O2 -I/usr/include/spooles -I$INC -DARCH=\"Linux\" -DSPOOLES -DARPACK -DMATRIXSTORAGE -DUSE_MT=1 -DPARDISO -fopenmp" \
     FFLAGS="-O2 -fallow-argument-mismatch -fopenmp" LIBS="-lspooles -lpthread -larpack -llapack -lblas -lm -lc $MKL -Wl,-rpath,$(dirname $MKL)" > build.log 2>&1
$SUDO cp ccx_2.21 /usr/local/bin/ccx_pardiso && echo "ccx_pardiso installerad (MKL: $MKL)"
