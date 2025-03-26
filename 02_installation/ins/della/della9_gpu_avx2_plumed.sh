#!/bin/bash

module purge
module load cmake/3.30.8       # DO NOT INCLUDE IN SLURM SCRIPT
module load gcc-toolset/14     # DO NOT INCLUDE IN SLURM SCRIPT
module load openmpi/gcc/4.1.6
module load cudatoolkit/12.8

#############################################################
# PLUMED
#############################################################
git clone https://github.com/plumed/plumed2.git
cd plumed2
git checkout v2.9
./configure --enable-modules=all
make -j 10
source sourceme.sh
plumed -h
cd ..

#############################################################
# set the version of GMX (plumed does not support all versions of GMX)
#############################################################
version=2024.3

#############################################################
# you probably don't need to make changes below
#############################################################
wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-${version}.tar.gz
tar zxvf gromacs-${version}.tar.gz
cd gromacs-${version}
plumed patch -p -e gromacs-${version} --runtime
mkdir build && cd build

OPTFLAGS="-O3 -DNDEBUG"

which cmake
which nvcc
which gcc
which mpicc

cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_C_FLAGS_RELEASE="$OPTFLAGS" \
    -DCMAKE_CXX_COMPILER=mpic++ \
    -DCMAKE_CXX_FLAGS_RELEASE="$OPTFLAGS" \
    -DGMX_BUILD_OWN_FFTW=ON \
    -DGMX_MPI=ON \
    -DGMX_OPENMP=ON \
    -DGMX_OPENMP_MAX_THREADS=192 \
    -DGMX_SIMD=AVX2_256 -DGMX_DOUBLE=OFF \
    -DGMX_GPU=CUDA -DGMX_CUDA_TARGET_SM=80 \
    -DCMAKE_INSTALL_PREFIX=$HOME/.local \
    -DREGRESSIONTEST_DOWNLOAD=ON \
    -DGMX_DEFAULT_SUFFIX=OFF -DGMX_BINARY_SUFFIX=_d9_gpu_plumed -DGMX_LIBS_SUFFIX=_d9_gpu_plumed \
    -DGMX_COOL_QUOTES=OFF

make -j 8
make check
make install
