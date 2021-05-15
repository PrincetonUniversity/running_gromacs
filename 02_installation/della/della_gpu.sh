#!/bin/bash

#############################################################
# set the version
#############################################################
version=2021.2

wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-${version}.tar.gz
tar zxvf gromacs-${version}.tar.gz
cd gromacs-${version}
mkdir build && cd build

module purge
module load cmake/3.18.2
module load openmpi/gcc/4.1.0

OPTFLAGS="-O3 -DNDEBUG"

cmake .. -DCMAKE_BUILD_TYPE=Release \
-DCMAKE_C_COMPILER=mpicc \
-DCMAKE_C_FLAGS_RELEASE="$OPTFLAGS" \
-DCMAKE_CXX_COMPILER=mpic++ \
-DCMAKE_CXX_FLAGS_RELEASE="$OPTFLAGS" \
-DGMX_BUILD_OWN_FFTW=ON \
-DGMX_MPI=ON -DGMX_OPENMP=ON \
-DGMX_OPENMP_MAX_THREADS=128 \
-DGMX_SIMD=AVX2_256 -DGMX_DOUBLE=OFF \
-DGMX_GPU=CUDA -DGMX_CUDA_TARGET_SM=80 \
-DCMAKE_INSTALL_PREFIX=$HOME/.local \
-DREGRESSIONTEST_DOWNLOAD=ON \
-DGMX_DEFAULT_SUFFIX=OFF -DGMX_BINARY_SUFFIX=_gpu -DGMX_LIBS_SUFFIX=_gpu \
-DGMX_COOL_QUOTES=OFF

make -j 32
make check
make install
