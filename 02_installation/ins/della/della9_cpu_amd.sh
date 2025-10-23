#!/bin/bash

#############################################################
# set the version
#############################################################
version=2025.3

#############################################################
# you probably don't need to make changes below
#############################################################
wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-${version}.tar.gz
tar zxvf gromacs-${version}.tar.gz
cd gromacs-${version}
mkdir build && cd build

module purge
module load cmake/3.30.8    # DO NOT INCLUDE IN SLURM SCRIPT
module load gcc-toolset/14
module load aocc/5.0.0
module load aocl/aocc/5.0.0
module load openmpi/aocc-5.0.0/4.1.6

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
    -DGMX_FFT_LIBRARY=fftw3 \
    -DFFTWF_INCLUDE_DIR=/opt/AMD/aocl/aocl-linux-aocc-5.0.0/aocc/include_LP64 \
    -DFFTWF_LIBRARY=/opt/AMD/aocl/aocl-linux-aocc-5.0.0/aocc/lib_LP64/libfftw3f.so \
    -DGMX_MPI=ON \
    -DGMX_GPU=OFF \
    -DGMX_OPENMP=ON \
    -DGMX_OPENMP_MAX_THREADS=192 \
    -DGMX_SIMD=AVX_512 \
    -DGMX_DOUBLE=OFF \
    -DCMAKE_INSTALL_PREFIX=$HOME/.local \
    -DREGRESSIONTEST_DOWNLOAD=ON \
    -DGMX_DEFAULT_SUFFIX=OFF -DGMX_BINARY_SUFFIX=_amd_cpu -DGMX_LIBS_SUFFIX=_amd_cpu \
    -DGMX_COOL_QUOTES=OFF

make -j 8
make check
make install
