#!/bin/bash
#############################################################
# set the version
#############################################################
version=2021.3

#############################################################
# you probably don't need to change anything below this line
#############################################################
wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-${version}.tar.gz
tar -zxvf gromacs-${version}.tar.gz
cd gromacs-${version}
mkdir build && cd build

module purge
module load openmpi/gcc/4.1.0

OPTFLAGS="-O3 -DNDEBUG"

cmake3 .. -DCMAKE_BUILD_TYPE=Release \
-DCMAKE_C_COMPILER=mpicc \
-DCMAKE_C_FLAGS_RELEASE="$OPTFLAGS" \
-DCMAKE_CXX_COMPILER=mpic++ \
-DCMAKE_CXX_FLAGS_RELEASE="$OPTFLAGS" \
-DGMX_BUILD_OWN_FFTW=ON \
-DGMX_MPI=ON -DGMX_OPENMP=ON \
-DGMX_SIMD=AVX_512 -DGMX_DOUBLE=OFF \
-DGMX_GPU=OFF \
-DCMAKE_INSTALL_PREFIX=$HOME/.local \
-DGMX_DEFAULT_SUFFIX=OFF -DGMX_BINARY_SUFFIX=_cpu -DGMX_LIBS_SUFFIX=_cpu \
-DGMX_COOL_QUOTES=OFF -DREGRESSIONTEST_DOWNLOAD=ON

make -j 8
make check
make install
