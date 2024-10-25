#!/bin/bash
#############################################################
# set the version
#############################################################
version=2024.3

#############################################################
# you probably don't need to change anything below this line
#############################################################

#wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-${version}.tar.gz
tar zxvf gromacs-${version}.tar.gz
cd gromacs-${version}
mkdir build && cd build

module purge
module load gcc-toolset/10
module load intel-oneapi/2024.2
module load intel-mpi/oneapi/2021.13
module load intel-mkl/2024.2

OPTFLAGS="-O3 -DNDEBUG"

cmake3 .. -DCMAKE_BUILD_TYPE=Release \
-DCMAKE_C_COMPILER=mpiicx \
-DCMAKE_C_FLAGS_RELEASE="$OPTFLAGS" \
-DCMAKE_CXX_COMPILER=mpiicpx \
-DCMAKE_CXX_FLAGS_RELEASE="$OPTFLAGS" \
-DGMX_MPI=ON \
-DGMX_OPENMP=ON \
-DGMX_DOUBLE=OFF \
-DGMX_FFT_LIBRARY=mkl \
-DGMX_GPU=OFF \
-DCMAKE_INSTALL_PREFIX=$HOME/.local \
-DGMX_COOL_QUOTES=OFF -DREGRESSIONTEST_DOWNLOAD=ON

make
make check
make install
