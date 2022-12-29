#!/bin/bash
#############################################################
# set the version
#############################################################
version=2022.4

#############################################################
# you probably don't need to change anything below this line
#############################################################

wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-${version}.tar.gz
tar zxvf gromacs-${version}.tar.gz
cd gromacs-${version}
mkdir build && cd build

#############################################################
# build gmx (stage 1)
#############################################################

module purge
module load anaconda3/2022.5
module load rh/devtoolset/9
module load cudatoolkit/11.3

# gromacs will add -march=core-avx2 to the next line
OPTFLAGS=""

cmake3 .. -DCMAKE_BUILD_TYPE=Release \
-DCMAKE_C_COMPILER=gcc -DCMAKE_C_FLAGS_RELEASE="$OPTFLAGS" \
-DCMAKE_CXX_COMPILER=g++ -DCMAKE_CXX_FLAGS_RELEASE="$OPTFLAGS" \
-DGMX_MPI=OFF -DGMX_OPENMP=ON \
-DGMX_SIMD=AVX2_256 -DGMX_DOUBLE=OFF \
-DGMX_FFT_LIBRARY=fftw3 \
-DGMX_GPU=CUDA -DGMX_CUDA_TARGET_SM=60 \
-DCMAKE_INSTALL_PREFIX=$HOME/.local \
-DGMX_COOL_QUOTES=OFF -DREGRESSIONTEST_DOWNLOAD=ON

make -j 8
make check
make install
