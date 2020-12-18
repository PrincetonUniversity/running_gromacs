#!/bin/bash
version=2019.6
wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-${version}.tar.gz
tar -zxvf gromacs-${version}.tar.gz
cd gromacs-${version}
mkdir build_stage1
cd build_stage1

#############################################################
# build gmx (for single-node jobs)
#############################################################

module purge
module load intel/19.0/64/19.0.1.144
module load cudatoolkit/10.2
module load rh/devtoolset/7

OPTFLAGS="-Ofast -xCORE-AVX512 -mtune=skylake-avx512 -DNDEBUG"

cmake3 .. -DCMAKE_BUILD_TYPE=Release \
-DCMAKE_C_COMPILER=icc -DCMAKE_C_FLAGS_RELEASE="$OPTFLAGS" \
-DCMAKE_CXX_COMPILER=icpc -DCMAKE_CXX_FLAGS_RELEASE="$OPTFLAGS" \
-DGMX_BUILD_MDRUN_ONLY=OFF -DGMX_MPI=OFF -DGMX_OPENMP=ON \
-DGMX_SIMD=AVX_512 -DGMX_DOUBLE=OFF \
-DGMX_FFT_LIBRARY=mkl \
-DGMX_GPU=ON -DGMX_CUDA_TARGET_SM=70 \
-DCMAKE_INSTALL_PREFIX=$HOME/.local \
-DGMX_COOL_QUOTES=OFF -DREGRESSIONTEST_DOWNLOAD=ON

make -j 10
make check
make install
