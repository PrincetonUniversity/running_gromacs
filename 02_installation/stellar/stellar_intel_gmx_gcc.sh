#!/bin/bash
#############################################################
# set the version
#############################################################
version=2021

#############################################################
# you probably don't need to change anything below this line
#############################################################

wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-${version}.tar.gz
tar -zxvf gromacs-${version}.tar.gz
cd gromacs-${version}
mkdir build_stage1
cd build_stage1

#############################################################
# build gmx (stage 1)
#############################################################

module purge
module load fftw/gcc/3.3.9

OPTFLAGS="-Ofast -march=native -mtune=native -DNDEBUG"

cmake .. -DCMAKE_BUILD_TYPE=Release \
-DCMAKE_C_COMPILER=gcc -DCMAKE_C_FLAGS_RELEASE="$OPTFLAGS" \
-DCMAKE_CXX_COMPILER=g++ -DCMAKE_CXX_FLAGS_RELEASE="$OPTFLAGS" \
-DGMX_BUILD_MDRUN_ONLY=OFF -DGMX_MPI=OFF -DGMX_OPENMP=ON \
-DGMX_SIMD=AVX_512 -DGMX_DOUBLE=OFF \
-DGMX_FFT_LIBRARY=fftw3 \
-DFFTWF_INCLUDE_DIR=/usr/local/fftw/gcc/3.3.9/include \
-DFFTWF_LIBRARY=/usr/local/fftw/gcc/3.3.9/lib64/libfftw3f.so \
-DGMX_GPU=OFF \
-DCMAKE_INSTALL_PREFIX=$HOME/.local \
-DGMX_COOL_QUOTES=OFF -DREGRESSIONTEST_DOWNLOAD=ON

make -j 10
make check
make install

#############################################################
# build mdrun_mpi (stage 2)
#############################################################

cd ..
mkdir build_stage2
cd build_stage2

module load openmpi/gcc/4.1.0

cmake .. -DCMAKE_BUILD_TYPE=Release \
-DCMAKE_C_COMPILER=gcc -DCMAKE_C_FLAGS_RELEASE="$OPTFLAGS" \
-DCMAKE_CXX_COMPILER=g++ -DCMAKE_CXX_FLAGS_RELEASE="$OPTFLAGS" \
-DGMX_BUILD_MDRUN_ONLY=ON -DGMX_MPI=ON -DGMX_OPENMP=ON \
-DGMX_SIMD=AVX_512 -DGMX_DOUBLE=OFF \
-DGMX_FFT_LIBRARY=fftw3 \
-DFFTWF_INCLUDE_DIR=/usr/local/fftw/gcc/3.3.9/include \
-DFFTWF_LIBRARY=/usr/local/fftw/gcc/3.3.9/lib64/libfftw3f.so \
-DGMX_GPU=OFF \
-DCMAKE_INSTALL_PREFIX=$HOME/.local \
-DGMX_COOL_QUOTES=OFF -DREGRESSIONTEST_DOWNLOAD=ON

make -j 10
source ../build_stage1/scripts/GMXRC
tests/regressiontests-${version}/gmxtest.pl all
make install
