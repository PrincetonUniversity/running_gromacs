#!/bin/bash
version_fftw=3.3.9
version_gmx=2021

#############################################################
# build a fast version of FFTW
#############################################################

wget ftp://ftp.fftw.org/pub/fftw/fftw-${version_fftw}.tar.gz
tar zxvf fftw-${version_fftw}.tar.gz
cd fftw-${version_fftw}

module purge

./configure CC=gcc CFLAGS="-Ofast -mcpu=power9 -mtune=power9 -DNDEBUG" --prefix=$HOME/.local \
--enable-shared --enable-single --enable-vsx --disable-fortran

make -j 10
make install
cd ..

#############################################################
# build gmx (for single node jobs)
#############################################################

wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-${version_gmx}.tar.gz
tar zxvf gromacs-${version_gmx}.tar.gz
cd gromacs-${version_gmx}
mkdir build_stage1
cd build_stage1

module purge
module load cmake/3.18.2
module load cudatoolkit/11.2

OPTFLAGS="-Ofast -mcpu=power9 -mtune=power9 -mvsx -DNDEBUG"

cmake3 .. -DCMAKE_BUILD_TYPE=Release \
-DCMAKE_C_COMPILER=gcc -DCMAKE_C_FLAGS_RELEASE="$OPTFLAGS" \
-DCMAKE_CXX_COMPILER=g++ -DCMAKE_CXX_FLAGS_RELEASE="$OPTFLAGS" \
-DGMX_BUILD_MDRUN_ONLY=OFF -DGMX_MPI=OFF -DGMX_OPENMP=ON \
-DGMX_SIMD=IBM_VSX -DGMX_DOUBLE=OFF \
-DGMX_FFT_LIBRARY=fftw3 \
-DFFTWF_INCLUDE_DIR=$HOME/.local/include \
-DFFTWF_LIBRARY=$HOME/.local/lib/libfftw3f.so \
-DGMX_GPU=CUDA -DGMX_CUDA_TARGET_SM=70 \
-DGMX_OPENMP_MAX_THREADS=128 \
-DCMAKE_INSTALL_PREFIX=$HOME/.local \
-DGMX_COOL_QUOTES=OFF -DREGRESSIONTEST_DOWNLOAD=ON

make -j 10
make check
make install

# two tests fail: #10 due to pinning and #42 which fails because
# it tries to run it using six GPU tasks but we have only
# four GPUs

#############################################################
# build mdrun_mpi (for multi-node jobs)
#############################################################

cd ..
mkdir build_stage2
cd build_stage2

module load openmpi/gcc/4.0.4/64

cmake3 .. -DCMAKE_BUILD_TYPE=Release \
-DCMAKE_C_COMPILER=gcc -DCMAKE_C_FLAGS_RELEASE="$OPTFLAGS" \
-DCMAKE_CXX_COMPILER=g++ -DCMAKE_CXX_FLAGS_RELEASE="$OPTFLAGS" \
-DGMX_BUILD_MDRUN_ONLY=ON -DGMX_MPI=ON -DGMX_OPENMP=ON \
-DGMX_SIMD=IBM_VSX -DGMX_DOUBLE=OFF \
-DGMX_FFT_LIBRARY=fftw3 \
-DFFTWF_INCLUDE_DIR=$HOME/.local/include \
-DFFTWF_LIBRARY=$HOME/.local/lib/libfftw3f.so \
-DGMX_GPU=CUDA -DGMX_CUDA_TARGET_SM=70 \
-DGMX_OPENMP_MAX_THREADS=128 \
-DCMAKE_INSTALL_PREFIX=$HOME/.local \
-DGMX_COOL_QUOTES=OFF -DREGRESSIONTEST_DOWNLOAD=ON

make -j 10
source ../build_stage1/scripts/GMXRC
tests/regressiontests-${version_gmx}/gmxtest.pl all
make install
