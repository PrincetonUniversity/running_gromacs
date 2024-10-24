#!/bin/bash
#############################################################
# set the version
#############################################################
version=2024.3

#############################################################
# you probably don't need to change anything below this line
#############################################################

wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-${version}.tar.gz
tar zxvf gromacs-${version}.tar.gz
cd gromacs-${version}
mkdir build_stage1
cd build_stage1

#############################################################
# build gmx (stage 1)
#############################################################

module purge
module load intel-oneapi/2024.2
module load intel-mkl/2024.2

OPTFLAGS="-Ofast -xsapphirerapids -DNDEBUG"

cmake3 .. -DCMAKE_BUILD_TYPE=Release \
-DCMAKE_C_COMPILER=icx \
-DCMAKE_C_FLAGS_RELEASE="$OPTFLAGS" \
-DCMAKE_CXX_COMPILER=icpx \
-DCMAKE_CXX_FLAGS_RELEASE="$OPTFLAGS" \
-DGMX_BUILD_MDRUN_ONLY=OFF \
-DGMX_MPI=OFF \
-DGMX_OPENMP=ON \
-DGMX_SIMD=AVX_512 \
-DGMX_DOUBLE=OFF \
-DGMX_FFT_LIBRARY=mkl \
-DGMX_GPU=OFF \
-DCMAKE_INSTALL_PREFIX=$HOME/.local \
-DGMX_DEFAULT_SUFFIX=OFF -DGMX_BINARY_SUFFIX=_cpu -DGMX_LIBS_SUFFIX=_cpu \
-DGMX_COOL_QUOTES=OFF -DREGRESSIONTEST_DOWNLOAD=ON

make -j 32
make check
make install

#############################################################
# build mdrun_mpi (stage 2)
#############################################################

cd ..
mkdir build_stage2
cd build_stage2

module load intel-mpi/oneapi/2021.13

cmake3 .. -DCMAKE_BUILD_TYPE=Release \
-DCMAKE_C_COMPILER=icx \
-DCMAKE_C_FLAGS_RELEASE="$OPTFLAGS" \
-DCMAKE_CXX_COMPILER=icpx \
-DCMAKE_CXX_FLAGS_RELEASE="$OPTFLAGS" \
-DGMX_BUILD_MDRUN_ONLY=ON \
-DGMX_MPI=ON \
-DGMX_OPENMP=ON \
-DGMX_SIMD=AVX_512 \
-DGMX_DOUBLE=OFF \
-DGMX_FFT_LIBRARY=mkl \
-DGMX_GPU=OFF \
-DCMAKE_INSTALL_PREFIX=$HOME/.local \
-DGMX_DEFAULT_SUFFIX=OFF -DGMX_BINARY_SUFFIX=_cpu_mpi -DGMX_LIBS_SUFFIX=_cpu_mpi \
-DGMX_COOL_QUOTES=OFF -DREGRESSIONTEST_DOWNLOAD=ON

make -j 16
source ../build_stage1/scripts/GMXRC
tests/regressiontests-${version}/gmxtest.pl all
make install
