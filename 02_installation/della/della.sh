#!/bin/bash
#############################################################
# set the version
#############################################################
version=2019.6

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
module load intel/19.0/64/19.0.5.281 rh/devtoolset/6

OPTFLAGS="-Ofast -xCORE-AVX2 -axCORE-AVX512 -DNDEBUG"

cmake3 .. -DCMAKE_BUILD_TYPE=Release \
-DCMAKE_C_COMPILER=icc -DCMAKE_C_FLAGS_RELEASE="$OPTFLAGS" \
-DCMAKE_CXX_COMPILER=icpc -DCMAKE_CXX_FLAGS_RELEASE="$OPTFLAGS" \
-DGMX_BUILD_MDRUN_ONLY=OFF -DGMX_MPI=OFF -DGMX_OPENMP=ON \
-DGMX_DOUBLE=OFF \
-DGMX_FFT_LIBRARY=mkl \
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

module load intel-mpi/intel/2018.3/64

cmake3 .. -DCMAKE_BUILD_TYPE=Release \
-DCMAKE_C_COMPILER=icc -DCMAKE_C_FLAGS_RELEASE="$OPTFLAGS" \
-DCMAKE_CXX_COMPILER=icpc -DCMAKE_CXX_FLAGS_RELEASE="$OPTFLAGS" \
-DGMX_BUILD_MDRUN_ONLY=ON -DGMX_MPI=ON -DGMX_OPENMP=ON \
-DGMX_DOUBLE=OFF \
-DGMX_FFT_LIBRARY=mkl \
-DGMX_GPU=OFF \
-DCMAKE_INSTALL_PREFIX=$HOME/.local \
-DGMX_COOL_QUOTES=OFF -DREGRESSIONTEST_DOWNLOAD=ON

make -j 10
source ../build_stage1/scripts/GMXRC
tests/regressiontests-${version}/gmxtest.pl all
make install
