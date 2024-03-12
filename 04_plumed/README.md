# PLUMED


## (updated on 2020)
See Pablo's page: [https://github.com/PabloPiaggi/running_gromacs/tree/master/04_plumed](https://github.com/PabloPiaggi/running_gromacs/tree/master/04_plumed)

## (updated on 2024) Instructions for Della-GPU

```
module purge
module load gcc-toolset/12
module load openmpi/gcc/4.1.2
module load cudatoolkit/12.2

#############################################################
### PLUMED
#############################################################
### plumed does not support all versions of gromacs

git clone https://github.com/plumed/plumed2.git
cd plumed2
git checkout v2.9
./configure --enable-modules=all
make -j10
source sourceme.sh
cd ..

#############################################################
### GROMACS
#############################################################
wget https://ftp.gromacs.org/gromacs/gromacs-2022.5.tar.gz
tar -xf gromacs-2022.5.tar.gz
cd gromacs-2022.5
plumed patch -p -e gromacs-2022.5 --runtime
# plumed patch -p -d ../plumed2/patches/gromacs-2022.5.diff -e gramacs-2022.5 --runtime
mkdir build 
cd build

OPTFLAGS="-O3 -DNDEBUG"

cmake3 .. -DCMAKE_BUILD_TYPE=Release \
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

make -j10

```
