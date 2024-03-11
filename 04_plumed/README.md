# PLUMED


## (updated on 2020)
See Pablo's page: [https://github.com/PabloPiaggi/running_gromacs/tree/master/04_plumed](https://github.com/PabloPiaggi/running_gromacs/tree/master/04_plumed)

## (updated on 2024) Instructions for Della

```
module purge
module load intel/2022.2.0
module load intel-mpi/intel/2021.7.0

#############################################################
### PLUMED
#############################################################
### plumed does not support all versions of gromacs

git clone https://github.com/plumed/plumed2.git
cd plumed2
git checkout v2.9
./configure --enable-modules=all CXX=mpiicpc
make -j10
source sourceme.sh
cd ..

#############################################################
### GROMACS
#############################################################
wget https://ftp.gromacs.org/gromacs/gromacs-2022.5.tar.gz
tar -xf gromacs-2022.5.tar.gz
cd gromacs-2022.5
plumed --no-mpi patch -p -e gromacs-2022.5 --runtime
mkdir build 
cd build

module load cudatoolkit/12.3
cmake3 .. -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=CUDA
make -j10

```
