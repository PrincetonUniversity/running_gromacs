# Tiger (CPU)

### GCC

```
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
mkdir build && cd build

module purge
module load openmpi/gcc/4.1.6
module load intel-mkl/2024.2

OPTFLAGS="-O3 -DNDEBUG"

cmake3 .. -DCMAKE_BUILD_TYPE=Release \
-DCMAKE_C_COMPILER=mpicc -DCMAKE_C_FLAGS_RELEASE="$OPTFLAGS" \
-DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_CXX_FLAGS_RELEASE="$OPTFLAGS" \
-DGMX_MPI=ON -DGMX_OPENMP=ON \
-DGMX_DOUBLE=OFF \
-DGMX_FFT_LIBRARY=mkl \
-DGMX_GPU=OFF \
-DCMAKE_INSTALL_PREFIX=$HOME/.local \
-DGMX_COOL_QUOTES=OFF -DREGRESSIONTEST_DOWNLOAD=ON

make
make check
make install
```

Below is an example script:

```
#!/bin/bash
#SBATCH --job-name=gmx           # create a short name for your job
#SBATCH --nodes=2                # node count
#SBATCH --ntasks=32               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=400G                 # memory per node (4G per cpu-core is default)
#SBATCH --time=00:10:00          # total run time limit (HH:MM:SS)

module purge
module load openmpi/gcc/4.1.6
module load intel-mkl/2024.2

export GMX_MAXBACKUP=-1

BCH=./rnase_cubic
gmx_mpi grompp -f $BCH/pme_verlet.mdp -c $BCH/conf.gro -p $BCH/topol.top -o bench.tpr
srun gmx_mpi mdrun -ntomp $SLURM_CPUS_PER_TASK -s bench.tpr
```


### Intel

The directions below use the Intel compilers and Intel MPI library:

```
$ ssh <NetID>@tiger3.princeton.edu
$ cd </path/to/software/directory>  # e.g., cd ~/software
$ mkdir cpu_version && cd cpu_version
$ wget https://raw.githubusercontent.com/PrincetonUniversity/running_gromacs/master/02_installation/ins/tigerCpu/tigerCpu.sh
# make modifications to tigercpu.sh if needed
$ bash tigerCpu.sh | tee build.log
```

For single-node jobs:

```bash
#!/bin/bash
#SBATCH --job-name=rnase         # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=8               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=4G                 # memory per node (4G per cpu-core is default)
#SBATCH --time=00:10:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all          # send email when job begins, ends and fails
#SBATCH --mail-user=<YourNetID>@princeton.edu

module purge
module load intel-oneapi/2024.2
module load intel-mkl/2024.2

gmx_cpu grompp -f pme_verlet.mdp -c conf.gro -p topol.top -o bench.tpr
gmx_cpu mdrun -ntmpi $SLURM_NTASKS -ntomp $SLURM_CPUS_PER_TASK -s bench.tpr
```

For multi-node runs:

```bash
#!/bin/bash
#SBATCH --job-name=rnase         # create a short name for your job
#SBATCH --nodes=3                # node count
#SBATCH --ntasks-per-node=16     # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=4G                 # memory per node (4G per cpu-core is default)
#SBATCH --time=00:10:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all          # send email when job begins, ends and fails
#SBATCH --mail-user=<YourNetID>@princeton.edu

module purge
module load intel-oneapi/2024.2
module load intel-mkl/2024.2
module load intel-mpi/oneapi/2021.13

gmx_cpu grompp -f $BCH/pme_verlet.mdp -c $BCH/conf.gro -p $BCH/topol.top -o bench.tpr
srun mdrun_cpu_mpi -ntomp $SLURM_CPUS_PER_TASK -s bench.tpr
```
