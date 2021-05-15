# Installation

The installation procedure depends on the cluster.

## Della-GPU

The GPU nodes on Della have AMD processors with AVX2 as the highest instruction set. The system version of GCC is 8.3.1. Below is a build procedure:

```
$ ssh <YourNetID>@della-gpu.princeton.edu
$ cd software  # or another directory
$ wget https://raw.githubusercontent.com/PrincetonUniversity/running_gromacs/main/02_installation/della/della_gpu.sh
$ bash della_gpu.sh | tee gmx.log
```

## TigerGPU

The code can be built by following these commands:

```
$ ssh <YourNetID>@tigergpu.princeton.edu
$ cd </path/to/your/software/directory>  # e.g., cd ~/software
$ wget https://raw.githubusercontent.com/PrincetonUniversity/running_gromacs/master/02_installation/tigerGpu/tigerGpu.sh
# make modifications to tigerGpu.sh if needed (e.g., choose a different version)
$ bash tigerGpu.sh | tee build.log
```

The `tee` command is used to have the output go to the terminal as well as a file.

For single-node jobs:

```bash
#!/bin/bash
#SBATCH --job-name=gmx           # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=8        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=4G                 # memory per node (4G per cpu-core is default)
#SBATCH --time=01:00:00          # total run time limit (HH:MM:SS)
#SBATCH --gres=gpu:1             # number of gpus per node
#SBATCH --mail-type=all          # send email when job begins, ends and fails
#SBATCH --mail-user=<YourNetID>@princeton.edu

module purge
module load intel/19.0/64/19.0.5.281
module load cudatoolkit/10.2

gmx grompp -f pme_verlet.mdp -c conf.gro -p topol.top -o bench.tpr
gmx mdrun -ntmpi $SLURM_NTASKS -ntomp $SLURM_CPUS_PER_TASK -s bench.tpr
```

Note that Gromacs only uses serial FFT routines so you will not benefit from a multithreaded or MPI FFT library.

For multi-node MPI jobs:

```bash
#!/bin/bash
#SBATCH --job-name=gmx           # create a short name for your job
#SBATCH --nodes=4                # node count
#SBATCH --ntasks-per-node=12     # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=1G         # memory per cpu-core (4G per cpu-core is default)
#SBATCH --time=01:00:00          # total run time limit (HH:MM:SS)
#SBATCH --gres=gpu:1             # number of gpus per node
#SBATCH --mail-type=all          # send email when job begins, ends and fails
#SBATCH --mail-user=<YourNetID>@princeton.edu

module purge
module load intel/19.0/64/19.0.5.281
module load cudatoolkit/10.2
module load intel-mpi/intel/2019.5/64

gmx grompp -f pme_verlet.mdp -c conf.gro -p topol.top -o bench.tpr
srun mdrun_mpi -ntomp $SLURM_CPUS_PER_TASK -s bench.tpr
```

Our clusters are configured to use `srun`. Please do not use `mpirun`.

## Traverse

Traverse is composed of 46 IBM POWER9 nodes with 4 NVIDIA V100 GPUs per node. Each node has two sockets each with a 16-core CPU. Each of the 16 cores has 4 hardware threads. Hence each node has 128 logical cpu-cores or hardware threads.

![smt-core](http://3s81si1s5ygj3mzby34dq6qf-wpengine.netdna-ssl.com/wp-content/uploads/2016/08/ibm-hot-chips-power9-smt4-core.jpg)

```
$ ssh <YourNetID>@traverse.princeton.edu
$ cd </path/to/your/software/directory>  # e.g., cd ~/software
$ wget https://raw.githubusercontent.com/PrincetonUniversity/running_gromacs/master/02_installation/traverse/traverse.sh
# make modifications to traverse.sh if needed (e.g., choose a different version)
$ bash traverse.sh | tee build.log
```

Traverse will build successfully but it fails a test:

```
[ RUN      ] HardwareTopologyTest.NumaCacheSelfconsistency
/home/jdh4/sw/gromacs-2019.4/src/gromacs/hardware/tests/hardwaretopology.cpp:179: Failure
Expected: (n.memory) > (0), actual: 0 vs 0
/home/jdh4/sw/gromacs-2019.4/src/gromacs/hardware/tests/hardwaretopology.cpp:179: Failure
Expected: (n.memory) > (0), actual: 0 vs 0
[  FAILED  ] HardwareTopologyTest.NumaCacheSelfconsistency (126 ms)
[----------] 4 tests from HardwareTopologyTest (599 ms total)

[----------] Global test environment tear-down
[==========] 5 tests from 2 test cases ran. (603 ms total)
[  PASSED  ] 4 tests.
[  FAILED  ] 1 test, listed below:
[  FAILED  ] HardwareTopologyTest.NumaCacheSelfconsistency

 1 FAILED TEST
```

In the Gromacs output one will see:

```
NOTE: Affinity setting for 10/16 threads failed.
NOTE: Thread affinity was not set.
```

For single-node jobs:

```bash
#!/bin/bash
#SBATCH --job-name=gmx           # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=16       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --threads-per-core=1     # setting to 1 turns off SMT (max value is 4)
#SBATCH --mem=4G                 # memory per node (4G per cpu-core is default)
#SBATCH --time=01:00:00          # total run time limit (HH:MM:SS)
#SBATCH --gres=gpu:1             # number of gpus per node
#SBATCH --gpu-bind=closest       # use the closest gpu
#SBATCH --mail-type=all          # send email when job begins, ends and fails
#SBATCH --mail-user=<YourNetID>@princeton.edu

module purge
module load cudatoolkit/10.2

gmx grompp -f pme_verlet.mdp -c conf.gro -p topol.top -o bench.tpr
gmx mdrun -pin on -ntmpi $SLURM_NTASKS -ntomp $SLURM_CPUS_PER_TASK -s bench.tpr
```

The above Slurm leads to poor performance but the below script seems to work:

```
#!/bin/bash
#SBATCH --job-name=rnase         # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=64               # total number of tasks across all nodes
# SBATCH --ntasks-per-node=16               # total number of tasks across all nodes
#SBATCH --ntasks-per-socket=64               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1       # cpu-cores per task (>1 if multi-threaded tasks)
# SBATCH --threads-per-core=1     # setting to 1 turns off SMT (max value is 4)
# SBATCH --cores-per-socket=16     # setting to 1 turns off SMT (max value is 4)
#SBATCH --mem=4G                 # memory per cpu-core (4G is default)
#SBATCH --gres=gpu:1             # number of gpus per node
#SBATCH --gpu-bind=closest
# SBATCH --cpu-bind=verbose,cores
# SBATCH --gres-flags=enforce-binding
#SBATCH --time=00:10:00          # total run time limit (HH:MM:SS)
# SBATCH --nodelist=traverse-k01g7
# SBATCH --exclude=traverse-k01g4

hostname
echo $$
date
env | grep SLURM | sort
#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export OMP_NUM_THREADS=$SLURM_NTASKS
#export OMP_PLACES=cores
#export OMP_PROC_BIND=TRUE
#export OMP_CPU_AFFINITY=0
export GMX_MAXBACKUP=-1

module purge
#module load openmpi/devtoolset-8/4.0.1/64
module load cudatoolkit/10.2

BCH=../gpu_bench/rnase_cubic
#srun gmx_188 grompp -f $BCH/pme_verlet.mdp -c $BCH/conf.gro -p $BCH/topol.top -o bench.tpr
#srun mdrun_188 -pin on -ntomp $SLURM_CPUS_PER_TASK -s bench.tpr -c conf.gro
gmx grompp -f $BCH/pme_verlet.mdp -c $BCH/conf.gro -p $BCH/topol.top -o bench.tpr
date
numactl -s
taskset -c -p $$
date

#--gres-flags=enforce-binding
#--gres=gpu:1
export OMP_DISPLAY_ENV=TRUE
export OMP_PLACES="{11,111}"

gmx mdrun -pin on -ntomp 16 -s bench.tpr -c conf.gro
#gmx mdrun -pin on -ntomp $SLURM_CPUS_PER_TASK -s bench.tpr -c conf.gro
#gmx mdrun -pin on -ntomp $SLURM_NTASKS -s bench.tpr -c conf.gro
#srun -n 1 gmx mdrun -ntomp $SLURM_NTASKS -v -pin on -s bench.tpr
#srun -n 1 --cpu-bind=verbose,none --mem-bind=verbose,local gmx mdrun -v -pin on -s bench.tpr
#gmx mdrun -pin on -gputasks 01 -nb gpu -pme gpu -ntmpi $SLURM_NTASKS -ntomp $SLURM_CPUS_PER_TASK -s bench.tpr -c conf.gro

hostname
nvidia-smi --query-gpu=index,pci.bus_id,name --format=csv,noheader
```

The ulimit -s unlimited line makes the stack size the maximum possible. This is important if your code dynamically allocates a large amount of memory. By purging the modules you can be sure nothing has been unintentionally loaded. The module list statement is useful because it writes out the explicit module versions. This important if you later need to know exactly which modules you used. Lastly, all the SLURM environment variables are outputted. One can examine the values to see if they are as expected.

The default values for SLURM for a cluster are found here: /etc/slurm/slurm.conf

To see the run time limits for a cluster, look at: /etc/slurm/job_submit.lua

Note that `rh/devtoolset/8` cannot be used to compile Gromacs on Traverse. The IBM xlc/C compilers are not supported by Gromacs 2019.

The shared library dependencies of `gmx` are:

```
$ ldd gmx
linux-vdso64.so.1 =>  (0x0000200000050000)
libgromacs.so.4 => /home/jdh4/.local/bin/./../lib64/libgromacs.so.4 (0x0000200000070000)
libstdc++.so.6 => /lib64/libstdc++.so.6 (0x0000200001710000)
libm.so.6 => /lib64/libm.so.6 (0x00002000018a0000)
libgomp.so.1 => /lib64/libgomp.so.1 (0x0000200001990000)
libgcc_s.so.1 => /lib64/libgcc_s.so.1 (0x00002000019f0000)
libpthread.so.0 => /lib64/libpthread.so.0 (0x0000200001a30000)
libc.so.6 => /lib64/libc.so.6 (0x0000200001a70000)
libdl.so.2 => /lib64/libdl.so.2 (0x0000200001c60000)
librt.so.1 => /lib64/librt.so.1 (0x0000200001c90000)
libcufft.so.10 => /usr/local/cuda-10.1/targets/ppc64le-linux/lib/libcufft.so.10 (0x0000200001cc0000)
libhwloc.so.5 => /lib64/libhwloc.so.5 (0x000020000a050000)
libopenblas.so.0 => /lib64/libopenblas.so.0 (0x000020000a0c0000)
/lib64/ld64.so.2 (0x0000200000000000)
libnuma.so.1 => /lib64/libnuma.so.1 (0x000020000acf0000)
libltdl.so.7 => /lib64/libltdl.so.7 (0x000020000ad20000)
libgfortran.so.3 => /lib64/libgfortran.so.3 (0x000020000ad50000)
```

ESSL cannot be used:

```
-DGMX_EXTERNAL_BLAS=ON -DGMX_BLAS_USER=/usr/lib64/libessl.so \
-DGMX_EXTERNAL_LAPACK=ON -DGMX_LAPACK_USER=/usr/lib64/libessl.so \
```

Here are the undefined references:

```
[100%] Linking CXX executable ../../bin/template
../../lib/libgromacs.so.4.0.0: undefined reference to `ssteqr_'
../../lib/libgromacs.so.4.0.0: undefined reference to `dsteqr_'
../../lib/libgromacs.so.4.0.0: undefined reference to `slasr_'
../../lib/libgromacs.so.4.0.0: undefined reference to `dlartg_'
../../lib/libgromacs.so.4.0.0: undefined reference to `dlacpy_'
../../lib/libgromacs.so.4.0.0: undefined reference to `slapy2_'
../../lib/libgromacs.so.4.0.0: undefined reference to `slascl_'
../../lib/libgromacs.so.4.0.0: undefined reference to `slaset_'
../../lib/libgromacs.so.4.0.0: undefined reference to `dlapy2_'
../../lib/libgromacs.so.4.0.0: undefined reference to `dlaev2_'
../../lib/libgromacs.so.4.0.0: undefined reference to `dlarnv_'
../../lib/libgromacs.so.4.0.0: undefined reference to `dlanst_'
../../lib/libgromacs.so.4.0.0: undefined reference to `dlascl_'
../../lib/libgromacs.so.4.0.0: undefined reference to `slacpy_'
../../lib/libgromacs.so.4.0.0: undefined reference to `dorm2r_'
../../lib/libgromacs.so.4.0.0: undefined reference to `dgeqr2_'
../../lib/libgromacs.so.4.0.0: undefined reference to `dlasr_'
../../lib/libgromacs.so.4.0.0: undefined reference to `slarnv_'
../../lib/libgromacs.so.4.0.0: undefined reference to `slartg_'
../../lib/libgromacs.so.4.0.0: undefined reference to `sgeqr2_'
../../lib/libgromacs.so.4.0.0: undefined reference to `slaev2_'
../../lib/libgromacs.so.4.0.0: undefined reference to `slanst_'
../../lib/libgromacs.so.4.0.0: undefined reference to `sorm2r_'
../../lib/libgromacs.so.4.0.0: undefined reference to `dlaset_'
collect2: error: ld returned 1 exit status
make[2]: *** [share/template/CMakeFiles/template.dir/build.make:85: bin/template] Error 1
make[1]: *** [CMakeFiles/Makefile2:2311: share/template/CMakeFiles/template.dir/all] Error 2
make: *** [Makefile:163: all] Error 2
```

```
$ ls -lL /lib64/*essl*
-rw-r--r--. 1 bin bin 45719787 Mar 29  2018 /lib64/libessl6464.so
-rw-r--r--. 1 bin bin 45719787 Mar 29  2018 /lib64/libessl6464.so.1
-rw-r--r--. 1 bin bin 45719787 Mar 29  2018 /lib64/libessl6464.so.1.10
-rw-r--r--. 1 bin bin 53379191 Mar 29  2018 /lib64/libesslsmp6464.so
-rw-r--r--. 1 bin bin 53379191 Mar 29  2018 /lib64/libesslsmp6464.so.1
-rw-r--r--. 1 bin bin 53379191 Mar 29  2018 /lib64/libesslsmp6464.so.1.10
-rw-r--r--. 1 bin bin 54737430 Mar 29  2018 /lib64/libesslsmpcuda.so
-rw-r--r--. 1 bin bin 54737430 Mar 29  2018 /lib64/libesslsmpcuda.so.1
-rw-r--r--. 1 bin bin 54737430 Mar 29  2018 /lib64/libesslsmpcuda.so.1.10
-rw-r--r--. 1 bin bin 53925425 Mar 29  2018 /lib64/libesslsmp.so
-rw-r--r--. 1 bin bin 53925425 Mar 29  2018 /lib64/libesslsmp.so.1
-rw-r--r--. 1 bin bin 53925425 Mar 29  2018 /lib64/libesslsmp.so.1.10
-rw-r--r--. 1 bin bin 46826939 Mar 29  2018 /lib64/libessl.so
-rw-r--r--. 1 bin bin 46826939 Mar 29  2018 /lib64/libessl.so.1
-rw-r--r--. 1 bin bin 46826939 Mar 29  2018 /lib64/libessl.so.1.10
```

## Della

Della is good for single node jobs. You should not be running small jobs on Tiger.

```bash
$ ssh <YourNetID>@della.princeton.edu
$ cd </path/to/your/software/directory>  # e.g., cd ~/software
$ wget https://raw.githubusercontent.com/PrincetonUniversity/running_gromacs/master/02_installation/della/della.sh
# make modifications to della.sh if needed (e.g., choose a different version)
$ bash della.sh | tee build.log
```

For single-node jobs:

```bash
#!/bin/bash
#SBATCH --job-name=gmx           # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=8               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G per cpu-core is default)
#SBATCH --time=01:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all          # send email when job begins, ends and fails
#SBATCH --mail-user=<YourNetID>@princeton.edu
#SBATCH --constraint=haswell|broadwell  # exclude ivy nodes

module purge
module load intel/19.0/64/19.0.5.281

gmx grompp -f pme_verlet.mdp -c conf.gro -p topol.top -o bench.tpr
gmx mdrun -ntmpi $SLURM_NTASKS -ntomp $SLURM_CPUS_PER_TASK -s bench.tpr
```

For multi-node MPI jobs:

```bash
#!/bin/bash
#SBATCH --job-name=gmx           # create a short name for your job
#SBATCH --nodes=2                # node count
#SBATCH --ntasks-per-node=32     # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G per cpu-core is default)
#SBATCH --time=01:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all          # send email when job begins, ends and fails
#SBATCH --mail-user=<YourNetID>@princeton.edu
#SBATCH --constraint=haswell|broadwell  # exclude ivy nodes

module purge
module load intel/19.0/64/19.0.5.281
module load intel-mpi/intel/2018.3/64

gmx grompp -f pme_verlet.mdp -c conf.gro -p topol.top -o bench.tpr
srun mdrun_mpi -ntomp $SLURM_CPUS_PER_TASK -s bench.tpr
```

## Perseus

Perseus is a lage cluster composed of Intel Broadwell CPUs. It is mainly used for astrophysics research. Run these commands to build GROMACS on perseus:

```bash
$ ssh <YourNetID>@perseus.princeton.edu
$ cd </path/to/your/software/directory>  # e.g., cd ~/software
$ wget https://raw.githubusercontent.com/PrincetonUniversity/running_gromacs/master/02_installation/perseus/perseus.sh
# make modifications to perseus.sh if needed (e.g., choose a different version)
$ bash perseus.sh | tee build.log
```

## TigerCPU

If you have an account on Tiger then consider building only a GPU version:

```
$ ssh <NetID>@tigercpu.princeton.edu
$ cd </path/to/software/directory>  # e.g., cd ~/software
$ mkdir cpu_version && cd cpu_version
$ wget https://raw.githubusercontent.com/PrincetonUniversity/running_gromacs/master/02_installation/tigerCpu/tigerCpu.sh
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
module load intel/19.0/64/19.0.5.281

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
module load intel/19.0/64/19.0.5.281
module load intel-mpi/intel/2019.5/64

gmx_cpu grompp -f $BCH/pme_verlet.mdp -c $BCH/conf.gro -p $BCH/topol.top -o bench.tpr
srun mdrun_cpu_mpi -ntomp $SLURM_CPUS_PER_TASK -s bench.tpr
```

## Adroit (GPU)

Gromacs can be built for GPU use on Adroit. Because it has four V100 GPUs, one can run comparison tests against Traverse. Gromacs can be built by running the following commands:

```bash
$ ssh <NetID>@adroit.princeton.edu
$ cd </path/to/your/software/directory>  # e.g., cd ~/software
$ wget https://raw.githubusercontent.com/PrincetonUniversity/running_gromacs/main/02_installation/adroit/adroit.sh
# make modifications to adroit.sh if needed (e.g., choose a different version)
$ bash adroit.sh | tee build.log
```

Below is a sample Slurm script for single-node jobs:

```
#!/bin/bash
#SBATCH --job-name=rnase         # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=4G                 # memory per node (4G per cpu-core is default)
#SBATCH --gres=gpu:tesla_v100:1  # number of gpus per node
#SBATCH --time=00:10:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all          # send email when job begins, ends and fails
#SBATCH --mail-user=<YourNetID>@princeton.edu

module purge
module load intel/19.0/64/19.0.1.144
module load cudatoolkit/10.2

gmx grompp -f pme_verlet.mdp -c conf.gro -p topol.top -o bench.tpr
gmx mdrun -ntmpi $SLURM_NTASKS -ntomp $SLURM_CPUS_PER_TASK -s bench.tpr
```

One should find the optimal values of `ntasks` and `cpus-per-task`. For a job that uses a single GPU, the product of those two parameters should not exceed 10.

## Adroit (CPU)

Adroit is a small cluster so the directions here produce a single-node version:

```bash
$ ssh <NetID>@adroit.princeton.edu
$ cd </path/to/your/software/directory>  # e.g., cd ~/software
$ mkdir cpu_version && cd cpu_version
$ wget https://raw.githubusercontent.com/PrincetonUniversity/running_gromacs/main/02_installation/adroit/adroit_cpu.sh
# make modifications to adroit_cpu.sh if needed (e.g., choose a different version)
$ bash adroit_cpu.sh | tee build.log
```

Below is a sample Slurm script for single-node jobs:

```
#!/bin/bash
#SBATCH --job-name=rnase         # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=4        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=4G                 # memory per node (4G per cpu-core is default)
#SBATCH --time=00:10:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all          # send email when job begins, ends and fails
#SBATCH --mail-user=<YourNetID>@princeton.edu

module purge
module load intel/19.0/64/19.0.1.144

gmx_cpu grompp -f pme_verlet.mdp -c conf.gro -p topol.top -o bench.tpr
gmx_cpu mdrun -ntmpi $SLURM_NTASKS -ntomp $SLURM_CPUS_PER_TASK -s bench.tpr
```

One should find the optimal values of `ntasks` and `cpus-per-task`.
