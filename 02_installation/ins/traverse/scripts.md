## Traverse

Traverse is composed of 46 IBM POWER9 nodes with 4 NVIDIA V100 GPUs per node. Each node has two sockets each with a 16-core CPU. Each of the 16 cores has 4 hardware threads. Hence each node has 128 logical cpu-cores or hardware threads.

Run the commands below to install Gromacs:

```
$ ssh <YourNetID>@traverse.princeton.edu
$ cd </path/to/your/software/directory>  # e.g., cd ~/software
$ wget https://raw.githubusercontent.com/PrincetonUniversity/running_gromacs/main/02_installation/ins/traverse/traverse_rhel8_v2021.sh
$ bash traverse_rhel8_2021.sh | tee build.log
```

Below is a sample Slurm script:

```bash

#!/bin/bash
#SBATCH --job-name=gmx           # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=8               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=32G                # memory per node (4G per cpu-core is default)
#SBATCH --time=01:00:00          # total run time limit (HH:MM:SS)
#SBATCH --gres=gpu:1             # number of gpus per node
#SBATCH --mail-type=all          # send email when job begins, ends and fails
#SBATCH --mail-user=<YourNetID>@princeton.edu

module purge
module load cudatoolkit/11.4
module load openmpi/gcc/4.1.1/64

BCH=./ADH/adh_cubic
gmx_mpi grompp -f $BCH/pme_verlet.mdp -c $BCH/conf.gro -p $BCH/topol.top -o bench.tpr
srun gmx_mpi mdrun -update gpu -pin on -ntomp $SLURM_CPUS_PER_TASK -s bench.tpr
```

## Poor Performance

We have seen poor performance from Gromacs on Traverse as have others on the POWER9 systems. You are encouraged to take entire nodes by requesting 4 GPUs. This prevents job sharing which helps with the performance issues.

## Other

The content below is outdated and is only left here to provides ideas for current troubleshooting.

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
