
# Della-GPU (OS Version 9)

The GPU nodes on Della have range of instruction sets. Here we choose AVX2 for compatibility. Below is a build procedure:

```
$ ssh <YourNetID>@della9.princeton.edu
$ cd software  # or another directory
$ wget https://raw.githubusercontent.com/PrincetonUniversity/running_gromacs/main/02_installation/ins/della/della9_gpu_avx2.sh
$ bash della9_gpu_avx2.sh | tee gmx.log
```

The above produces an executable called `gmx_d9_gpu`.

For additional performance, one might try with cuFFTMp, using AVX-512 and then only running on the Intel CPUs with 80 GB A100's, and using openmpi/cuda-12.6/gcc/4.1.6 for Open MPI.

# Della-GPU (OS Version 8)

The GPU nodes on Della have AMD processors with AVX2 as the highest instruction set. The system version of GCC is 8.3.1. Below is a build procedure:

```
$ ssh <YourNetID>@della-gpu.princeton.edu
$ cd software  # or another directory
$ wget https://raw.githubusercontent.com/PrincetonUniversity/running_gromacs/main/02_installation/ins/della/della_gpu.sh
$ bash della_gpu.sh | tee gmx.log
```

The above procedure will create `gmx_gpu`. Below is a sample Slurm script:

```bash
#!/bin/bash
#SBATCH --job-name=gmx           # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=8        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=32G                # memory per node (4G per cpu-core is default)
#SBATCH --time=01:00:00          # total run time limit (HH:MM:SS)
#SBATCH --gres=gpu:1             # number of gpus per node
#SBATCH --mail-type=all          # send email when job begins, ends and fails
#SBATCH --mail-user=<YourNetID>@princeton.edu

module purge
module load cudatoolkit/11.4
module load openmpi/gcc/4.1.0

BCH=./ADH/adh_cubic
gmx_gpu grompp -f $BCH/pme_verlet.mdp -c $BCH/conf.gro -p $BCH/topol.top -o bench.tpr
srun gmx_gpu mdrun -update gpu -pin on -ntomp $SLURM_CPUS_PER_TASK -s bench.tpr
```

#### NGC Container

As an alternative to the procedure above, you may consider using the [NVIDIA container](https://ngc.nvidia.com/catalog/containers/hpc:gromacs):

```
$ singularity pull docker://nvcr.io/hpc/gromacs:2021
```

The container appears to be restricted to a single node. So you can you use up to 2 GPUs and 128 CPU-cores. You will need to conduct a scaling analysis to find the optimal numbers. However, one GPU and a few CPU-cores should give excellent performance for most systems.

```
#!/bin/bash
#SBATCH --job-name=gmx           # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=12       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=4G                 # memory per node (4G per cpu-core is default)
#SBATCH --time=00:10:00          # total run time limit (HH:MM:SS)
#SBATCH --gres=gpu:1             # number of gpus per node

module purge
SIF=$(pwd)/gromacs_2021.sif
SINGULARITY="singularity run --nv -B ${PWD}:/host_pwd --pwd /host_pwd ${SIF}"

BCH=./rnase_cubic
${SINGULARITY} gmx grompp -f $BCH/pme_verlet.mdp -c $BCH/conf.gro -p $BCH/topol.top -o bench.tpr
${SINGULARITY} gmx mdrun -ntmpi $SLURM_NTASKS -ntomp $SLURM_CPUS_PER_TASK -update gpu -s bench.tpr
```

`md.log` shows that the GPU was correctly found:

```
GPU info:
    Number of GPUs detected: 1
    #0: NVIDIA NVIDIA A100-PCIE-40GB, compute cap.: 8.0, ECC: yes, stat: compatible
```

# Della (CPU)

Della is good for single node jobs. You should not be running small jobs on Tiger.

```bash
$ ssh <YourNetID>@della.princeton.edu
$ cd </path/to/your/software/directory>  # e.g., cd ~/software
$ wget https://raw.githubusercontent.com/PrincetonUniversity/running_gromacs/main/02_installation/ins/della/della.sh
# make modifications to della.sh if needed (e.g., choose a different version)
$ bash della.sh | tee build.log
```

Below is a sample Slurm script:

```bash
#!/bin/bash
#SBATCH --job-name=gmx           # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=4               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G per cpu-core is default)
#SBATCH --time=01:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=<YourNetID>@princeton.edu
#SBATCH --exclude=della-r4c[1-4]n[1-16],della-r1c[3,4]n[1-16]

module purge
module load intel-mkl/2022.2.0
module load openmpi/gcc/4.1.2

BCH=rnase_cubic
gmx_mpi grompp -f $BCH/pme_verlet.mdp -c $BCH/conf.gro -p $BCH/topol.top -o bench.tpr
srun gmx_mpi mdrun -ntomp $SLURM_CPUS_PER_TASK -s bench.tpr
```
