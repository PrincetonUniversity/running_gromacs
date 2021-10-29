# Della-GPU

The GPU nodes on Della have AMD processors with AVX2 as the highest instruction set. The system version of GCC is 8.3.1. Below is a build procedure:

```
$ ssh <YourNetID>@della-gpu.princeton.edu
$ cd software  # or another directory
$ wget https://raw.githubusercontent.com/PrincetonUniversity/running_gromacs/main/02_installation/ins/della/della_gpu.sh
$ bash della_gpu.sh | tee gmx.log
```

The above procedure will create `gmx_gpu`.

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
