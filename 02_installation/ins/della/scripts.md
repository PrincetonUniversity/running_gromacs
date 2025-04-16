
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

# Della (CPU for OS Version 9)

A CPU-only build for the AMD nodes:

```bash
$ ssh <YourNetID>@della.princeton.edu
$ cd </path/to/your/software/directory>  # e.g., cd ~/software
$ wget https://raw.githubusercontent.com/PrincetonUniversity/running_gromacs/main/02_installation/ins/della/della9_cpu_amd.sh
# make modifications to della_cpu_amd.sh if needed (e.g., choose a different version)
$ bash della_cpu_amd.sh | tee build.log
```

Below is a sample Slurm script:

```bash
#!/bin/bash
#SBATCH --job-name=gmx           # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=8               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G per cpu-core is default)
#SBATCH --time=01:00:00          # total run time limit (HH:MM:SS)
#SBATCH --constraint=amd
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=<YourNetID>@princeton.edu

module purge
module load gcc-toolset/14
module load aocc/5.0.0
module load aocl/aocc/5.0.0
module load openmpi/aocc-5.0.0/4.1.6

BCH=rnase_cubic
gmx_amd_cpu grompp -f $BCH/pme_verlet.mdp -c $BCH/conf.gro -p $BCH/topol.top -o bench.tpr
srun gmx_amd_cpu mdrun -ntomp $SLURM_CPUS_PER_TASK -s bench.tpr
```

There is also a [GCC build](della9_cpu_gcc.sh) but the above showed a performance advantage.

The Della AMD nodes provide 192 CPU-cores. You may consider building GMX without MPI and instead using the built-in thread MPI if single-node runs are sufficient for your work.
