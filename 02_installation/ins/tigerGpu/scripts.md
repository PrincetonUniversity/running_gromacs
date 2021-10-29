# TigerGPU

The code can be built by following these commands:

```
$ ssh <YourNetID>@tigergpu.princeton.edu
$ cd </path/to/your/software/directory>  # e.g., cd ~/software
$ wget https://raw.githubusercontent.com/PrincetonUniversity/running_gromacs/master/02_installation/ins/tigerGpu/tigerGpu.sh
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
