# TigerCPU

If you have an account on Tiger then consider building only a GPU version:

```
$ ssh <NetID>@tigercpu.princeton.edu
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
