## Adroit (GPU)

Gromacs can be built for the A100 GPU node of Adroit:

```bash
$ ssh <NetID>@adroit.princeton.edu
$ cd </path/to/your/software/directory>  # e.g., cd ~/software
$ wget https://raw.githubusercontent.com/PrincetonUniversity/running_gromacs/main/02_installation/ins/adroit/adroit.sh
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
#SBATCH --gres=gpu:1             # number of gpus per node
#SBATCH --time=00:10:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all          # send email when job begins, ends and fails
#SBATCH --mail-user=<YourNetID>@princeton.edu
#SBATCH --constraint=a100

module purge
module load cudatoolkit/11.4

gmx grompp -f pme_verlet.mdp -c conf.gro -p topol.top -o bench.tpr
gmx mdrun -ntmpi $SLURM_NTASKS -ntomp $SLURM_CPUS_PER_TASK -s bench.tpr
```

One should find the optimal values of `ntasks` and `cpus-per-task`. For a job that uses a single GPU, the product of those two parameters should not exceed 10.
