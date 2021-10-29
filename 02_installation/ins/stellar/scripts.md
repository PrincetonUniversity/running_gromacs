# Stellar-Intel

```
$ ssh <YourNetID>@stellar-intel.princeton.edu
$ cd </path/to/your/software/directory>  # e.g., cd ~/software
$ https://raw.githubusercontent.com/PrincetonUniversity/running_gromacs/main/02_installation/ins/stellar/stellar_intel_gmx_gcc.sh
# make modifications to stellar_intel_gmx_gcc.sh if needed (e.g., choose a different version)
$ bash stellar_intel_gmx_gcc.sh | tee build.log
```

Below is a sample Slurm script:

```
#!/bin/bash
#SBATCH --job-name=gmx           # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=96              # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:10:00          # total run time limit (HH:MM:SS)

module purge
module load openmpi/gcc/4.1.0

BCH=rnase_cubic
gmx_mpi grompp -f $BCH/pme_verlet.mdp -c $BCH/conf.gro -p $BCH/topol.top -o bench.tpr
srun gmx_mpi mdrun -ntomp $SLURM_CPUS_PER_TASK -s bench.tpr
```
