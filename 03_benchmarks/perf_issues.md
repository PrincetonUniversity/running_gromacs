# Performance Issues

`srun -t 1 -l -n 4 -N 1 numactl -s | grep physc | sort`

```
#!/bin/bash
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=16       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --threads-per-core=1     # setting to 1 turns off SMT (max value is 4)
#SBATCH --mem=4G                 # memory per cpu-core (4G is default)
#SBATCH --gres=gpu:1             # number of gpus per node
#SBATCH --gpu-bind=closest
#SBATCH --time=00:10:00          # total run time limit (HH:MM:SS)
# SBATCH --nodelist=traverse-k02g4
# SBATCH --exclude=traverse-k01g2

hostname

# OMP_PLACES
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export GMX_MAXBACKUP=-1

module purge
module load openmpi/devtoolset-8/4.0.1/64
module load cudatoolkit/10.2

BCH=../gpu_bench/rnase_cubic
srun gmx grompp -f $BCH/pme_verlet.mdp -c $BCH/conf.gro -p $BCH/topol.top -o bench.tpr

date
taskset -c -p $$
date

srun gmx mdrun -pin on -ntomp $SLURM_CPUS_PER_TASK -s bench.tpr -c conf.gro
```

## When no other jobs running on the node

```
$ cat slurm-37863.out
traverse-k01g2
...
Wed Dec 11 23:17:01 EST 2019
pid 31313's current affinity list: 0-63
Wed Dec 11 23:17:01 EST 2019
...
Overriding thread affinity set outside gmx mdrun
starting mdrun 'RNASE ZF-1A in water'
10000 steps,     20.0 ps.

Writing final coordinates.

               Core t (s)   Wall t (s)        (%)
       Time:       61.892        3.870     1599.4
                 (ns/day)    (hour/ns)
Performance:      446.590        0.054
```

## When another job running on the node (GPU bus id = 2 via nvidia-smi)

```
$ cat slurm-37864.out
traverse-k02g4
...
Wed Dec 11 23:19:01 EST 2019
pid 22302's current affinity list: 96-111
Wed Dec 11 23:19:01 EST 2019
...
Overriding thread affinity set outside gmx mdrun
NOTE: Affinity setting for 14/16 threads failed.

NOTE: Thread affinity was not set.
starting mdrun 'RNASE ZF-1A in water'
10000 steps,     20.0 ps.

Writing final coordinates.

               Core t (s)   Wall t (s)        (%)
       Time:      150.780        9.424     1599.9
                 (ns/day)    (hour/ns)
Performance:      183.377        0.131
```
