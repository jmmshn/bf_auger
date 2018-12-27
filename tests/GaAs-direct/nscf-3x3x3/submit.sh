#!/bin/bash -l

#SBATCH -J GaAs
#SBATCH -p debug
#SBATCH -C knl,quad,cache
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -A m934

export OMP_PROC_BIND=true
export OMP_PLACES=threads

export OMP_NUM_THREADS=4
export OMP_STACKSIZE=256m
export MKL_FAST_MEMORY_LIMIT=0

module load craype-hugepages2M

srun -n 68 -c4 --cpu_bind=cores ~/bin/Auger/src-v2/auger_eeh.x < auger.in > auger_eeh.out
srun -n 68 -c4 --cpu_bind=cores ~/bin/Auger/src-v2/auger_hhe.x < auger.in > auger_hhe.out
