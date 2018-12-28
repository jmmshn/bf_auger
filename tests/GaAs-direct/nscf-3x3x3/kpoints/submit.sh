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

for((i=1;i<=27;i++))
do
cd k_$i
srun -n 68 -c4 --cpu_bind=cores ~/espresso/QE-6.2/qe-6.2/bin/pw.x < nscf.in > nscf.out
srun -n 68 -c4 --cpu_bind=cores ~/espresso/QE-6.2/qe-6.2/bin/pw2wannier90.x < pw2wannier.in > pw2wannier.out
rm GaAs.wfc*
cd ..
done
