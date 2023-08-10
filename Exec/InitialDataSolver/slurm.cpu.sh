#!/bin/bash

## specify your allocation (with the _g) and that you want GPU nodes
#SBATCH -A mp111
#SBATCH -C cpu

## the job will be named "Z4c2" in the queue and will save stdout to z4c2_[job ID].out
#SBATCH -J Z4c2
#SBATCH -o z4c2_%j.out

#SBATCH --qos=regular

## set the max walltime
#SBATCH -t 03:00:00

## specify the number of nodes you want
#SBATCH -N 16

## we use the same number of MPI ranks per node as GPUs per node
#SBATCH --ntasks-per-node=64

export OMP_PLACES=threads
export OMP_PROC_BIND=true

export OMP_NUM_THREADS=8

# the -n argument is (--ntasks-per-node) * (-N) = (number of MPI ranks per node) * (number of nodes)
srun -n 1024 --cpu_bind=cores -c 4 ./main3d.gnu.x86-milan.MPI.OMP.ex inputs
