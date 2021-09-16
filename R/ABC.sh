#!/bin/sh

#SBATCH -t 168:00:00

#SBATCH -J large_free

#SBATCH -N 1

#SBATCH --tasks-per-node=19

#SBATCH --exclusive

#SBATCH --mem-per-cpu=3100

#SBATCH --mail-user=olivia@kth.se

#SBATCH --mail-type=BEGIN,END

module load GCC/8.3.0  CUDA/10.1.243  OpenMPI/3.1.4
module load GCC/8.3.0  OpenMPI/3.1.4
module load R/3.6.2

R CMD BATCH --vanilla runABCMCMC_Matt_LU_aut.R ABC_output.Rout


